/**
 * PhLynx half of the module pipeline test: builds a model from this library's modules with
 * PhLynx's *real* code (a checkout at PHLYNX_DIR) and exports the .omex PhLynx would send to
 * CUFLynx.
 *
 *   node tools/phlynx_bridge/export_omex.mjs <job.json> <out.omex>
 *
 * job.json:
 *   { "cellml": [paths], "units": [paths], "configs": [paths],       library files to load
 *     "instances": [{name, module_type, module_subtype, inp_instances, out_instances}],
 *     "parameters": [{variable_name, units, value, data_reference}],
 *     "sim_time": 2.0, "dt": 0.01 }
 *
 * Prints a JSON report on stdout: nodes and edges built, the connections PhLynx refused, the
 * modules it could not find, warnings, and the export error if there was one. Exits 0 even when
 * PhLynx fails, so the Python test can report why; exits 3 when the bridge itself can't run.
 *
 * Everything is PhLynx's own code. Parameters go through PhLynx's applyParametersToNodes
 * (src/utils/parameters.js); for an older checkout without it, applyParameters below mirrors
 * loadParametersData from src/views/WorkspaceArea.vue.
 *
 * DOM: jsdom, NOT happy-dom: happy-dom's selector engine doesn't match tag names with an
 * underscore (map_variables), which silently empties PhLynx's CellML connection parsing
 * (see CUFLynx apps/api/tests/phlynx_bridge/roundtrip.mjs).
 * Needs node >= 22.15 (module.registerHooks) and jsdom, installed next to this script
 * (npm install in tools/phlynx_bridge) or found from JSDOM_DIR.
 */
import { readFile, writeFile } from 'node:fs/promises'
import { createRequire, registerHooks } from 'node:module'
import path from 'node:path'
import { pathToFileURL } from 'node:url'

const PHLYNX_DIR = process.env.PHLYNX_DIR
const JSDOM_DIR = process.env.JSDOM_DIR || path.dirname(new URL(import.meta.url).pathname)
const [, , jobPath, outPath] = process.argv

if (!PHLYNX_DIR || !jobPath || !outPath) {
  console.error('BRIDGE_UNSUPPORTED: usage PHLYNX_DIR=... node export_omex.mjs <job.json> <out.omex>')
  process.exit(3)
}
if (typeof registerHooks !== 'function') {
  console.error('BRIDGE_UNSUPPORTED: node is too old for module.registerHooks (need >= 22.15)')
  process.exit(3)
}

const fromPhlynx = createRequire(`${PHLYNX_DIR}/package.json`)
const fromJsdom = createRequire(`${JSDOM_DIR}/package.json`)

// Load a package by its ESM entry, the same instance PhLynx's own imports get (a dual-published
// package's CommonJS entry is a different instance: pinia's setActivePinia would not be the one
// PhLynx's stores see).
async function load(req, name) {
  const pkgPath = req.resolve(`${name}/package.json`)
  const pkg = JSON.parse(await readFile(pkgPath, 'utf8'))
  const entry = pkg.exports?.['.']?.import?.default ?? pkg.exports?.['.']?.import ?? pkg.module ?? pkg.main
  const base = pathToFileURL(pkgPath)
  return import(typeof entry === 'string' ? new URL(entry, base).href : pathToFileURL(req.resolve(name)).href)
}

// PhLynx is a Vite app with extensionless relative imports; add '.js' when Node can't resolve.
registerHooks({
  resolve(specifier, context, nextResolve) {
    try {
      return nextResolve(specifier, context)
    } catch (err) {
      const looksLikePath = specifier.startsWith('.') || specifier.startsWith('file:')
      if (looksLikePath && !/\.[mc]?js$/.test(specifier)) return nextResolve(`${specifier}.js`, context)
      throw err
    }
  },
})

const warnings = []
const origWarn = console.warn
console.warn = (...args) => {
  warnings.push(args.map(String).join(' '))
}
const origLog = console.log
const logs = []
console.log = (...args) => {
  logs.push(args.map(String).join(' '))
}

const report = { ok: false, stage: 'setup' }
try {
  const { JSDOM } = await load(fromJsdom, 'jsdom')
  const dom = new JSDOM('<!doctype html><html></html>')
  globalThis.DOMParser = dom.window.DOMParser
  globalThis.Node = dom.window.Node
  globalThis.XMLSerializer = dom.window.XMLSerializer
  globalThis.document = dom.window.document
  globalThis.__APP_VERSION__ = JSON.parse(await readFile(`${PHLYNX_DIR}/package.json`, 'utf8')).version

  const { setActivePinia, createPinia } = await load(fromPhlynx, 'pinia')
  setActivePinia(createPinia())
  const createLibCellML = (await load(fromPhlynx, 'libcellml.js')).default
  const wasm = path.join(path.dirname(fromPhlynx.resolve('libcellml.js/package.json')), 'libcellml.wasm')
  const libcellml = await createLibCellML({ locateFile: (f, prefix) => (f.endsWith('.wasm') ? wasm : prefix + f) })

  const cellml = await import(`${PHLYNX_DIR}/src/utils/cellml.js`)
  cellml.initLibCellML(libcellml)
  const { useLibraryStore } = await import(`${PHLYNX_DIR}/src/stores/libraryStore.js`)
  const { buildWorkflowGraph } = await import(`${PHLYNX_DIR}/src/services/import/buildWorkflow.js`)
  const { generateOmexArchive } = await import(`${PHLYNX_DIR}/src/services/compress.js`)
  const JSZip = (await load(fromPhlynx, 'jszip')).default

  const job = JSON.parse(await readFile(jobPath, 'utf8'))
  const library = useLibraryStore()

  // --- load the library, as hydrateCellmlAndDependents / loadCellMLData do ---
  report.stage = 'load library'
  report.load_errors = []
  for (const file of [...(job.cellml || []), ...(job.units || [])]) {
    const result = cellml.processCellMLData(await readFile(file, 'utf8'))
    const name = path.basename(file)
    if (result.type !== 'success') {
      report.load_errors.push(`${name}: ${result.message ?? JSON.stringify(result).slice(0, 300)}`)
      continue
    }
    if ((result.components?.length ?? 0) > 0) library.addMathFile(name, result.components)
    if (result.units?.count > 0) library.addUnitsFile({ componentFile: name, model: result.units.model })
  }
  report.configs_loaded = 0
  for (const file of job.configs || []) {
    report.configs_loaded += library.addConfigFile(path.basename(file), JSON.parse(await readFile(file, 'utf8')))
  }

  // --- build the graph from the instance array, as useLoadFromInstanceArray does ---
  report.stage = 'build graph'
  const wanted = job.instances.map((r) => `${r.module_type}:${r.module_subtype}`)
  report.missing_modules = [...new Set(wanted.filter((ref) => !library.availableModules.has(ref)))]
  report.stub_modules = [...new Set(wanted.filter((ref) => library.availableModules.get(ref)?.isStub))]
  const { pendingInstances: nodes, pendingEdges: edges } = buildWorkflowGraph(
    job.instances, library.availableModules, [])
  report.nodes = nodes.length
  report.edges = edges.length
  report.expected_edges = job.instances.reduce((n, r) => n + String(r.out_instances || '').split(/\s+/).filter(Boolean).length, 0)
  report.refused_edges = warnings.filter((w) => w.startsWith('[buildEdges]'))
  // multi_port warnings are printed while the whole library loads; keep the ones about modules
  // this model uses
  report.multiport_warnings = warnings.filter((w) => /multi.?port/i.test(w) && !w.startsWith('[buildEdges]')
    && wanted.some((ref) => w.includes(`"${ref}"`)))

  // --- parameters: PhLynx's own applyParametersToNodes (src/utils/parameters.js) where it exists,
  // else a mirror of loadParametersData for older PhLynx checkouts ---
  report.stage = 'apply parameters'
  let phlynxParameters = null
  try {
    phlynxParameters = await import(`${PHLYNX_DIR}/src/utils/parameters.js`)
  } catch {
    phlynxParameters = null
  }
  if (phlynxParameters?.applyParametersToNodes) {
    report.parameters_applied = phlynxParameters.applyParametersToNodes(nodes, job.parameters || [], library).totalUpdated
    report.parameters_via = 'PhLynx applyParametersToNodes'
  } else {
    report.parameters_applied = applyParameters(nodes, job.parameters || [], library)
    report.parameters_via = 'bridge mirror of loadParametersData'
  }

  // --- export, as generateOmexArchiveAction (src/composables/useImportExportSend.js) does ---
  report.stage = 'flatten'
  const blob = await cellml.generateFlattenedModel(nodes, edges, library, [])
  const modelText = await blob.text()
  report.stage = 'extract parameters'
  const extractedData = cellml.extractVoiAndParametersFromModel(modelText, {})
  report.stage = 'omex'
  const simTime = Number(job.sim_time ?? 1)
  const dt = Number(job.dt ?? 0.01)
  // The same CellML as text: Node's Blob isn't a type JSZip reads (CUFLynx's bridge does the same).
  const returned = await generateOmexArchive(
    { blob: modelText },
    JSON.stringify({ nodeData: nodes, edgeData: edges, mathLibrary: {} }),
    { simulationSettings: { initialPoint: 0, startingPoint: 0, endingPoint: simTime, pointInterval: dt },
      plotConfig: [], parameterScanConfig: {} },
    { extractedData, modified: true, cellmlFileName: 'model.cellml' },
  )
  const out = returned?.blob ?? returned
  const bytes = out instanceof Uint8Array ? out : new Uint8Array(await out.arrayBuffer())
  await writeFile(outPath, bytes)
  report.omex_members = Object.keys((await JSZip.loadAsync(bytes)).files)
  report.flattened_bytes = modelText.length
  report.stage = 'done'
  report.ok = true
} catch (err) {
  report.error = String(err?.stack ?? err).slice(0, 4000)
} finally {
  report.warnings = warnings.slice(0, 200)
  report.log_tail = logs.slice(-40)
  console.warn = origWarn
  console.log = origLog
  process.stdout.write(JSON.stringify(report) + '\n')
}

function applyParameters(nodes, content, library) {
  // Mirrors loadParametersData: a node's variable <var> takes the row named <var>_<instance>;
  // rows naming a variable some node has, with no instance suffix, become global constants.
  const catalogue = new Set()
  let total = 0
  for (const node of nodes) {
    const instance = node.data.name
    for (const v of node.data.variables) catalogue.add(v.name.trim())
    const rows = content
      .filter((e) => e.variable_name.trimEnd().endsWith(instance))
      .map((e) => ({ ...e, name: e.variable_name.trimEnd().slice(0, -instance.length).replace(/_+$/, '') }))
    const byName = new Map(rows.map((p) => [p.name.trim(), p]))
    node.data.variables = node.data.variables.map((v) => {
      const m = byName.get(v.name.trim())
      if (!m) return v
      total++
      if (m.units.trim() !== v.units.trim()) {
        console.warn(`Unit mismatch for "${v.name}": node has "${v.units.trim()}", parameter has "${m.units.trim()}"`)
      }
      return { ...v, value: String(m.value).trim(), data_reference: String(m.data_reference ?? '').trim(), type: 'constant' }
    })
  }
  for (const e of content.filter((e) => catalogue.has(e.variable_name.trimEnd()))) {
    library.assignGlobalConstant(e.variable_name.trimEnd(), String(e.value), e.units, e.data_reference)
    total++
  }
  return total
}
