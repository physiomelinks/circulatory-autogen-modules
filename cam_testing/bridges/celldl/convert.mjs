/**
 * Makes CellML from a CellDL bond graph with PhLynx's converter, for cam_testing's bondgraph_cellml_test:
 *
 *   node cam_testing/bridges/celldl/convert.mjs <document.svg> <out.cellml>
 *
 * The converter is @physiomelinks/bg2cellml (the celldl2cellml repo), through convert(), with the BG-RDF framework the
 * CellDL editor bundles (@celldl/editor/bgrdf), at the versions PhLynx pins (package.json). A document whose
 * transformers end at a voltage source, or that sets its own RT (bgf:constantValue), needs converter 0.6.3: set
 * CELLDL_BRIDGE_DIR to a build of it until it is released. To convert with an unreleased BG-RDF, set BGRDF_DIR to a
 * BG-RDF checkout: its schema/ontology.ttl replaces the bundled ontology, and each templates/<name>.ttl there replaces
 * the bundled template of that name or adds it (bundled templates it doesn't have are kept). Prints a JSON report on stdout: { ok, issues, warnings,
 * generator }, where issues are the converter's errors (ok is false with any) and warnings don't stop it. Exits 0 when the converter ran (ok says whether it made CellML), 3 when the bridge
 * can't run (its packages aren't installed: npm ci here, or set CELLDL_BRIDGE_DIR to a directory where they are).
 */
import { readdir, readFile, writeFile } from 'node:fs/promises'
import { createRequire } from 'node:module'
import path from 'node:path'
import { pathToFileURL } from 'node:url'

const BRIDGE_DIR = process.env.CELLDL_BRIDGE_DIR || path.dirname(new URL(import.meta.url).pathname)
const [, , documentPath, outPath] = process.argv
if (!documentPath || !outPath) {
  console.error('BRIDGE_UNSUPPORTED: usage node convert.mjs <document.svg> <out.cellml>')
  process.exit(3)
}

const req = createRequire(path.join(BRIDGE_DIR, 'package.json'))

// A package by its ESM entry (or one of its exports), resolved from the bridge's node_modules.
async function load(name, subpath = '.') {
  // a package whose exports hide its package.json is found in node_modules directly
  let pkgPath = path.join(BRIDGE_DIR, 'node_modules', name, 'package.json')
  try {
    pkgPath = req.resolve(`${name}/package.json`)
  } catch {
    try {
      await readFile(pkgPath)
    } catch {
      console.error(`BRIDGE_UNSUPPORTED: ${name} is not installed in ${BRIDGE_DIR} (npm ci there, or set CELLDL_BRIDGE_DIR)`)
      process.exit(3)
    }
  }
  const pkg = JSON.parse(await readFile(pkgPath, 'utf8'))
  const exp = typeof pkg.exports === 'object' ? pkg.exports[subpath] : undefined
  const entry = exp?.import?.default ?? exp?.import ?? exp?.default ?? (subpath === '.' ? pkg.module ?? pkg.main : null)
  return import(new URL(typeof entry === 'string' ? entry : 'index.js', pathToFileURL(pkgPath)).href)
}

const { DOMParser, XMLSerializer } = await load('@xmldom/xmldom')
globalThis.DOMParser ??= DOMParser
globalThis.XMLSerializer ??= XMLSerializer
const [{ default: createLibCellML }, oxigraph, { Converter }, { bgRdfFramework }] = await Promise.all([
  load('libcellml.js'), load('oxigraph'), load('@physiomelinks/bg2cellml'), load('@celldl/editor', './bgrdf'),
])

const { version } = JSON.parse(await readFile(req.resolve('@physiomelinks/bg2cellml/package.json'), 'utf8'))
const converter = new Converter({ libcellml: await createLibCellML(), oxigraph })

// The bundled BG-RDF framework, or with BGRDF_DIR set, that checkout's ontology and templates over it.
async function framework() {
  const bundled = bgRdfFramework(), dir = process.env.BGRDF_DIR
  if (!dir) return bundled
  const templates = { ...bundled.templates }
  for (const name of (await readdir(path.join(dir, 'templates'))).filter((n) => n.endsWith('.ttl')))
    templates[name] = await readFile(path.join(dir, 'templates', name), 'utf8')
  return { ontology: await readFile(path.join(dir, 'schema', 'ontology.ttl'), 'utf8'), templates }
}

const frameworkIssues = converter.loadFramework(await framework())
const document = await readFile(documentPath, 'utf8')
const result = await converter.convert(document, { sourceUri: pathToFileURL(path.resolve(documentPath)).href })
const text = (issue) => [issue.severity, issue.stage, issue.elementId].filter(Boolean).join(' ') + `: ${issue.message}`
const issues = [...(frameworkIssues || []).map(String), ...(result.errors || []).map(text)]
if (result.exception) issues.push(`exception: ${String(result.exception).split('\n')[0]}`)
if (result.cellml) await writeFile(outPath, result.cellml)
console.log(JSON.stringify({ ok: Boolean(result.ok) && !issues.length, issues, warnings: (result.warnings || []).map(text),
                             generator: { name: '@physiomelinks/bg2cellml', version } }))
