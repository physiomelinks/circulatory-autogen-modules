/**
 * Makes CellML from a CellDL bond graph with PhLynx's converter, for cam_testing's bondgraph_cellml_test:
 *
 *   node cam_testing/bridges/celldl/convert.mjs <document.svg> <out.cellml>
 *
 * The converter is @physiomelinks/bg2cellml (the celldl2cellml repo), with the BG-RDF framework the CellDL editor
 * bundles (@celldl/editor/bgrdf), at the versions PhLynx pins (package.json). Prints a JSON report on stdout:
 * { ok, issues, generator }. Exits 0 when the converter ran (ok says whether it made CellML), 3 when the bridge
 * can't run (its packages aren't installed: npm ci here, or set CELLDL_BRIDGE_DIR to a directory where they are).
 */
import { readFile, writeFile } from 'node:fs/promises'
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
const frameworkIssues = converter.loadFramework(bgRdfFramework())
const document = await readFile(documentPath, 'utf8')
const result = await converter.celldl2cellml(document, { sourceUri: pathToFileURL(path.resolve(documentPath)).href })
const issues = [...(frameworkIssues || []), ...(result.issues || [])].map(String)
if (result.exception) issues.push(`exception: ${String(result.exception).split('\n')[0]}`)
if (result.cellml) await writeFile(outPath, result.cellml)
console.log(JSON.stringify({ ok: Boolean(result.cellml) && !issues.length, issues,
                             generator: { name: '@physiomelinks/bg2cellml', version } }))
