"""The CellDL documents of the modules' bond graphs (tools/celldl.py): specs for the CellDL editor's render tool, and
the rendered documents installed in the versions.

    python tools/build_celldl.py specs OUT_DIR [module ...]
    (cd CellDLEditorComponent && bun run render-celldl RENDER_DIR OUT_DIR/<module>.json ...)
    python tools/build_celldl.py install OUT_DIR RENDER_DIR [module ...]

`specs` writes OUT_DIR/<module>.json, the module at rest (its default parameters; the lumped cytosol at its export
point), and for the detailed modules OUT_DIR/<module>_run.json, with the boundary concentrations of the verification
config's run_parameters (tools/compare_celldl.py runs both against the module). `install` writes each version's
<module_type>_<version>_bondgraph.ttl (the document's Turtle) and _bondgraph.svg (the rendered CellDL document,
RENDER_DIR/<module>.celldl.svg).
"""
import contextlib
import io
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_metabolism_detailed as detailed  # noqa: E402
from celldl import celldl_spec  # noqa: E402
from module_writer import REPO  # noqa: E402


def _detailed(module, rest_of):
    def build():
        with contextlib.redirect_stdout(io.StringIO()):
            bg, kappa = detailed.BUILDERS[module]()
        return bg, kappa, rest_of(bg)
    return build


def _cytosol():
    import build_bondgraph_lumped as lumped
    with contextlib.redirect_stdout(io.StringIO()):
        bg, kappa, _, conc = lumped.build_cytosol()
    return bg, kappa, conc


# module -> (module_type, version, build -> (bond graph, kappa [mM/s], concentrations [mM]), volume [m3], detailed)
MODULES = {
    'glycolysis': ('glycolysis', detailed.VERSION, _detailed('glycolysis', lambda bg: dict(detailed.GLY_REST)),
                   detailed.CYTO['V_c'], True),
    'mitochondrion': ('mitochondrion', detailed.VERSION, _detailed('mitochondrion', detailed.mitochondrion_rest),
                      detailed.MITO['V_c'], True),
    'glycogen': ('glycogen', detailed.VERSION, _detailed('glycogen', detailed.glycogen_rest), detailed.CYTO['V_c'], True),
    'gluconeogenesis': ('gluconeogenesis', detailed.VERSION, _detailed('gluconeogenesis', detailed.gluconeogenesis_rest),
                        detailed.MITO['V_c'], True),
    'energy_metabolism_adaptive': ('energy_metabolism', 'adaptive_ArgusUNPUBLISHED_v01', _cytosol, detailed.CYTO['V_c'],
                                   False),
}


def version_dir(module):
    module_type, version = MODULES[module][:2]
    return os.path.join(REPO, 'modules', 'cell', module_type, 'versions', version), f'{module_type}_{version}'


def run_parameters(module):
    d, stem = version_dir(module)
    with open(os.path.join(d, f'{stem}_verification_config.json')) as f:
        spec = json.load(f)
    return {k[2:]: v for k, v in (spec.get('run_parameters') or {}).items() if k.startswith('C_')}


def specs(module):
    '''{name: spec}: the module at rest, and for a detailed module at its run_parameters.'''
    _, _, build, volume, is_detailed = MODULES[module]
    bg, kappa, rest = build()
    out = {module: celldl_spec(bg, rest, kappa, volume, module=module)}
    if is_detailed:
        out[f'{module}_run'] = celldl_spec(bg, dict(rest, **run_parameters(module)), kappa, volume, module=module)
    return out


def install(spec_dir, render_dir, module):
    d, stem = version_dir(module)
    with open(os.path.join(spec_dir, f'{module}.json')) as f:
        turtle = json.load(f)['turtle']
    with open(os.path.join(d, f'{stem}_bondgraph.ttl'), 'w') as f:
        f.write(turtle)
    shutil.copyfile(os.path.join(render_dir, f'{module}.celldl.svg'), os.path.join(d, f'{stem}_bondgraph.svg'))


if __name__ == '__main__':
    command, out, *rest = sys.argv[1:]
    if command == 'specs':
        os.makedirs(out, exist_ok=True)
        for module in rest or MODULES:
            for name, spec in specs(module).items():
                with open(os.path.join(out, f'{name}.json'), 'w') as f:
                    json.dump(spec, f, indent=1)
            print(f'{module}: {len(spec["elements"])} elements, {len(spec["connections"])} connections')
    elif command == 'install':
        render_dir, *modules = rest
        for module in modules or MODULES:
            install(out, render_dir, module)
            print(f'{module}: installed in {version_dir(module)[0]}')
    else:
        sys.exit(__doc__)
