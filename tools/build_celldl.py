"""Writes the CellDL specs of the detailed metabolism modules' bond graphs (tools/celldl.py), for the CellDL editor's render
tool to draw:

    python tools/build_celldl.py OUT_DIR [module ...]

OUT_DIR/<module>.json is the module at rest (its default parameters); OUT_DIR/<module>_run.json has its boundary
concentrations at the verification config's run_parameters (the run the module tests start off rest), for comparing
the CellML bg2cellml makes from it against the module.
"""
import contextlib
import io
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_metabolism_detailed as detailed  # noqa: E402
from celldl import celldl_spec  # noqa: E402

MODULES = {
    'glycolysis': (lambda bg: dict(detailed.GLY_REST), detailed.CYTO['V_c']),
    'mitochondrion': (detailed.mitochondrion_rest, detailed.MITO['V_c']),
    'glycogen': (detailed.glycogen_rest, detailed.CYTO['V_c']),
    'gluconeogenesis': (detailed.gluconeogenesis_rest, detailed.MITO['V_c']),
}


def run_parameters(module):
    with open(os.path.join(detailed.REPO, 'modules', 'cell', module, 'versions', detailed.VERSION,
                           f'{module}_{detailed.VERSION}_verification_config.json')) as f:
        spec = json.load(f)
    return {k[2:]: v for k, v in (spec.get('run_parameters') or {}).items() if k.startswith('C_')}


def specs(module):
    '''(rest spec, run spec) of a module.'''
    rest_of, volume = MODULES[module]
    with contextlib.redirect_stdout(io.StringIO()):
        bg, kappa = detailed.BUILDERS[module]()
    rest = rest_of(bg)
    run = dict(rest, **run_parameters(module))
    return celldl_spec(bg, rest, kappa, volume), celldl_spec(bg, run, kappa, volume)


if __name__ == '__main__':
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    for module in sys.argv[2:] or MODULES:
        rest, run = specs(module)
        for name, spec in ((module, rest), (f'{module}_run', run)):
            with open(os.path.join(out, f'{name}.json'), 'w') as f:
                json.dump(spec, f, indent=1)
        print(f'{module}: {len(rest["elements"])} elements, {len(rest["connections"])} connections; folded: '
              f'{len(rest["folded"])}')
