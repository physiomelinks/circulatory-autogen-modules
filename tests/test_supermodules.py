"""
Supermodules (modules/supermodules/<name>/): a module_format "supermodule" config entry whose
submodules libcuflynx expands in a model.

    pytest tests/test_supermodules.py

  supermodule_structure_test    every submodule is a library (module_type, module_subtype), internal
                                connections name sibling submodules, and the default parameters
                                name submodules (or are globals)
  supermodule_equivalence_test  each spec.tests.equivalent entry: a system model that uses the
                                supermodule reproduces a system model with the submodules written
                                out, output for output (submodule outputs renamed <instance>_<sub>)
"""
import csv

import pytest

from cam_testing import library, supermodule as sm, system as systems

pytestmark = pytest.mark.system_model

SUPERMODULES = sm.supermodules()
EQUIVALENT = [pytest.param(s, e, id=f'{s.name}/{e["model"]}')
              for s in SUPERMODULES for e in ((s.spec.get('tests') or {}).get('equivalent') or [])]


@pytest.mark.parametrize('s', [pytest.param(s, id=s.name) for s in SUPERMODULES])
def supermodule_structure_test(s):
    known = {(e['vessel_type'], e['BC_type']) for n in library.module_names() for e in library.load_module(n).config}
    # a submodule may itself be a supermodule (nesting)
    known |= {(o.entry['module_type'], o.entry['module_subtype']) for o in SUPERMODULES}
    names = set(s.submodule_names)
    problems = []
    if len(names) != len(s.submodules):
        problems.append('duplicate submodule names')
    for sub in s.submodules:
        rec = library.normalise_config_entry({'module_type': sub['module_type'], 'module_subtype': sub['module_subtype'],
                                              'component_file': '-', 'component_type': '-'})
        if (rec['vessel_type'], rec['BC_type']) not in known:
            problems.append(f'{sub["name"]}: ({sub["module_type"]}, {sub["module_subtype"]}) is not in the module library')
        for n in sub.get('inp_instances', []) + sub.get('out_instances', []):
            if n not in names:
                problems.append(f'{sub["name"]}: {n} is not a submodule (external connections belong in the host\'s '
                                'per_submodule_inputs / per_submodule_outputs)')
    if s.defaults_path:
        with open(s.defaults_path) as f:
            for row in csv.DictReader(f):
                var = (row.get('variable_name') or '').strip()
                owner = next((n for n in sorted(names, key=len, reverse=True) if var.endswith('_' + n)), None)
                if owner is None and var not in (s.spec.get('globals') or ['T', 'rho', 'l_eff']):
                    problems.append(f'default parameter {var} names no submodule and is not a declared global')
    assert not problems, '\n'.join(problems)


@pytest.mark.parametrize('s, e', EQUIVALENT)
def supermodule_equivalence_test(s, e, tmp_path):
    model = systems.load_system(e['model'])
    target = systems.load_system(e['reproduces'])
    solver_info = e.get('solver_info') or {'rtol': 1e-10, 'atol': 1e-12}
    t_ref, ref = systems.simulate(systems.generate(target, str(tmp_path / 'target')), target.spec, solver_info)
    t_new, new = systems.simulate(systems.generate(model, str(tmp_path / 'model')), model.spec, solver_info)
    # submodule outputs renamed <instance>_<sub>; an explicit output_map in the spec (for nested
    # supermodules, or variables that moved between components) takes precedence
    output_map = {**sm.prefixed_output_map(ref, e.get('instance', s.name), s.submodule_names), **(e.get('output_map') or {})}
    ignore = e.get('ignore') or {}
    ref = {k: v for k, v in ref.items() if k not in ignore}
    rows, missing = systems.compare(ref, new, output_map, float(e.get('tol', 1e-9)))
    bad = [r for r in rows if not r['ok']]
    assert rows and not bad and not missing, (
        f'{len(bad)} of {len(rows)} outputs differ; missing {missing[:10]}\n'
        + '\n'.join(f"{r['reference']} -> {r['model']}: {r['difference']:.3g}" for r in bad[:10]))
