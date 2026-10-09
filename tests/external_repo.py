"""
A tiny repo laid out like this one but outside it, for the tests of cam_testing in another repository
(tests/test_external_repo.py): one module_type, ext_outlet_pressure (this library's
outlet_pressure/constant under a new name), whose harness puts it downstream of this library's
inlet_flow and arterial_simple; and a system model of the same three vessels, whose reference is a
ready CellML model. The library's modules/ is an extra library (pyproject.toml [tool.cam_testing]).
"""
import json
import os
import shutil

from cam_testing import paths

SOURCE = os.path.join('boundary_conditions', 'outlet_pressure', 'versions', 'constant')
MT, VERSION = 'ext_outlet_pressure', 'constant'
OLD_COMPONENT = 'constant_pressure_BC_type__outlet_pressure_nn_constant'
NEW_COMPONENT = 'constant_pressure_BC_type__ext_outlet_pressure_constant'
SYSTEM_CATEGORY, SYSTEM = 'demo', 'ext_chain'

SYSTEM_RECORDS = [
    {'name': 'flow_in', 'module_type': 'inlet_flow', 'module_subtype': 'constant', 'parameterisation': 'default',
     'inp_instances': [], 'out_instances': ['va']},
    {'name': 'va', 'module_type': 'arterial_simple', 'module_subtype': 'vp', 'parameterisation': 'default',
     'inp_instances': ['flow_in'], 'out_instances': ['out']},
    {'name': 'out', 'module_type': MT, 'module_subtype': VERSION, 'parameterisation': 'default',
     'inp_instances': ['va'], 'out_instances': []},
]


def library_modules_dir():
    '''This library's modules/ (the checkout cam_testing is installed from).'''
    return os.path.join(paths.PACKAGE_CHECKOUT, 'modules')


def make(root, library=None, pyproject=True):
    '''Writes the repo at ``root``; returns root. ``library``: the extra library's modules/ (default this one's).'''
    library = library or library_modules_dir()
    src = os.path.join(library, SOURCE)
    vdir = os.path.join(root, 'modules', 'boundary_conditions', MT, 'versions', VERSION)
    os.makedirs(os.path.join(vdir, 'parameterisations', 'default'))
    with open(os.path.join(root, 'modules', 'README.md'), 'w') as f:
        f.write('# modules\n')

    def text(name):
        with open(os.path.join(src, name.replace(MT, 'outlet_pressure'))) as f:
            return f.read().replace(OLD_COMPONENT, NEW_COMPONENT)

    stem = f'{MT}_{VERSION}'
    for suffix in ('modules.cellml', 'units.cellml', 'references.bib'):
        with open(os.path.join(vdir, f'{stem}_{suffix}'), 'w') as f:
            f.write(text(f'{stem}_{suffix}'))
    config = json.loads(text(f'{stem}_modules_config.json'))
    config[0].update({'module_type': MT, 'component_file': f'{stem}_modules.cellml'})
    with open(os.path.join(vdir, f'{stem}_modules_config.json'), 'w') as f:
        json.dump(config, f, indent=2)
    verification = json.loads(text(f'{stem}_verification_config.json'))
    verification['module_type'] = MT
    # the harness: this library's inlet_flow and arterial_simple, upstream of the version under test
    assert [r[2] for r in verification['harness']['module_array']] == ['inlet_flow', 'arterial_simple', 'outlet_pressure']
    verification['harness']['module_array'][2][2] = MT
    with open(os.path.join(vdir, f'{stem}_verification_config.json'), 'w') as f:
        json.dump(verification, f, indent=2)
    with open(os.path.join(vdir, f'{stem}_tests.yaml'), 'w') as f:
        f.write(f'module_type: {MT}\nversion: {VERSION}\nreviewed: true\nnotes: a test copy of outlet_pressure/constant\n')
    shutil.copy(os.path.join(src, 'parameterisations', 'default', 'default_parameters.csv'),
                os.path.join(vdir, 'parameterisations', 'default', 'default_parameters.csv'))

    sdir = os.path.join(root, 'system_models', SYSTEM_CATEGORY, SYSTEM)
    os.makedirs(os.path.join(sdir, 'reference'))
    with open(os.path.join(sdir, f'{SYSTEM}_module_array.json'), 'w') as f:
        json.dump(SYSTEM_RECORDS, f, indent=1)
    with open(os.path.join(sdir, f'{SYSTEM}_parameters.csv'), 'w') as f:
        f.write('variable_name,units,value,data_reference\n')
        for name, units, value, ref in verification['harness']['parameters']:
            f.write(f'{name},{units},{value},{ref}\n')
    with open(os.path.join(sdir, f'{SYSTEM}_system.yaml'), 'w') as f:
        f.write(f'model: {SYSTEM}\ncategory: {SYSTEM_CATEGORY}\nreviewed: false\nsim_time: 2.0\npre_time: 0.0\n'
                f'dt: 0.01\nequivalence:\n  reference_cellml: reference/{SYSTEM}.cellml\n  tol: 1.0e-06\n'
                f'  solver_info: {{rtol: 1.0e-10, atol: 1.0e-12}}\n')
    if pyproject:
        with open(os.path.join(root, 'pyproject.toml'), 'w') as f:
            f.write(f'[tool.cam_testing]\nmodule_library_dirs = [{json.dumps(os.path.dirname(library))}]\n')
    return root


def write_reference_cellml(system, work_dir):
    '''A ready CellML reference for a system model, as PhLynx would export it: the model generated from
    the libraries, flattened into one file with each vessel's component named by the vessel
    ("va", not libcuflynx's "va_module").'''
    import libcellml
    from cam_testing import system as systems
    path = systems.generate(system, work_dir)
    parser, importer = libcellml.Parser(False), libcellml.Importer(False)
    with open(path) as f:
        model = parser.parseModel(f.read())
    importer.resolveImports(model, os.path.dirname(path) + os.sep)
    flat = importer.flattenModel(model)
    def components(parent):
        for i in range(parent.componentCount()):
            c = parent.component(i)
            yield c
            yield from components(c)

    every = list(components(flat))
    names = {c.name() for c in every}
    for c in every:     # the vessel's inner component first, then its wrapper takes the vessel's name
        if c.name() + '_module' in names:
            c.setName(c.name() + '_inner')
    for c in every:
        if c.name().endswith('_module'):
            c.setName(c.name()[:-len('_module')])
    out = os.path.join(system.dir, 'reference', f'{system.name}.cellml')
    with open(out, 'w') as f:
        f.write(libcellml.Printer().printModel(flat))
    return out
