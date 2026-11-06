"""A one-state module outside the library, for the ad hoc harness tests: dx/dt = k (|u| - x)."""
import json
import textwrap

CELLML = textwrap.dedent('''\
    <?xml version='1.0' encoding='UTF-8'?>
    <model name="relaxer_v01" xmlns="http://www.cellml.org/cellml/1.1#" xmlns:cellml="http://www.cellml.org/cellml/1.1#">
        <component name="relaxer">
            <variable name="t" public_interface="in" units="second"/>
            <variable name="k" public_interface="in" units="dimensionless"/>
            <variable name="u" public_interface="in" units="dimensionless"/>
            <variable name="x_init" public_interface="in" units="dimensionless"/>
            <variable initial_value="x_init" name="x" public_interface="out" units="dimensionless"/>
            <math xmlns="http://www.w3.org/1998/Math/MathML">
                <apply><eq/><apply><diff/><bvar><ci>t</ci></bvar><ci>x</ci></apply>
                    <apply><divide/><apply><times/><ci>k</ci><apply><minus/><apply><abs/><ci>u</ci></apply><ci>x</ci></apply></apply>
                    <cn cellml:units="second">1</cn></apply></apply>
            </math>
        </component>
    </model>
    ''')


def entry():
    return {'module_type': 'relaxer', 'module_subtype': 'v01', 'module_format': 'cellml',
            'component_file': 'relaxer_v01_modules.cellml', 'component_type': 'relaxer',
            'entrance_ports': [], 'exit_ports': [], 'general_ports': [],
            'variables_and_units': [['k', 'dimensionless', 'access', 'constant'],
                                    ['u', 'dimensionless', 'access', 'boundary_condition'],
                                    ['x_init', 'dimensionless', 'access', 'constant'],
                                    ['x', 'dimensionless', 'access', 'variable']]}


def write_system(tmp_path, u):
    """The module's files and a one-instance system "relax" (instance r, k = 2, input u, x_init = 0);
    returns (external modules dir, module array, parameters file)."""
    ext = tmp_path / 'ext'
    ext.mkdir()
    (ext / 'relaxer_v01_modules.cellml').write_text(CELLML)
    (ext / 'relaxer_v01_modules_config.json').write_text(json.dumps([entry()]))
    array = tmp_path / 'relax_module_array.json'
    array.write_text(json.dumps([{'name': 'r', 'module_type': 'relaxer', 'module_subtype': 'v01',
                                  'inp_instances': [], 'out_instances': []}]))
    params = tmp_path / 'relax_parameters.json'
    params.write_text(json.dumps({'schema_version': 1, 'r': {
        'k': {'value': 2, 'units': 'dimensionless', 'data_reference': 'test'},
        'u': {'value': u, 'units': 'dimensionless', 'data_reference': 'test'},
        'x_init': {'value': 0, 'units': 'dimensionless', 'data_reference': 'test'}}}))
    return ext, array, params
