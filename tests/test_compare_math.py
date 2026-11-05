"""cam_testing.compare_math: two libraries' versions of a module compared equation by equation."""
import json
import os

import yaml

from cam_testing import compare_math

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

# dx/dt = k (S - x) / (K + S), y = 2 x
CELLML = '''<?xml version='1.0' encoding='UTF-8'?>
<model name="m" xmlns="http://www.cellml.org/cellml/1.1#" xmlns:cellml="http://www.cellml.org/cellml/1.1#">
    <component name="comp">
        <variable name="t" public_interface="in" units="second"/>
        <variable name="{k}" public_interface="in" units="per_second"/>
        <variable name="S" public_interface="in" units="mM"/>
        <variable name="K" public_interface="in" units="mM"/>
        <variable name="{x}" initial_value="1" public_interface="out" units="mM"/>
        <variable name="y" public_interface="out" units="mM"/>
        <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><eq/><apply><diff/><bvar><ci>t</ci></bvar><ci>{x}</ci></apply>{rate}</apply>
            {y_eq}
        </math>
    </component>
</model>
'''
RATE = ('<apply><divide/><apply><times/><ci>{k}</ci><apply><minus/><ci>S</ci><ci>{x}</ci></apply></apply>'
        '<apply><plus/><ci>K</ci><ci>S</ci></apply></apply>')
# the same rate, rearranged: (k S - k x) / (S + K)
RATE_REARRANGED = ('<apply><divide/><apply><minus/><apply><times/><ci>{k}</ci><ci>S</ci></apply>'
                   '<apply><times/><ci>{x}</ci><ci>{k}</ci></apply></apply>'
                   '<apply><plus/><ci>S</ci><ci>K</ci></apply></apply>')
# a different law: k (S - x)
RATE_LINEAR = '<apply><times/><ci>{k}</ci><apply><minus/><ci>S</ci><ci>{x}</ci></apply></apply>'
Y_EQ = '<apply><eq/><ci>y</ci><apply><times/><cn cellml:units="dimensionless">2.0</cn><ci>{x}</ci></apply></apply>'


def _version(root, path, rate=RATE, y=True, k='k', x='x', params=None):
    vdir = compare_math.version_dir(str(root), path)
    os.makedirs(os.path.join(vdir, 'parameterisations', 'default'))
    names = {'k': k, 'x': x}
    body = CELLML.format(rate=rate.format(**names), y_eq=Y_EQ.format(**names) if y else '', **names)
    with open(os.path.join(vdir, 'm_modules.cellml'), 'w') as f:
        f.write(body)
    with open(os.path.join(vdir, 'm_modules_config.json'), 'w') as f:
        json.dump([{'component_file': 'm_modules.cellml', 'component_type': 'comp'}], f)
    with open(os.path.join(vdir, 'parameterisations', 'default', 'default_parameters.csv'), 'w') as f:
        f.write('variable_name,units,value,data_reference,sourced\n')
        for name, value in (params or {k: 1.0, 'K': 0.5}).items():
            f.write(f'{name},mM,{value},design,no\n')
    return vdir


def _pairs(tmp_path, **pair):
    path = tmp_path / 'pairs.yaml'
    path.write_text(yaml.safe_dump({'pairs': [{'cam': 'cell/m/v1', 'here': 'cell/m/v1', **pair}]}))
    return str(path)


def test_the_same_maths_under_other_names_and_term_orders_has_no_differences(tmp_path):
    _version(tmp_path / 'a', 'cell/m/v1')
    _version(tmp_path / 'b', 'cell/m/v1', rate=RATE_REARRANGED, k='k_cm', x='x_cm', params={'k_cm': 3.0, 'K': 0.5})
    (r,) = compare_math.compare(_pairs(tmp_path, rename={'k_cm': 'k', 'x_cm': 'x'}), tmp_path / 'a', tmp_path / 'b')
    assert r['differences'] == []
    # a parameter value that differs is reported for information only
    assert r['parameters'] == [{'name': 'k', 'a': 1.0, 'b': 3.0}]


def test_a_different_rate_law_is_unexplained_until_listed(tmp_path):
    _version(tmp_path / 'a', 'cell/m/v1')
    _version(tmp_path / 'b', 'cell/m/v1', rate=RATE_LINEAR)
    (r,) = compare_math.compare(_pairs(tmp_path), tmp_path / 'a', tmp_path / 'b')
    (d,) = r['differences']
    assert (d['equation'], d['kind'], d['status']) == ('d/dt x', 'differs', 'unexplained')
    assert d['a'] == 'k*(S - x)/(K + S)' and d['b'] == 'k*(S - x)'
    listed = [{'variable': 'x', 'reason': 'no saturation in the cardiomyocyte'}]
    (r,) = compare_math.compare(_pairs(tmp_path, divergences=listed), tmp_path / 'a', tmp_path / 'b')
    assert [(d['status'], d['reason']) for d in r['differences']] == [('explained', 'no saturation in the cardiomyocyte')]


def test_an_equation_only_one_side_has_and_a_stale_divergence_are_reported(tmp_path):
    _version(tmp_path / 'a', 'cell/m/v1')
    _version(tmp_path / 'b', 'cell/m/v1', y=False)
    listed = [{'variable': 'v_gone', 'reason': 'was different once'}]
    (r,) = compare_math.compare(_pairs(tmp_path, divergences=listed), tmp_path / 'a', tmp_path / 'b')
    assert [(d['equation'], d['kind'], d['status']) for d in r['differences']] == [
        ('y', 'only in cell/m/v1', 'unexplained'), ('v_gone', 'listed, no longer differs', 'stale')]


def test_the_report_never_fails_unless_strict(tmp_path):
    _version(tmp_path / 'a', 'cell/m/v1')
    _version(tmp_path / 'b', 'cell/m/v1', rate=RATE_LINEAR)
    args = ['--pairs', _pairs(tmp_path), '--a-root', str(tmp_path / 'a'), '--b-root', str(tmp_path / 'b'),
            '--markdown', str(tmp_path / 'report.md'), '--json', str(tmp_path / 'report.json'),
            '--a-label', 'CAM', '--b-label', 'cardiometabolism']
    assert compare_math.main(args) == 0
    assert compare_math.main(args + ['--strict']) == 1
    report = (tmp_path / 'report.md').read_text()
    assert '1 unexplained, 0 explained, 0 stale' in report and 'cardiometabolism: `k*(S - x)`' in report
    assert json.loads((tmp_path / 'report.json').read_text())[0]['differences'][0]['equation'] == 'd/dt x'


def test_a_library_module_compared_with_itself_has_no_differences(tmp_path):
    '''A library version (Paci2013's cardiomyocyte Ca handling, with piecewise and exponential terms) parses whole
    and matches itself.'''
    path = 'cell/Ca_handling/cardiomyocyte_Paci2013_v01'
    pairs = tmp_path / 'pairs.yaml'
    pairs.write_text(yaml.safe_dump({'pairs': [{'a': path, 'b': path}]}))
    (r,) = compare_math.compare(str(pairs), REPO, REPO)
    assert r['differences'] == [] and r['parameters'] == []
    assert len(compare_math.version_equations(compare_math.version_dir(REPO, path))) > 5


def test_a_missing_version_is_named(tmp_path):
    pairs = tmp_path / 'pairs.yaml'
    pairs.write_text(yaml.safe_dump({'pairs': [{'a': 'cell/none/v1', 'b': 'cell/none/v1'}]}))
    try:
        compare_math.compare(str(pairs), str(tmp_path), str(tmp_path))
    except FileNotFoundError as e:
        assert 'cell/none/versions/v1' in str(e)
    else:
        raise AssertionError('no error for a missing version')
