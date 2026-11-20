"""
Reduced versions (cam_testing/reduction.py): the annotation rules, and reduction_test on a tiny repo
(tests/external_repo.py) whose version has a copy declared as reduced from it. The copy has the parent's
equations, so it matches until the parent is given other values.
"""
import json
import os
import shutil
from types import SimpleNamespace

import pytest

import external_repo
from cam_testing import checks, paths, reduction
from cam_testing.library import load_version

MT, PARENT, REDUCED = external_repo.MT, external_repo.VERSION, 'constant_reduced'


# ---- the annotations ------------------------------------------------------------------------------

def _version(key='m/reduced', parameterisations=('default',), **spec):
    return SimpleNamespace(key=key, spec=spec, parameterisation_names=lambda: list(parameterisations))


def _library(**versions):
    def load(mt, v):
        return versions[v]
    return load


ENTRY = {'module_type': 'm', 'version': 'full', 'outputs': ['y'], 'rtol': 1e-6}


def _problems(load=None, **spec):
    return reduction.annotation_problems(_version(**spec), load or _library(full=_version(key='m/full')))


def test_a_version_with_no_annotation_has_no_problem_and_is_not_reduced():
    assert _problems() == []
    assert not reduction.is_reduced(_version())


def test_complete_annotations_have_no_problem():
    assert _problems(reduction_method='symbolic_steady_state', reduced_from=[ENTRY]) == []
    other = {'name': 'balanced_truncation', 'description': 'the least controllable and observable states removed'}
    assert _problems(reduction_method=other, reduced_from=[ENTRY]) == []


@pytest.mark.parametrize('spec, message', [
    ({'reduction_method': 'operator_inference'}, 'reduction_method without reduced_from'),
    ({'reduced_from': [ENTRY]}, 'reduced_from without reduction_method'),
    ({'reduction_method': 'by_eye', 'reduced_from': [ENTRY]}, "'by_eye' is not one of"),
    ({'reduction_method': {'name': 'by_eye'}, 'reduced_from': [ENTRY]}, "'by_eye' is not one of"),
    ({'reduction_method': 'manual_term_removal', 'reduced_from': [{'version': 'full'}]}, 'needs module_type and version'),
    ({'reduction_method': 'manual_term_removal', 'reduced_from': [dict(ENTRY, version='reduced')]},
     'not reduced from itself'),
    ({'reduction_method': 'manual_term_removal', 'reduced_from': [dict(ENTRY, version='gone')]}, 'no such version'),
    ({'reduction_method': 'manual_term_removal', 'reduced_from': [{k: v for k, v in ENTRY.items() if k != 'rtol'}]},
     'needs rtol'),
    ({'reduction_method': 'manual_term_removal', 'reduced_from': [{k: v for k, v in ENTRY.items() if k != 'outputs'}]},
     'needs outputs'),
    ({'reduction_method': 'manual_term_removal', 'reduced_from': [dict(ENTRY, tolerance=1)]}, 'unknown key(s) tolerance'),
    ({'reduction_method': 'manual_term_removal', 'reduced_from': [dict(ENTRY, compare='by_eye')]}, 'compare is'),
    ({'reduction_method': 'manual_term_removal', 'reduced_from': [dict(ENTRY, parameterisation='fitted')]},
     "parameterisation 'fitted' is not a parameterisation of m/full"),
])
def test_incomplete_annotations_say_what_is_missing(spec, message):
    problems = _problems(**spec)
    assert any(message in p for p in problems), problems


def test_a_cycle_of_reductions_is_a_problem():
    full = _version(key='m/full', reduction_method='manual_term_removal',
                    reduced_from=[dict(ENTRY, version='reduced')])
    problems = _problems(_library(full=full, reduced=_version()), reduction_method='manual_term_removal',
                         reduced_from=[ENTRY])
    assert any('a cycle' in p for p in problems), problems


def test_the_default_comparison_follows_the_method():
    steady = _version(reduction_method='symbolic_steady_state')
    assert reduction._compare(steady, {}) == 'steady_state'
    assert reduction._compare(_version(reduction_method='operator_inference'), {}) == 'time_course'
    assert reduction._compare(_version(reduction_method={'name': 'x', 'description': 'y'}), {}) == 'time_course'
    assert reduction._compare(steady, {'compare': 'time_course'}) == 'time_course'


def test_every_listed_method_is_documented():
    with open(os.path.join(paths.PACKAGE_CHECKOUT, 'modules', 'README.md')) as f:
        readme = f.read()
    for method in reduction.METHODS:
        assert f'| `{method}` |' in readme, method


# ---- the test, on a repo with a version and its copy ------------------------------------------------

pytest.importorskip('libcuflynx')


def _copy_version(root, annotations):
    '''A copy of the repo's version under the name REDUCED, with ``annotations`` in its verification config.'''
    versions = os.path.join(root, 'modules', 'boundary_conditions', MT, 'versions')
    src, dst = os.path.join(versions, PARENT), os.path.join(versions, REDUCED)
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns('results', 'plots', '*.html'))
    for name in os.listdir(dst):
        if name.startswith(f'{MT}_{PARENT}_'):
            os.rename(os.path.join(dst, name), os.path.join(dst, name.replace(f'{MT}_{PARENT}_', f'{MT}_{REDUCED}_')))
    stem = os.path.join(dst, f'{MT}_{REDUCED}')
    with open(f'{stem}_modules_config.json') as f:
        config = json.load(f)
    config[0].update(module_subtype=REDUCED, component_file=f'{MT}_{REDUCED}_modules.cellml')
    with open(f'{stem}_modules_config.json', 'w') as f:
        json.dump(config, f, indent=2)
    with open(f'{stem}_verification_config.json') as f:
        verification = json.load(f)
    verification['version'] = REDUCED
    verification['harness']['module_array'][2][1] = REDUCED
    verification.update(annotations)
    with open(f'{stem}_verification_config.json', 'w') as f:
        json.dump(verification, f, indent=2)
    with open(f'{stem}_tests.yaml', 'w') as f:
        f.write(f'module_type: {MT}\nversion: {REDUCED}\nreviewed: true\nnotes: a copy of {PARENT}\n')


def _shifted_parent(root, pressure):
    '''A parameterisation `shifted` of the parent, with another outlet pressure.'''
    pdir = os.path.join(root, 'modules', 'boundary_conditions', MT, 'versions', PARENT, 'parameterisations')
    os.makedirs(os.path.join(pdir, 'shifted'), exist_ok=True)
    with open(os.path.join(pdir, 'default', 'default_parameters.json')) as f:
        parameters = json.load(f)
    parameters['P']['value'] = pressure
    with open(os.path.join(pdir, 'shifted', 'shifted_parameters.json'), 'w') as f:
        json.dump(parameters, f, indent=1)


@pytest.fixture
def repo(tmp_path):
    root = str(tmp_path / 'repo')
    external_repo.make(root)
    paths.configure(repo_root=root, library_dirs=[external_repo.library_modules_dir()], export=False)
    yield root
    paths.reset()


def _check(root, tmp_path, method='symbolic_steady_state', **entry):
    from cam_testing import library
    _copy_version(root, {'reduction_method': method,
                         'reduced_from': [dict({'module_type': MT, 'version': PARENT, 'outputs': ['va/u'],
                                                'rtol': 1e-6}, **entry)]})
    library.clear_caches()
    return reduction.reduction_check(load_version(MT, REDUCED), str(tmp_path / 'work'))


@pytest.mark.integration
def test_a_version_that_matches_its_parent_passes_at_every_point(repo, tmp_path):
    result = _check(repo, tmp_path, points=[{}, {'P': 2000.0}], sweep={'ranges': {'P': {'values': [500.0, 4000.0]}}})
    assert result.status == checks.PASSED, (result.message, result.details)
    parent = result.metrics['parents'][0]
    assert (parent['parent'], parent['compare'], parent['points']) == (f'{MT}/{PARENT}', 'steady_state', 4)
    # the points reached both models: the vessel settles at the outlet pressure plus its resistive drop
    settled = [row['outputs']['va/u'] for row in parent['rows']]
    assert [round(s['parent']) for s in settled] == [2066, 3000, 1500, 5000]
    assert all(s['reduced'] == pytest.approx(s['parent'], rel=1e-9) for s in settled)
    assert result.metrics['method'] == 'symbolic_steady_state'
    saved = checks.load(load_version(MT, REDUCED), 'reduction_test')
    assert saved.status == checks.PASSED and saved.plots == ['plots/reduction.png']


@pytest.mark.integration
def test_a_sweep_starts_from_the_first_point(repo, tmp_path):
    # the factors of a sweep with no range multiply the first point's value, not the parameterisation's
    result = _check(repo, tmp_path, points=[{'P': 2000.0}], sweep={'ranges': {'P': {'values': [500.0]}}})
    assert result.status == checks.PASSED, (result.message, result.details)
    assert [row['point'] for row in result.metrics['parents'][0]['rows']] == [{'P': 2000.0}, {'P': 500.0}]


@pytest.mark.integration
def test_a_time_course_is_compared_at_every_output_time(repo, tmp_path):
    result = _check(repo, tmp_path, method='operator_inference', sim_time=0.5, dt=0.01)
    assert result.status == checks.PASSED, (result.message, result.details)
    assert result.metrics['parents'][0]['compare'] == 'time_course'


@pytest.mark.integration
def test_a_version_that_differs_from_its_parent_fails_and_says_where(repo, tmp_path):
    _shifted_parent(repo, 1166.0)
    result = _check(repo, tmp_path, parameterisation='shifted')
    assert result.status == checks.FAILED
    # 100 Pa apart at about 2166 Pa
    assert 'va/u differs from the parent\'s va/u by 100' in result.details[0], result.details
    assert result.metrics['parents'][0]['worst_relative_difference'] == pytest.approx(100 / 2166, rel=1e-3)
    # a tolerance that allows the difference passes: the comparison, not the names, decides
    assert _check(repo, tmp_path, parameterisation='shifted', rtol=0.05).status == checks.PASSED


@pytest.mark.integration
def test_a_steady_state_comparison_needs_settled_runs(repo, tmp_path):
    result = _check(repo, tmp_path, sim_time=0.05, dt=0.001)
    assert result.status == checks.FAILED
    assert any('has not settled' in d for d in result.details), result.details


@pytest.mark.integration
def test_a_version_that_is_not_reduced_is_not_applicable(repo, tmp_path):
    result = reduction.reduction_check(load_version(MT, PARENT), str(tmp_path / 'work'))
    assert (result.status, result.message) == (checks.NOT_APPLICABLE, reduction.NOT_REDUCED)


@pytest.mark.integration
def test_wrong_annotations_fail_the_test(repo, tmp_path):
    result = _check(repo, tmp_path, version='nowhere')
    assert result.status == checks.FAILED and 'no such version' in result.message
