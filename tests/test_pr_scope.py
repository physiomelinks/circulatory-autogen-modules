'''tools/pr_scope.py: a module PR ("Module:" title) changes only module content.'''
import importlib.util
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location('pr_scope', os.path.join(ROOT, 'tools', 'pr_scope.py'))
pr_scope = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pr_scope)

MODULE_FILES = [
    'modules/heart/chamber/versions/Liang2009_v01/chamber_Liang2009_v01_modules.cellml',
    'modules/heart/chamber/versions/Liang2009_v01/parameterisations/lv/lv_obs_data.json',
    'system_models/closed_loop_cvs/3compartment/3compartment_module_config.json',
    'manifests/all.json',
    'reviews/BG_review.yaml',
]


def test_module_pr_with_only_module_content_passes():
    assert pr_scope.out_of_scope('Module: heart chamber review', MODULE_FILES) == []


def test_module_pr_lists_repository_files():
    bad = ['cam_testing/checks.py', 'modules/README.md', 'modules/directory_schema.json',
           '.github/workflows/module-tests.yml', 'README.md', 'tests/test_modules.py']
    assert pr_scope.out_of_scope('Module: heart chamber review', MODULE_FILES + bad) == sorted(bad)


def test_other_pr_lists_module_content():
    repo = ['cam_testing/report.py', 'modules/README.md', 'modules/directory_schema.json', 'Makefile']
    assert pr_scope.out_of_scope('Reports: partial in blue', repo) == []
    assert pr_scope.out_of_scope('Reports: partial in blue', repo + MODULE_FILES) == sorted(MODULE_FILES)


def test_module_prefix_only_at_the_start():
    # the prefix only counts at the start of the title
    assert not pr_scope.is_module_pr('Fix the Module: label')
    assert pr_scope.is_module_pr('  Module: valve review')
