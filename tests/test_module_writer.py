"""tools/module_writer.py: what a rebuild of a generated module version keeps from the version it replaces."""
import json
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'tools'))
from module_writer import kept_record_keys  # noqa: E402


def test_a_rebuild_keeps_the_licence_and_creator_given_in_review(tmp_path):
    path = tmp_path / 'm_v_modules_config.json'
    path.write_text(json.dumps([{'module_type': 'm', 'licence': 'CC-BY-4.0', 'creator': ['A. Author']}]))
    assert kept_record_keys(str(path)) == {'licence': 'CC-BY-4.0', 'creator': ['A. Author']}


def test_a_new_version_gets_the_default_licence_and_no_creator(tmp_path):
    assert kept_record_keys(str(tmp_path / 'missing_modules_config.json')) == {'licence': 'CC0-1.0', 'creator': []}
