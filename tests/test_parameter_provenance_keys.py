"""A parameterisation's parameters may say what set each value: libcuflynx's origin, and reason."""
import json

import pytest

pytest.importorskip('libcuflynx.utilities.parameters_file')
from cam_testing.suite.test_structure import _parameters_file_problems  # noqa: E402


def _write(tmp_path, leaf):
    path = tmp_path / 'p_parameters.json'
    path.write_text(json.dumps({'schema_version': 1, 'C': {'value': 1, 'units': 'dim', 'data_reference': 'ref',
                                                           'sourced': 'yes', **leaf}}))
    return str(path)


def test_origin_and_reason_are_allowed(tmp_path):
    assert _parameters_file_problems(_write(tmp_path, {'origin': 'measured'})) == []
    assert _parameters_file_problems(_write(tmp_path, {'origin': 'set_by_hand', 'reason': 'why'})) == []


def test_other_keys_are_not(tmp_path):
    problems = _parameters_file_problems(_write(tmp_path, {'note': 'x'}))
    assert len(problems) == 1 and "['note']" in problems[0]
