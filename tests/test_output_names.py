"""
Output names written with a nested vessel's path (SN/soma/membrane/V) become the model's flat
names (SN__soma__membrane/V) wherever cam_testing reads them; the flat form is unchanged.
"""
import numpy as np

from cam_testing import harness, system


def test_a_path_output_becomes_the_model_name():
    assert harness.output_name('SN/soma/membrane/V') == 'SN__soma__membrane/V'
    assert harness.output_name('SN__soma__membrane/V') == 'SN__soma__membrane/V'
    assert harness.output_name('pvn/u') == 'pvn/u'
    assert harness.output_name('V') == f'{harness.VESSEL}/V'


def test_maps_and_ignore_lists_take_the_path_form():
    assert harness.model_output_map({'soma_SN/V': 'SN/soma/membrane/V'}) == {'soma_SN/V': 'SN__soma__membrane/V'}
    # an ignore dict's values are reasons, left as written
    assert harness.model_output_keys({'SN/soma/t': 'time (s/ms)'}) == {'SN__soma/t': 'time (s/ms)'}
    assert harness.model_output_keys(['heart/clock/mt']) == ['heart__clock/mt']


def test_compare_maps_a_reference_output_to_a_path():
    ref = {'soma_SN/V': np.array([1.0, 2.0]), 'soma_SN/t': np.array([0.0, 1.0])}
    new = {'SN__soma__membrane/V': np.array([1.0, 2.0])}
    rows, missing = system.compare(ref, new, {'soma_SN/V': 'SN/soma/membrane/V'}, 1e-12,
                                   ignore={'soma_SN/t': 'time'})
    assert missing == [] and [(r['model'], r['ok']) for r in rows] == [('SN__soma__membrane/V', True)]
