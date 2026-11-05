"""
Unit tests of the stability test's time limits (cam_testing.checks): runs stopped at their
wall-clock limit, and the step ladder ending within its budgets. The runs are faked.
"""
import multiprocessing as mp
import time

import pytest

from cam_testing import checks

needs_fork = pytest.mark.skipif('fork' not in mp.get_all_start_methods(), reason='runs are stopped through fork')


@needs_fork
def test_a_run_past_its_limit_is_stopped():
    start = time.perf_counter()
    with pytest.raises(checks.RunTimeout):
        checks.run_limited(lambda: time.sleep(60), 0.5)
    assert time.perf_counter() - start < 5.0


@needs_fork
def test_a_run_within_its_limit_returns_and_failures_come_back():
    assert checks.run_limited(lambda: (1, {'x': [2.0]}), 10.0) == (1, {'x': [2.0]})

    def fails():
        raise ValueError('simulation failed')
    with pytest.raises(RuntimeError, match='simulation failed'):
        checks.run_limited(fails, 10.0)


def _record(rows):
    def record(cfg, fn, limit):
        start = time.perf_counter()
        works, timed_out = False, False
        try:
            works = checks.run_limited(fn, limit)
        except checks.RunTimeout:
            timed_out = True
        except RuntimeError:
            pass
        row = {'config': cfg, 'label': str(cfg), 'works': works, 'timed_out': timed_out,
               'seconds': time.perf_counter() - start}
        rows.append(row)
        return row
    return record


@needs_fork
def test_quick_failures_do_not_let_a_slow_run_go_on():
    # An explicit solver on a stiff model fails at once at large steps, then at a small step runs
    # far longer than the budget (the fermenter's RK45 at 1e-6 s took 47 min): the run is stopped.
    rows = []

    def run(step):
        if step is None or step > 1e-5:
            raise ValueError('simulation failed')
        time.sleep(60)
        return True
    start = time.perf_counter()
    last, stop = checks.refine_family(_record(rows), lambda s: {'max_step': s}, run, 1e-2, 10.0, True,
                                      budget=1.0, min_step=1e-7, deadline=time.perf_counter() + 60.0)
    assert time.perf_counter() - start < 10.0
    assert rows[-1]['timed_out'] and not last['works']
    assert 'limit' in stop


@needs_fork
def test_an_unconstrained_run_past_its_limit_ends_the_ladder():
    # the first try has no step (unconstrained); stopped at its limit, the reason names it so
    rows = []

    def run(step):
        time.sleep(60)
        return True
    last, stop = checks.refine_family(_record(rows), lambda s: {'max_step': s}, run, 1e-2, 10.0, True,
                                      budget=0.5, min_step=1e-7, deadline=time.perf_counter() + 60.0)
    assert len(rows) == 1 and rows[0]['timed_out'] and rows[0]['config']['max_step'] is None
    assert stop.startswith('stopped: run (unconstrained) reached the')


def test_the_ladder_ends_at_the_tests_deadline():
    rows = []
    last, stop = checks.refine_family(_record(rows), lambda s: {'max_step': s}, lambda s: False, 1e-2, 10.0, True,
                                      budget=10.0, min_step=1e-12, deadline=time.perf_counter() - 1.0)
    assert rows == [] and last is None
    assert 'time budget' in stop


def test_the_first_working_step_ends_the_ladder():
    rows = []
    last, stop = checks.refine_family(_record(rows), lambda s: {'max_step': s}, lambda s: s is not None and s <= 1e-3,
                                      1e-2, 10.0, True, budget=10.0, min_step=1e-7, deadline=time.perf_counter() + 60.0)
    assert last['works'] and stop == ''
    assert [r['config']['max_step'] for r in rows] == [None, 1e-2, 1e-3]


def _timed_out(rows):
    def record(cfg, fn, limit):
        rows.append({'config': cfg, 'works': False, 'timed_out': True, 'seconds': limit})
        return rows[-1]
    return record


def test_a_timed_out_unconstrained_run_is_reported_not_raised():
    # The first run has no step limit (step None); its timeout is a message, not a TypeError.
    rows = []
    last, stop = checks.refine_family(_timed_out(rows), lambda s: {'max_step': s}, None, 1e-2, 10.0, True,
                                      budget=1.0, min_step=1e-7, deadline=time.perf_counter() + 60.0)
    assert [r['config']['max_step'] for r in rows] == [None] and last['timed_out']
    assert stop == 'stopped: run (unconstrained) reached the 1 s limit'


def test_a_timed_out_stepped_run_names_its_step():
    rows = []
    last, stop = checks.refine_family(_timed_out(rows), lambda s: {'max_step': s}, None, 1e-2, 10.0, False,
                                      budget=1.0, min_step=1e-7, deadline=time.perf_counter() + 60.0)
    assert stop == 'stopped: run (step 0.01) reached the 1 s limit'
    assert checks.step_label(None) == 'unconstrained' and checks.step_label(1e-3) == 'step 0.001'
