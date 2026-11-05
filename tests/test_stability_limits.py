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
