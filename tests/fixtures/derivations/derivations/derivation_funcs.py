"""Derivation functions of the derivation-checks fixture (tests/test_derivation_checks.py)."""


def death_rate(inputs, units):
    """The predator death rate: a baseline rate plus the prey birth rate.

    Not a model of anything: a chained Python derivation for the tests.
    """
    return {'gamma': inputs['d0'] + inputs['alpha']}
