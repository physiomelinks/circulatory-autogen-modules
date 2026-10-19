"""
cam_testing's test suite, for this library and any repo laid out the same way (modules/README.md).

    test_structure       static checks: directory schema, config <-> CellML, units, uniqueness, manifests
    test_modules         per-version verification and per-instance validation (libcuflynx)
    test_systems         system models (system_models/): run, reproduce the reference, invariants
    test_phlynx          the PhLynx -> .omex -> CUFLynx pipeline, per version
    test_instance_omex   each instance's COMBINE archive in CUFLynx

Run them with ``pytest --pyargs cam_testing.suite`` (or one module of it), or from a test file of
your own: ``from cam_testing.suite.test_modules import *``. The plugin cam_testing.pytest_plugin
(registered on install) provides the options, parametrisation and fixtures; cam_testing.paths says
which repo is tested.
"""
