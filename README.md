# Circulatory Autogen Modules

The module library for [circulatory_autogen / libcuflynx](https://github.com/physiomelinks/circulatory_autogen) and Phlynx. Every module version is verified, and every instance validated, with libcuflynx; each module_type and each version gets an HTML report covering its equations, parameters and test results.

Phlynx loads the versions listed in `manifests/index.json` (and `vitalworkshop.json`). `manifests/all.json` lists every version. Rebuild the manifests with `make manifests` after adding a version.

## Layout

**[`modules/README.md`](modules/README.md) documents the library in full**: the layout, the placement and naming rules, versions and instances, supermodule versions, what is tested per version and per instance, and the directory schema (`modules/directory_schema.json`, checked by `tests/test_structure.py`). In short:

```
modules/<category>/<module_type>/                    (or modules/<module_type>/, e.g. heart)
  <module_type>.html                                 report across the versions (generated)
  <nested module_type>/                              module_types used only within this one (same layout)
  versions/<version>/                                version == the config entry's module_subtype
    <module_type>_<version>_modules.cellml / _modules_config.json / _units.cellml
    <module_type>_<version>_verification_config.json   what the checks run (JSON)
    <module_type>_<version>_tests.yaml               review and record-keeping
    <module_type>_<version>_references.bib           BibTeX entries the parameters cite
    risk/                                            joint failure-risk samples + classifier (committed)
    instances/<instance>/                            parameter sets; instance == obs_data_name
      <instance>_parameters.csv                      variable_name,units,value,data_reference,sourced
      <instance>_obs_data.json, _params_for_id.csv, ...  calibration data and results
system_models/<category>/<model>/                    whole models built from the library
```

To use the library from libcuflynx, set these in `user_inputs.yaml` (on circulatory_autogen master; not yet in a release):

```yaml
module_library_dirs: [/path/to/circulatory-autogen-modules/modules]
use_builtin_modules: false
```

## Reports

`python -m cam_testing.report` (`make report` / `make site`, also run by CI) builds the HTML at two levels from repo files alone: `<module_type>/<module_type>.html` lists every version with its test results and instances, and `versions/<v>/<module_type>_<v>.html` is the version's page:
- equations from the CellML;
- ports and variables from the config;
- values, references, sourced flags and ranges from the default instance's parameters;
- the bibliography from `<module_type>_<v>_references.bib`;
- invariants and validation settings from `<module_type>_<v>_verification_config.json`, notes and review from `<module_type>_<v>_tests.yaml`;
- an instance table with each instance's data, baseline and calibration results;
- results and plots from the last test run.

The site index (`make site`) groups the module_types by category (heart, directly under `modules/`, is its own group), lists a nested module_type by its path (e.g. `neuron/soma/SN_membrane_soma` under `cell`) and links both levels.

Nothing else is needed to regenerate them.

## Tests

Each version (one config entry: `module_type`/`module_subtype`) is generated alone with libcuflynx at its default instance's parameters. Every boundary condition becomes a parameter. The version then goes through the verification tests; each instance goes through the validation tests, for whatever data it has (calibration is not applicable to an instance without obs_data; the version's Calibration in the report overview fails when no instance has calibration data; see `modules/README.md`):

| Test | What it checks |
|---|---|
| `run_test` | generates and simulates the version; every output finite |
| `verification_test_invariants` | the spec's invariants hold (exact solutions, conservation laws, bounds, delays) |
| `verification_test_BC` | sweeps each boundary condition and constant over a range (a supermodule: its globals and open boundary conditions); every run finite and invariants hold |
| `verification_test_timestep` | fixed-step integration of the generated right-hand side at halved steps converges at the scheme's order, and agrees with CVODE at tight tolerances |
| `stability_test` | a matrix of solvers, tolerances and timesteps (CVODE, SciPy `solve_ivp`, fixed-step Euler/Heun/RK4); reports which work, and requires the declared-supported ones to work |
| `validation_test_baseline` | per instance: model vs baseline data (NRMSE threshold) or scalar targets |
| `validation_test_calibrate` | per instance: calibrate to its obs_data with libcuflynx parameter identification, predict held-out data, and write `<instance>_calibrated_parameters.csv` / `<instance>_calibration.json` |

Supermodule versions run the same tests, generated alone and expanded into their submodules, plus `supermodule_structure_test` and `supermodule_equivalence_test` (the system model using them reproduces the one with the submodules written out). Their BC sweep varies only their global constants and the boundary conditions no internal connection closes (`bc_sweep.sweep: globals_and_bcs`, the supermodule default; `all` sweeps every parameter): each submodule's own parameters are swept in that submodule version's tests. PhLynx has no supermodule support, so their PhLynx export test fails as a recorded known issue.

Every version's report has the same test columns (Run, Invariants, BC sweep, Timestep, Stability, Calibration, PhLynx export, CUFLynx simulate, PhLynx equivalence, Supermodule structure, Supermodule reproduces). A column that doesn't apply shows N/A with the reason, one that applies but hasn't run shows Not run, and only a test that ran and passed shows Passed. A version with no calibration data of its own used inside supermodule versions shows "Pass in super" / "Fail in super" for Calibration, from those supermodules.

`tests/test_structure.py` holds fast static checks: the directory schema, config ↔ CellML consistency, units, library-wide uniqueness and manifest paths. Known defects are listed under `known_issues` in a version's spec. They show as xfail until fixed.

A version whose spec has `reviewed: false` still runs in CI, but its failures don't fail the build.

## Parameter ranges

Each parameter's **tested range** is set in `<module_type>_<version>_verification_config.json`:
- `bc_sweep.ranges.<param>` as `[min, max]`, or `{min, max, points, scale: log}`;
- otherwise `bc_sweep.factors` × the nominal value.

`verification_test_BC` varies each parameter over its tested range with the others at nominal, and fails if any run fails. `make risk` samples all parameters together over the same ranges. The report shows the tested range and the **validated values** for each parameter; the validated values are the spread of values that matched data in passing validation tests.

To check a system model's parameters against the modules it uses:

```bash
python -m cam_testing.ranges check --module-array my_module_array.json --parameters my_parameters.csv
```

This lists each value as inside or outside its tested range and validated spread (over all the version's instances), and estimates the failure risk near the values from the `make risk` samples. It exits non-zero if any value is outside a tested range.

## Failure risk over the joint parameter space

The BC sweep varies one parameter at a time, so it doesn't cover combinations. `make risk MODULE=<module_type or category>` samples every parameter together, using a scrambled Sobol sequence over the same box the sweep uses. Each sample counts as a failure if the run errors, an output is non-finite, or an invariant is violated.

The analysis reports, for each version:
- the overall failure probability with a 95% interval; if no sample fails, an upper bound of 3/N;
- a failure-risk curve for each parameter;
- each parameter's η²: the share of the failure variance it explains on its own, which ranks what drives failures;
- a breakdown by failure mode;
- a corner plot of pairwise failure rates.

The results are committed in each version's `risk/` (`<module_type>_<version>_risk.json` and `_risk_samples.npz`, tens of KB each), so they're available on any clone. The JSON includes a failure classifier: logistic regression on quadratic terms of the box coordinates, stored as plain coefficients, with its 5-fold cross-validated AUC. It gives P(fail) at any parameter set in the box.

The version report shows all of this, regenerating the plots from the committed samples. The site index lists each module_type's highest version risk. `ranges check` also estimates the failure risk near a system model's parameter values from the stored samples. On GitHub, run the workflow manually with "risk" ticked to include it in the published reports.

## System models

`system_models/<category>/<model>/` (at the repo root) holds circulatory_autogen's CellML models rebuilt from this library, with the heart split into cardiac clock, chamber and valve vessels (`heart/cardiac_clock`, `heart/chamber`, `heart/valve`). System models are not modules or supermodules: they use the library. Every module-array record names a module_type, a version (`module_subtype`) and an instance (`"instance": "default"`); the model's own parameters file wins over the instance values. Each directory contains:
- the module array and parameters;
- the original under `reference/`: circulatory_autogen's files with only the bib keys changed to this library's (`tools/bib_rekey_map.json`, applied by `tools/import_systems.py`; checked against the originals by `tests/test_structure.py`);
- a `<model>_system.yaml` spec.

The reference can instead be a ready CellML model, simulated as it is: for example the flattened model PhLynx exports, for a model first built in PhLynx. Set `equivalence.reference_cellml: reference/<file>.cellml` in the spec (a path relative to the model's directory), or put it at `reference/<model>.cellml` with no vessel array beside it. Its outputs are `<component>/<variable>` for every component except `environment`, `global_parameters` and `instance_parameters` (`equivalence.ignore_components` changes the list); `output_map`, `ignore`, `wrapped`, `reciprocal` and `tol` work as for an original. PhLynx names a vessel's component by the vessel's name, so its outputs match the library model's `<vessel>/<variable>` with no `output_map`.

`tests/test_systems.py` checks three things for each model:
- it runs;
- it reproduces the original, with every logged output within 1e-6 at the spec's CVODE tolerances;
- its invariants hold.

Models whose original doesn't generate in circulatory_autogen carry `expected_failures` or `skip` with the reason. `tools/import_systems.py --ca-dir ../circulatory_autogen` re-imports them and never overwrites a spec.

Some system models are examples built from this library, with no circulatory_autogen original (equivalence `not_applicable`); their invariants check exact solutions instead:
- `system_models/diffusion/`: finite-volume diffusion meshes (`tools/build_diffusion_examples.py`).

Module arrays are JSON (`<model>_module_array.json`): a list of records in PhLynx's key names, e.g. `{"name": "venous_svc", "module_type": "venous", "module_subtype": "vp", "instance": "default", "inp_instances": ["systemic_T"], "out_instances": ["heart", "volume_sum"]}`. libcuflynx reads them, and still reads the older CSV layout by converting each row to the same record. `tools/convert_module_arrays.py` converts CSV files; `cam_testing/module_array.py` reads either.

Supermodules are versions: `heart` version `Argus2026_v01` is the cardiac clock, four chambers and four valves, and `cell/neuron`, its nested `soma` and `varicosity` have `sympathetic` supermodule versions (see [`modules/README.md`](modules/README.md#supermodule-versions)). A model uses one through one record:

```json
{"name": "heart", "module_type": "heart", "module_subtype": "Argus2026_v01", "instance": "default",
 "per_submodule_inputs":  {"ra": ["venous_svc"], "la": ["pvn"]},
 "per_submodule_outputs": {"puv": ["par"], "aov": ["aortic_root"], "ra": ["volume_sum"], "rv": ["volume_sum"], "la": ["volume_sum"], "lv": ["volume_sum"]}}
```

The host's own vessels name `heart` in their `inp_instances` / `out_instances`. libcuflynx expands the record: each submodule becomes `heart_<submodule>` (`heart_ra`, `heart_aov`, …), the listed host modules are coupled to it, and the instance's parameters are renamed to match (`E_A_ra` becomes `E_A_heart_ra`), with the model's own values taking precedence. `system_models/closed_loop_cvs/3compartment_supermodules` is 3compartment built this way; its version-level test (`pytest tests/test_modules.py -k supermodule`, also `make systems`) checks that it reproduces `3compartment` to 1e-9, and `cellular/SN_simple_supermodules` reproduces `SN_simple`.

## PhLynx → CUFLynx pipeline

`tests/test_phlynx.py` checks that each version works end to end in the web tooling. Its three tests:
- `phlynx_export_test`: builds the version's test network in PhLynx, using PhLynx's own code from a checkout, run headlessly under Node and jsdom by `cam_testing/bridges/phlynx/export_omex.mjs`. It loads this library's modules and parameters, checks every connection was made, and exports the `.omex` PhLynx sends to CUFLynx.
- `cuflynx_simulate_test`: imports that archive into a released CUFLynx binary through its HTTP API (`cam_testing/bridges/cuflynx/simulate_omex.py`) and simulates it.
- `phlynx_equivalence_test`: compares CUFLynx's run with libcuflynx's model of the same network, within 1e-6 at the output times both runs share.

A version passing all three gets the "PhLynx / CUFLynx compatible" tick in its report, on its module_type's page and in the site index.

```
make pipeline-setup                                  # jsdom for the bridge
make pipeline MODULE=coupling                        # a category or module_type; PHLYNX_DIR=../phlynx, CUFLYNX_BIN=~/software/CUFLynx
make pipeline PHLYNX_DIR=/path/to/phlynx CUFLYNX_BIN=/path/to/CUFLynx
```

The tests skip when Node (22.15 or newer), the PhLynx checkout (with `yarn install` done) or the CUFLynx binary isn't available. CI checks out PhLynx at `PHLYNX_REF` and downloads the latest CUFLynx Ubuntu release.

Parameters are applied with PhLynx's own `applyParametersToNodes` (`src/utils/parameters.js`, from PhLynx #595). For an older PhLynx checkout without it, the bridge falls back to a copy of `loadParametersData`.

## Using cam_testing in another repository

`cam_testing` tests, reports and builds manifests for any repo laid out like this one: a `modules/` directory with the same layout and files (`modules/<category>/<module_type>/versions/<version>/...`, see [`modules/README.md`](modules/README.md)), and optionally `system_models/`, `manifests/` and `reviews/` beside it. That repo's modules (and system models) may use this library's modules, e.g. a harness downstream of `inlet_flow` and `arterial_simple`, by listing this library as an **extra module library**.

**Install** (libcuflynx needs SUNDIALS and MPI, see "Running locally"). The `libcuflynx` extra installs libcuflynx at the git ref `requirements.txt` pins:

```bash
pip install "cam-testing[libcuflynx] @ git+https://github.com/physiomelinks/circulatory-autogen-modules@<branch or tag>"
```

The package holds the code, the test suite, the PhLynx/CUFLynx bridges and a copy of the directory schema, not this library's modules. To use them, check out this repo too (a git submodule, or a clone at the same ref in CI) and point to its `modules/` as an extra library.

**Configure** in the repo's `pyproject.toml` (paths relative to that file):

```toml
[tool.cam_testing]
# root = "."                                      # the repo under test (default: this file's directory)
module_library_dirs = ["../circulatory-autogen-modules/modules"]   # a modules/ directory, or a repo holding one
# [tool.cam_testing.manifests]                   # curated PhLynx manifests besides all.json (optional)
# "index.json" = ["vessels", "control"]
```

or with environment variables, which win over `pyproject.toml`: `CAM_REPO_ROOT=/path/to/repo` and `CAM_MODULE_LIBRARY_DIRS=/path/to/circulatory-autogen-modules/modules` (several separated by `:`), or the pytest options `--cam-root` and `--cam-library` (repeatable), which win over both. With none of them, the repo is the nearest directory at or above the working directory that has `modules/`. `python -m cam_testing.paths` prints what is in use and where it came from; pytest prints it in its header.

With extra libraries:
- only the repo's own `modules/` is tested, reported and put in its manifests;
- libcuflynx gets `module_library_dirs: [<repo>/modules, <extra libraries>...]` for every model cam_testing generates (a version alone or in its harness, system models, the C++ run, the PhLynx pipeline), so harnesses, supermodules and system models can name versions from any of them. A module_type may be in both (e.g. your own `arterial_simple` versions): each version is found in the library that has it;
- the structure checks look across all of them: component names, `(module_type, version)` pairs and units must not clash with the extra libraries' (libcuflynx merges them), and system-model records and harnesses may name their versions.

The repo's own `modules/directory_schema.json` is used when it has one, else the copy in the package (`cam_testing/data/directory_schema.json`).

**Tests.** The plugin `cam_testing.pytest_plugin` is registered on install, with the options `--module`, `--component`, `--include-unreviewed`, `--quick-unreviewed`, `--cam-root` and `--cam-library` and the test names (`run_test`, `verification_test_BC`, ...). Run the suite straight from the package:

```bash
pytest --pyargs cam_testing.suite.test_structure                      # static checks
pytest --pyargs cam_testing.suite.test_modules -m "not slow" --module control   # V&V (MODULE as here)
pytest --pyargs cam_testing.suite.test_systems                        # system models
pytest --pyargs cam_testing.suite                                     # everything, incl. the PhLynx pipeline
```

or from one-line test files, which keep `pytest tests/...`, `-k` and `testpaths` working as for your own tests:

```python
# tests/test_cam_modules.py (likewise test_cam_structure.py, test_cam_systems.py, ...)
from cam_testing.suite.test_modules import *  # noqa: F401,F403
```

The suite's modules are `test_structure`, `test_modules`, `test_systems`, `test_phlynx` and `test_instance_omex`. A system model's reference can be a ready CellML model, e.g. the model PhLynx exports (see "System models").

**Reports and manifests** for the configured repo: `python -m cam_testing.report [--module X] [--site]` writes the pages into its `modules/` (and `site/`), `python -m cam_testing.manifests` writes its `manifests/*.json`. For the PhLynx/CUFLynx pipeline, the bridges ship in the package: install jsdom next to `cam_testing/bridges/phlynx` (`npm ci` there) or anywhere and set `JSDOM_DIR`, and set `PHLYNX_DIR` and `CUFLYNX_BIN` as here.

A minimal Makefile and CI job:

```make
PYTHON ?= venv/bin/python
# a checkout of circulatory-autogen-modules (or set module_library_dirs in pyproject.toml instead)
CAM_LIB ?= ../circulatory-autogen-modules
export CAM_MODULE_LIBRARY_DIRS = $(CAM_LIB)/modules
MODULE_ARGS = $(if $(MODULE),--module $(MODULE),)
structure: ; $(PYTHON) -m pytest --pyargs cam_testing.suite.test_structure
test:      ; $(PYTHON) -m pytest --pyargs cam_testing.suite.test_modules -m "not slow" $(MODULE_ARGS)
systems:   ; $(PYTHON) -m pytest --pyargs cam_testing.suite.test_systems
report:    ; $(PYTHON) -m cam_testing.report $(MODULE_ARGS)
manifests: ; $(PYTHON) -m cam_testing.manifests
```

```yaml
    steps:
      - uses: actions/checkout@v4
      - uses: actions/checkout@v4
        with: { repository: physiomelinks/circulatory-autogen-modules, ref: <branch or tag>, path: _cam }
      - uses: actions/setup-python@v5
        with: { python-version: '3.11' }
      - run: sudo apt-get install -y -qq libopenmpi-dev openmpi-bin libsundials-dev build-essential
      - run: pip install "./_cam[libcuflynx]"     # cam_testing and libcuflynx at that checkout's ref
      - run: make structure test CAM_LIB=$PWD/_cam PYTHON=python
```

## Running locally

libcuflynx needs SUNDIALS and MPI. On Ubuntu:

```bash
sudo apt-get install libsundials-dev libopenmpi-dev openmpi-bin build-essential
```

Then:

```bash
make setup                                   # venv + libcuflynx (pinned git ref) + cam_testing
make setup LIBCUFLYNX=../circulatory_autogen # ...or a local libcuflynx checkout
make structure                               # static checks
make test MODULE=Lotka_Volterra              # V&V tests for one module_type and those nested in it (MODULE=cell: a category; make test: all)
make test-all MODULE=Lotka_Volterra          # including slow calibration
make systems                                 # system models + supermodule versions
make risk MODULE=Lotka_Volterra              # joint failure-risk analysis (on request)
make report MODULE=Lotka_Volterra            # -> modules/benchmarks/Lotka_Volterra/Lotka_Volterra.html + version pages
make serve                                   # build site/ as GitHub Pages serves it; http://localhost:8000
```

The reports open straight from disk. Plots load by relative path, and MathJax loads from a CDN.

## CI

`.github/workflows/module-tests.yml` runs the structural checks and then, per top-level directory of `modules/` (`vessels`, `heart`, `cell`, ...: a category, or a module_type directly under `modules/` with its nested module_types), one job per stage: `vv` (V&V tests), `pipeline` (PhLynx -> CUFLynx) and `omex` (CUFLynx archives). A stage with more work is split into shards of whole module_types, about ten minutes each (`ci/shards.json`; `cam_testing/shards.py` balances them by the per-module_type seconds in `ci/test_durations.json`, refreshed from the jobs' JUnit artifacts with `python -m cam_testing.shards durations <stage> junit-*.xml`; `python -m cam_testing.shards plan` shows the split). Locally, `PYTEST_ARGS="--shard 2/4 --shard-stage vv"` runs one shard. Each job uploads its results, and the site job builds every report from them; failures of versions not reviewed yet are warnings (`tools/ci_reviewed_failures.py`, over the V&V, PhLynx -> CUFLynx pipeline and CUFLynx archive tests). On `main`, the site is assembled and deployed to GitHub Pages (enable Pages with source "GitHub Actions" in the repository settings). The same `make` targets run locally.

## Importing from libcuflynx

`tools/import_from_libcuflynx.py --ca-dir ../circulatory_autogen` copies modules from a libcuflynx checkout into the pre-versions per-module layout and extracts each module's units; `tools/restructure_modules.py` then moves such a tree into versions and instances (driven by `tools/restructure_map.yaml`). It seeds parameter values from libcuflynx's example models. It never overwrites a module's test spec, and it only overwrites the parameters with `--reseed-parameters`.
