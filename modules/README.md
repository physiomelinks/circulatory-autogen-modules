# The module library

This directory is the module library that libcuflynx, PhLynx and CUFLynx load. It is organised by
**module_type**, each with its **versions** (the math) and each version's **instances** (parameter
sets). This file describes the layout, the naming and placement rules, how models use versions and
instances, and what is tested at each level. `directory_schema.json` is the machine-readable form
of the layout, and `tests/test_structure.py` checks the tree against it.

## Layout

```
modules/
  README.md                         this file
  directory_schema.json             the layout, machine-readable (tests/test_structure.py checks it)
  <category>/[<subcategory>/...]    physiological grouping only: subcategories and module_types, no files
    <module_type>/                  a directory with versions/ is a module_type
      <module_type>.html            report across the versions (generated)
      versions/
        <version>/                  version == the config entry's module_subtype
          <module_type>_<version>_modules.cellml         its CellML component
          <module_type>_<version>_modules_config.json    exactly one config entry, with "default_instance"
          <module_type>_<version>_units.cellml           the units it uses (and what they depend on)
          <module_type>_<version>_verification_config.json   what the checks run (see "The version spec")
          <module_type>_<version>_tests.yaml             review and record-keeping (see "The version spec")
          <module_type>_<version>_references.bib         BibTeX entries its parameters cite
          <module_type>_<version>_references_proposed.bib   citations proposed for review (optional)
          <module_type>_<version>.html                   version report (generated)
          risk/                                          joint failure-risk analysis (make risk; committed)
          plots/  results/                               written by the tests (not committed)
          instances/
            <instance>/                                  instance == its obs_data_name
              <instance>_parameters.csv                  variable_name,units,value,data_reference,sourced
              <instance>_obs_data.json                   calibration and held-out validation data (optional)
              <instance>_params_for_id.csv               parameters to calibrate, with bounds (optional)
              <instance>_calibrated_parameters.csv       written by calibration (committed)
              <instance>_calibration.json                calibration summary: method, cost, values, date (committed)
              SOURCES.md, raw data files                 where the data came from (e.g. hudson_bay_lynx_hare.csv)
```

The categories are:

| Category | Holds |
|---|---|
| `vessels/compartments`, `vessels/junctions`, `vessels/terminals`, `vessels/microvasculature`, `vessels/properties` | 0D vessel segments, junctions, terminal beds, microvascular networks, wall material laws |
| `cardiac` | the heart (monolithic and supermodule versions), cardiac clock, chamber, valve |
| `cell`, `cell/ion_channels`, `cell/neurons`, `cell/cardiomyocytes` | cell models and their parts |
| `respiratory` | lungs, gas exchange and transport |
| `control` | baroreflex, autonomic control, effectors, observers, PID control |
| `boundary_conditions` | inlet/outlet pressures and flows, generators, stimuli |
| `coupling` | 0D-1D coupling, volume sums |
| `transport` | tissue diffusion |
| `organs` | organ-level volume control |
| `benchmarks` | test problems (Lotka-Volterra, Van der Pol, FitzHugh-Nagumo, two minima, delay) |

`modules/poiseuille_transport/` and `modules/system/` are another session's work in the old layout and
are excluded from the schema (`directory_schema.json`, `excluded`) until they move.

## Placement rule

A module_type goes in **the most specific category that covers every context it is used in**, not as
shallow as possible:
- `soma`, `axon` and `varicosity` are only parts of neurons, so they are in `cell/neurons/`.
- A sarcoplasmic reticulum (or a generic Ca store) would belong to several cell types, so it would go
  in `cell/`, as `NKE_pump` does today.
- Every ion channel is in `cell/ion_channels/`, one module_type per channel kind (`i_Na`, `i_CaL`,
  ...), because channels of one kind are used across cell types. The cell type a version was built
  for is in its config entry's `"notes"`, e.g. `"Built for: sympathetic neuron (Tao et al. 2011)"`.
- `myocyte_membrane_voltage` is only used by the Paci cardiomyocyte, so it is in `cell/cardiomyocytes/`.

**No directory name appears twice** among the categories and module_types, so a category is never
named like a module_type: the category `cell/neurons/` holds the module_type `neuron`. The generic
names (`versions`, `instances`, `risk`, `plots`, `results`), version names (which repeat across
module_types: `nn`, `vp`, ...) and instance names (which name data sets, and `default`) are exempt.

## Versions

A version is one piece of math for a module_type: one CellML component and one config entry. **The
version's name is the config entry's `module_subtype`**, so a model picks a version with
`"module_subtype"`.
- **New versions** are named `<source>_vXX`: `Argus2026_v01`, `Paci2013_v01`, `Tao2011_v01`.
- **Existing subtypes** (`vp`, `pp_nonlinear`, `nn_constant`, ...) are kept as version names, except
  the ion channels, which became `Paci2013_v01`, `Tao2011_v01`, `Argus2026_v01` (the former `_SN`
  channels) and `Argus_unpublished_v01` (Piezo).
- **The port-letter constraint.** libcuflynx's vessel-port connection rule reads a module_subtype's
  first two letters as the inlet and outlet port kinds when they are `pp`, `pv`, `vp` or `vv`, and
  exempts names starting `nn`. So a version name starts with one of those four only when it follows
  that port convention. `Argus2026_v01` or `sympathetic` are safe; a version called `vp_Smith2027_v01`
  would be read as a v-in/p-out vessel.
- **Shared components.** libcuflynx finds CellML components by name across the whole library, so two
  versions can't define components with the same name. A component that several old entries shared
  was copied into each version and renamed `<component>__<module_type>_<version>` (in the CellML and
  in the config's `component_type`). The generated model names every component after its vessel, so
  results are unchanged.

## Instances

An instance is **a parameter set of a version: only the parameters change, never the math**.
- **Parameters.** `<instance>_parameters.csv` has one row per parameter, without a vessel suffix:
  libcuflynx renames a row `C` to `C_<vessel>` for the record that uses it. A row naming one of the
  version's `global_constant` variables keeps its plain name.
- **Naming.** An instance is named by its data set: its obs_data carries
  `"obs_data_name": "<instance>"`, which must equal the directory's name.
- **One data file.** `<instance>_obs_data.json` holds all of an instance's data. Its `data_items`
  are what calibration fits. Held-out validation data goes in the same file, as `prediction_items`
  that carry a `value` (with `data_type`, `std`, and `obs_dt` for a series), named
  `<variable>_validation`; the protocol runs to the end of the held-out data. There is no separate
  validation file. Each held-out series also has scalar held-out features: constant prediction_items
  with an `operation` (`<variable>_max_validation`, `_min_`, `_mean_`), their value and std taken
  from the data. They are validated like the series, and they are the prediction items that
  libcuflynx's sensitivity analysis and emulation can include as features
  (`include_prediction_items` in `sa_options` / `emulator_settings`; only items with an operation,
  i.e. a scalar, are used). libcuflynx validates the calibrated model against these items
  (`validation_results.json`), and CUFLynx shows the result.
  `python -m cam_testing.calibrate from-csv ... --window 0 15 --validation-window 0 20` makes one
  (with the max / min / mean features; `--validation-features` chooses them).
- **The default instance.** Every version has `default` (named in the config's `"default_instance"`):
  the values it has unless a model names another instance. The verification tests run at the
  default instance's parameters.
- **Data instances.** They currently hold the version's default values, plus their data, e.g.
  `benchmarks/Lotka_Volterra/versions/nn/instances/hudson_bay_lynx_hare/`. Calibration writes
  `<instance>_calibrated_parameters.csv` and `<instance>_calibration.json` beside them.

## How system models use versions and instances

System models live in `../system_models/<category>/<model>/` (not in `modules/`; they are not
modules or supermodules). A vessel-array record names a module_type, a version and an instance:

```json
{"name": "aortic_root", "module_type": "arterial_simple", "module_subtype": "vv", "instance": "default",
 "inp_instances": ["heart"], "out_instances": ["systemic_T", "volume_sum"]}
```

libcuflynx loads `<config dir>/instances/<instance>/<instance>_parameters.csv` for the record. With no
`"instance"`, it uses the config's `default_instance` if that file exists. **The model's own
`<model>_parameters.csv` wins** over instance values, so a model only lists what it changes. Every
record in `system_models/` names `"instance": "default"`, and the version specs' harness networks do
the same.

## Supermodule versions

A supermodule version is a version whose config entry has `"module_format": "supermodule"`. It is
built from other versions (its `submodules`, each with its own `module_type`, `module_subtype` and
`instance`), with internal connections only. It replaces the monolithic version of the same
module_type:

| Supermodule version | Built from |
|---|---|
| `cardiac/heart` `Argus2026_v01` | cardiac clock, four chambers, four valves |
| `cell/neurons/neuron` `sympathetic` | soma `sympathetic`, axon `sympathetic_monolithic_v01`, varicosity `sympathetic` |
| `cell/neurons/soma` `sympathetic` | membrane, Na/K, Ca handling and 15 `cell/ion_channels` versions |
| `cell/neurons/varicosity` `sympathetic` | membrane, Ca handling, channels, NE release |

The monolithic versions sit beside them. For example, `soma` also has `sympathetic_monolithic_v01`.

A model uses one record, e.g.
`{"name": "heart", "module_type": "heart", "module_subtype": "Argus2026_v01", "instance": "default",
"per_submodule_inputs": {"ra": ["venous_svc"]}, "per_submodule_outputs": {"aov": ["aortic_root"]}}`.
- **Expansion.** libcuflynx expands the record: each submodule becomes `<name>_<submodule>`, and the
  supermodule instance's rows `{var}_{submodule}` are renamed to match.
- **Precedence.** The model's parameters file wins, then the supermodule's instance, then each
  submodule's instance.
- **Where it lives.** Its `tests.yaml` has a `description`; its `verification_config.json` has a
  `supermodule` block. `globals` lists its global parameters. `equivalent` lists the system models
  it must reproduce: `system_models/closed_loop_cvs/3compartment_supermodules` reproduces
  `3compartment`, and `system_models/cellular/SN_simple_supermodules` reproduces `SN_simple`.

## What runs

**Per version** (`tests/test_modules.py`, ids `<module_type>/<version>`), at the default instance's
parameters:

| Test | Checks |
|---|---|
| `run_test` | generates the version alone (or in its spec's `harness` network) and simulates it; every output finite |
| `verification_test_invariants` | the spec's invariants hold at the run parameters |
| `verification_test_BC` | each boundary condition / constant swept over its tested range; runs finite, invariants hold. Points outside a valid range are skipped; failures are classified and a solver setting recommended (below) |
| `verification_test_timestep` | fixed-step integration converges at the scheme's order and matches CVODE |
| `stability_test` | the solver/tolerance/step matrix; the declared-supported configurations work |
| `phlynx_export_test`, `cuflynx_simulate_test`, `phlynx_equivalence_test` | the PhLynx → CUFLynx pipeline (`tests/test_phlynx.py`) |
| `supermodule_structure_test`, `supermodule_equivalence_test` | supermodule versions only, in place of the above |

**Per instance** (ids `<module_type>/<version>/<instance>`), for whatever data the instance has:

| Test | Checks |
|---|---|
| `validation_test_baseline` | the model at the instance's parameters against its baseline data or scalar targets |
| `validation_test_calibrate` | calibrates to `<instance>_obs_data.json` with libcuflynx, predicts the held-out data (the obs_data's `prediction_items` that carry a value), and writes the calibrated files |

How the boundary-condition sweep treats impossible values and failures:
- **Valid ranges.** A sweep value outside its parameter's valid range is not run; it is listed as
  "skipped: outside the valid range". The valid range is `bc_sweep.bounds[param]` when given. Without
  one, two generic rules apply: a dimensionless `<x>_init` whose state `x` an invariant keeps in
  [0, 1] (the invariant has both `x >= 0` or `x >= -eps` and `x <= 1` or `x <= 1 + eps`, eps up to
  1e-3) is a gate initial value and lies in [0, 1]; a concentration, amount, volume, conductance,
  capacitance, resistance, time or temperature (by its units) with a non-negative value stays >= 0.
  A point that breaks a `bc_sweep.constraints` relation (e.g. `B_Ca_init < B_Ca_total`) is skipped too.
  The report's tested range and the joint risk analysis use the clipped range.
- **Failure classes.** A failed run is rerun once with a tight reference: the reference tolerances
  (`reference_solver_info`, default rtol 1e-10) with MaximumStep dt/10. If the reference fails too it
  is a **parameter** failure (the model breaks there), and the test fails. If it passes it is a
  **numerical** failure (the version's `solver_info` is inadequate there). An invariant violation's size
  relative to `atol + rtol*|bound|` is recorded as a hint.
- **Recommended settings.** For numerical failures the test searches CVODE settings: the configurations
  the stability test found working plus a ladder of rtol 1e-6/1e-8/1e-10 by MaximumStep none/1e-4/1e-5,
  keeping only those at least as strict as the version's and accurate at the nominal point. They are
  ordered by measured run time, and a binary search finds the cheapest that passes the hard points (the
  numerical failures, the nominal point and each parameter's extremes, up to 10) within the stability
  test's `tol` of the reference. The winner is confirmed over every sweep point. Runs predicted over
  20 s per state variable are not made. The test passes when there are no parameter failures and every
  numerical failure passes with the recommended settings; the report shows the recommendation (settings,
  mean run time, failures n/N), which you may adopt into the version's `solver_info`.

Status rules (decided in the Lotka_Volterra review, 2026-10-01):
- An instance without obs_data has calibration **not applicable** ("no calibration data in this
  instance"); likewise an instance without baseline data has the baseline not applicable. The instance
  stays usable in models.
- The version's **Calibration** column in the report's test overview fails when **no instance** has
  calibration data, fails when any instance's calibration fails, and passes when they all pass.
- A check blocked by an external release (e.g. the released CUFLynx is too old for a feature) is a
  **failure with a known issue** (`expected_failures`), not not-applicable.
- In the stability matrix a configuration that is too inaccurate at a step (e.g. explicit Euler at a
  coarse dt) is shown as not working at that step; the test fails only if a declared-supported
  configuration fails.

The version's spec is described in the next section.

```bash
pytest tests/test_modules.py --module Lotka_Volterra          # a module_type
pytest tests/test_modules.py --module cell/neurons            # every module_type in a category
pytest tests/test_modules.py --component Lotka_Volterra/nn    # one version (and its instances)
pytest tests/test_modules.py --module Lotka_Volterra -m slow  # calibration
python -m cam_testing.report --module Lotka_Volterra          # Lotka_Volterra.html + the version pages
```

## The version spec: tests.yaml and verification_config.json

A version's spec is two files, merged into one spec when the tests run (the verification config
wins on the keys it owns):
- `<module_type>_<version>_verification_config.json` holds **what the checks run**. It is
  machine-read, pretty-printed JSON (indent 2) with its keys in the order of the table.
- `<module_type>_<version>_tests.yaml` holds **the review and record-keeping**, written by people.

Both files carry `module_type` and `version`, equal to the directory names. Any other key belongs in
exactly one file. A key in the wrong file, or a key in neither list, fails
`tests/test_structure.py` (the lists are `spec_keys` in `directory_schema.json`).

| Key | File | Meaning |
|---|---|---|
| `module_type`, `version` | both | the version's identity (the directory names) |
| `time_label` | verification_config | axis label for the time in plots (default `time [s]`) |
| `sim_time`, `pre_time`, `dt` | verification_config | simulated time, time run before logging, output interval [s] |
| `solver` | verification_config | libcuflynx solver for the runs (`CVODE_myokit`) |
| `solver_info` | verification_config | solver settings for the runs (`rtol`, `atol`, `MaximumStep`, ...) |
| `reference_solver_info` | verification_config | tolerances of the tight CVODE reference (timestep and stability tests) |
| `outputs` | verification_config | variables to log and check (`var`, or `vessel/var` for a harness neighbour); default: every `variable` of the config |
| `run_parameters` | verification_config | `{parameter: value}`: the operating point of the run and invariant tests |
| `harness` | verification_config | the test network: `vessel_array` rows `[name, module_subtype, module_type, inp, out, instance]` (the version under test is `mod`) and `parameters` rows `[name, units, value, source]` for the neighbours |
| `invariants` | verification_config | numpy expressions (`expr`, with `description` and `applies`) that must hold |
| `bc_sweep` | verification_config | the parameter sweep: `sweep` (`all` or `bcs`), `factors`, `points`, `ranges` (`[min, max]` or `{min, max, points, scale}` or `{values}`), `bounds` (`{param: [min, max]}`, `null` for no limit: the valid range; values outside it are skipped), `constraints` (numpy expressions of the parameters by variable name, e.g. `B_Ca_init < B_Ca_total`; a point that breaks one is skipped), `exclude`, `extra_parameters`, `plot_output`, `rationale` |
| `timestep` | verification_config | the convergence test: `scheme`, `dts`, `t_end`, `min_order`, `tol`, `cvode_tol`, `roundoff`, `wrapped` |
| `stability` | verification_config | the solver matrix: `supported` (must work), `cvode`, `solve_ivp`, `fixed_step`, `max_step_start`, `min_step`, `time_budget`, `t_end`, `tol` |
| `parameter_ranges` | verification_config | published parameter intervals, `{parameter: [lo, hi]}` (reported as validated values; usually inside a baseline block) |
| `validation` | verification_config | per instance, `validation.<instance>.baseline` / `.calibrate`: `status`, `source`, data references relative to the version directory (`data: instances/<i>/<file>.csv`, `obs_data`, `params_for_id`), variables, targets, metric, threshold, optimiser settings |
| `supermodule` | verification_config | supermodule versions: `globals` (parameters that name no submodule) and `equivalent` (system models it must reproduce: `model`, `reproduces`, `instance`, `tol`, `solver_info`, `output_map`, `ignore`) |
| `reviewed` | tests | whether the version has been reviewed (unreviewed versions run in CI without blocking it) |
| `description` | tests | what the version is (supermodule versions) |
| `notes` | tests | notes on the version, shown in its report |
| `skip` | tests | a reason not to run the version's tests |
| `known_issues` | tests | substrings of structure problems that are known (xfail) |
| `expected_failures` | tests | `{test: reason}`: known failures, recorded as failed and xfail in pytest |
| `reference_proposals` | tests | proposed parameter references (and values for TODO parameters), awaiting review |
| `review` | tests | the review: `summary`, `comment`, `questions`, `proposed_fixes`, `findings`. A review shared by several versions (a whole old module's review) is kept once, in `reviews/<name>_review.yaml` at the repo root, and `review:` gives that path instead |
| `review_scope` | tests | which versions a review copied from a whole old module covers |

## CUFLynx archives (.omex) per instance

Every instance can be opened in CUFLynx as one file: `instances/<instance>/<module_type>_<version>_<instance>.omex`, a COMBINE archive.

The archive is **generated** from the files above by `make omex` (`tools/build_instance_omex.py`; `cam_testing/omex.py`). It is never edited, and it isn't committed (it's gitignored). So there is only ever one place to change anything: the version's CellML, config, units or verification config, or the instance's parameters and data. Rebuild with `make omex MODULE=<module_type or category>`. CI builds the archives, and the version report links each instance's archive for download.

What an archive holds, in this order:

| Member | What it is |
|---|---|
| `<module_type>_<version>_<instance>.cellml` | **Master model** (marked in `manifest.xml`): the flattened model of the version's test network (the version alone, or its `harness` network) with this instance's parameters, as libcuflynx generates it. This is what CUFLynx opens and simulates. |
| `<instance>_obs_data.json` | The instance's data, if it has any (CUFLynx's observations): calibration data, and held-out validation data as `prediction_items`. |
| `<instance>_params_for_id.csv` | The parameters to identify, if any (CUFLynx's calibration setup). |
| `<instance>_parameters.csv`, raw data files, `SOURCES.md`, `<instance>_calibrated_parameters.csv`, `<instance>_calibration.json` | The instance's own files, unchanged. |
| `<module_type>_<version>_modules.cellml`, `_modules_config.json`, `_units.cellml` | The version's math, unchanged. |
| `<module_type>_<version>_verification_config.json` | The version's verification settings, in preparation for CUFLynx running them. |
| `<module_type>_<version>_<instance>_vessel_array.json`, `_model_parameters.csv` | The test network and the parameters file the model was generated from. |

CUFLynx picks members by name: the first `.json` with "obs" in its name is the obs_data, and the first `.csv` with "param" in its name is the params_for_id. So the study members come first.

An instance with no obs_data or params_for_id still loads and simulates, but CUFLynx then tries other members in those roles and shows two harmless "could not read" notes. That limitation is noted for a future change.

`tests/test_instance_omex.py` (`make omex-test`) loads every archive into a released CUFLynx and checks four things:
- CUFLynx takes the master model, the obs_data and the params_for_id;
- the run is finite;
- the outputs match libcuflynx's model of the same instance within 1e-6;
- the version report shows this as the instance's CUFLynx tick.

## Directory schema

`directory_schema.json` describes every level: `modules_root`, `category`, `module_type`, `versions`,
`version`, `instances`, `instance`, `risk`, and `system_models_root` / `system_category` /
`system_model`. For each level it gives the name pattern, the required and optional files (with
`{module_type}`, `{version}`, `{instance}`, `{model}` placeholders) and the allowed subdirectories.
`spec_keys` lists which spec keys go in which file. Its `rules` list the cross-level rules:
- the version is the module_subtype, with one config entry;
- the spec keys each sit in their own file (tests.yaml or verification_config.json);
- the default instance exists;
- an instance is its obs_data_name;
- no directory name repeats;
- no category is named like a module_type;
- where a module_type goes (placement);
- how versions are named;
- every system-model record names an existing version and instance.

`tests/test_structure.py` walks `modules/` and `system_models/` and checks all of this. It also checks
each version's CellML, config and units, and validates the JSON files against libcuflynx's schemas.

## Adding things

- **A new version** of an existing module_type: add
  `versions/<source>_vXX/` with the six files, one config entry (`module_subtype` = the directory name,
  `"default_instance": "default"`), a component name no other version uses, and
  `instances/default/default_parameters.csv`. Then run `make structure` and
  `pytest tests/test_modules.py --component <module_type>/<source>_vXX`.
- **A data instance**: add `instances/<name>/` with `<name>_parameters.csv`, the obs_data files
  (`"obs_data_name": "<name>"`), `<name>_params_for_id.csv`, the raw data and `SOURCES.md`. Then add
  `validation.<name>` to the version's verification_config.json.
- **A new module_type**: put it in the most specific category covering all its uses (create the
  category if none fits, with a name no module_type has). Then run `make manifests`.

`tools/restructure_modules.py` moved the old per-module layout into this one, driven by
`tools/restructure_map.yaml` (made by `tools/restructure_map.py`).
