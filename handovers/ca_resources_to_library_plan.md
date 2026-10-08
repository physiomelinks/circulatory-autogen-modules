# Plan: circulatory_autogen's `resources/` replaced by circulatory-autogen-modules

**Decision (user):** "circulatory_autogen resources will be removed and corresponding models in circulatory-autogen-modules will be used. Tests will also call to the modules repo. Make sure this is planned for in the set of circulatory_autogen PRs."

**Snapshot:** 2026-10-07.
- CA `upstream/master` is `a36bfa29`. #521, #534, #537 and #549 are merged.
- Open CA PRs read: #520, #522, #535, #536, #538, #540, #541, #547, #548, #550, #552.
- Library at `devel_module_tests` `fbd9b70`. Library PRs #3, #4, #5 and #6 are open; #7 and #8 are merged.
- Read-only. Nothing was changed in either repo.

---

## 0. Summary

| | Count |
|---|---|
| Files in CA `resources/` (master) | **126**: 39 models with a module array, 2 external CellML models, 4 helper files (2 scripts, 2 txt) |
| CA resource models with a library counterpart in `system_models/` | **35 of 39** (same names; none renamed) |
| ...of which generate and run in the library | **24** (4 of them with an xfailed equivalence test) |
| ...that don't generate (6) or are skipped (5) in the library | **11** |
| CA resource models with **no** library counterpart | **7**: `aortic_bif_1d`, `aortic_bif_hybrid_V1`, `aortic_bif_hybrid_V2`, `cvs_model_with_arm_hybrid` (C++/1D, excluded by `tools/import_systems.py`); `Goodwin`, `Teusink` (external CellML); plus `aortic_bif_0d`, added by #520 |
| Resource data files with no library copy | **29** files (§1.2) |
| Built-in module entries (`*_modules_config.json`) | **255** on master, 256 on #520/#538 |
| ...with no library counterpart (after rename maps) | **6** on master; **7** with #520's `FV1D_solver` (§1.3) |
| CA test files that read model inputs from `resources_dir` | **43** (52 mention `resources` at all) on the tip of the stack (`feat/tests-on-lumped-vessels`, #540) |
| CA test files that name a resource model | **69**, including tests that use committed pre-generated models |
| CA non-test files that read or document `resources/` | user_run_files 2, tutorial 13, benchmarks 4, `src/libcuflynx/scripts` 15, `external_testing` 1, README, CLAUDE.md, funcs_user READMEs 3 |

**Already done by open PRs.** #538 has already done most of the "tests call the modules repo" work:
- `CUFLYNX_MODULE_LIBRARY` is the default module source for configs that choose none.
- `tests/conftest.py` points the suite at the library, and CI checks the library out.
- CA-only test modules have moved to `tests/test_inputs/test_modules`.
- `utilities/junction_migration.py` exists, and `resources/` has been rewritten to library module names.

What is missing:
1. The **models** (`resources/`) still live in CA.
2. The CI library ref is a **branch** (`fix/consistent-junctions`), not a pinned commit.
3. That branch **predates the library's 2026-10-06 restructure**, so CA's renamed `resources/` uses names the library no longer has.
4. The built-in modules still ship in the wheel.
5. The library's own equivalence tests still **generate CA's originals with CA's built-in modules**. They break once the built-ins go (§2.4).

---

## 1. Inventory

### 1.1 CA resource models and their library counterparts

CA `resources/` holds 39 models with a module array (plus 2 BACKUP copies), 2 external CellML models (`Goodwin.cellml`, `Teusink.cellml`), and 4 helpers (`excel_to_obs_data_json.py`, `parameters_cellml_2_csv.py`, `data_references.txt`, `modifications.txt`). #520 adds `aortic_bif_0d` and edits the hybrid models.

Library status comes from `system_models/*/*/*_system.yaml` (`skip`, `expected_failures`). No model was renamed: each library model keeps CA's name.

| CA model | Library | Status in library |
|---|---|---|
| 3compartment | `closed_loop_cvs/3compartment` | generates; equivalent |
| 3compartment_extra_ops | `closed_loop_cvs/3compartment_extra_ops` | generates; equivalent |
| 3compartment_nonstiff | `closed_loop_cvs/3compartment_nonstiff` | generates; equivalent |
| simple_physiological | `closed_loop_cvs/simple_physiological` | generates; equivalent |
| physiological | `closed_loop_cvs/physiological` | generates; equivalent |
| cvs_model_0d | `closed_loop_cvs/cvs_model_0d` | generates; equivalent |
| cvs_model_with_arm_0d | `closed_loop_cvs/cvs_model_with_arm_0d` | generates; equivalent |
| neonatal | `closed_loop_cvs/neonatal` | generates; equivalent |
| generic_junction_test_closed_loop | `closed_loop_cvs/...` | generates; equivalent |
| generic_junction_test_open_loop | `open_loop_arterial/...` | generates; equivalent |
| control_parasymp | `controlled_cvs/control_parasymp` | generates; equivalent |
| parasympathetic_model | `controlled_cvs/parasympathetic_model` | generates; equivalent |
| test_fft | `open_loop_arterial/test_fft` | generates; equivalent |
| NKE_pump | `cellular/NKE_pump` | generates; equivalent |
| Lotka_Volterra, VanDerPol, FitzHugh_Nagumo, TwoMinima, delay_test, pid_control | `benchmarks/*` | generate; equivalent |
| control_phys | `controlled_cvs/control_phys` | generates; equivalence xfail (valve internals 5e-6 > 1e-6) |
| control_phys_asd | `controlled_cvs/control_phys_asd` | generates; equivalence xfail (same) |
| SN_simple | `cellular/SN_simple` | generates; equivalence xfail: the library's Nav1.6 and E_Ca fixes. CA #537 merged the same fixes on 2026-10-01, so `reference/` is probably stale and the xfail may now pass (re-import and check) |
| behdad_test1 | `open_loop_arterial/behdad_test1` | generates; equivalence xfail (the original fails in CA, #529) |
| FinalModel | `closed_loop_cvs/FinalModel` | **doesn't generate** (Delta_q_us_venous_* units; fails in CA too) |
| FTU_wCVS | `closed_loop_cvs/FTU_wCVS` | **doesn't generate** (EMPTY_MUST_BE_FILLED parameters; fails in CA too) |
| new_valve_p_est | `closed_loop_cvs/new_valve_p_est` | **doesn't generate** (same) |
| elic | `controlled_cvs/elic` | **doesn't generate** (`pressure_observer_terminal` isn't defined anywhere) |
| cerebral | `open_loop_arterial/cerebral` | **doesn't generate** (terminal2 `u_in` port variable undeclared; also in CA) |
| cerebral_elic | `open_loop_arterial/cerebral_elic` | **doesn't generate** (same) |
| generic_junction_test2_closed_loop | `closed_loop_cvs/...` | **skipped**: no parameters file. CA's tests reuse `generic_junction_test_closed_loop_parameters.csv` |
| generic_junction_test2_open_loop | `open_loop_arterial/...` | **skipped**: same, with the open-loop file |
| ports_test | `benchmarks/ports_test` | **skipped**: uses CA test modules |
| test_init_states | `benchmarks/test_init_states` | **skipped**: uses CA test modules |
| Simple_ODE_Benchmark | `benchmarks/Simple_ODE_Benchmark` | **skipped**: uses CA test modules |
| aortic_bif_1d | **none** | C++/1D, excluded from the import |
| aortic_bif_hybrid_V1 | **none** | C++/1D, excluded |
| aortic_bif_hybrid_V2 | **none** | C++/1D, excluded |
| cvs_model_with_arm_hybrid | **none** | C++/1D, excluded |
| aortic_bif_0d (#520) | **none** | added by #520 |
| Goodwin | **none** | external CellML, not generated |
| Teusink | **none** | external CellML, not generated |

The library also has 5 models CA doesn't: `3compartment_supermodules`, `SN_simple_flat`, `SN_simple_supermodules`, `diffusion_2d_triangles` and `diffusion_3d_boxes`.

**The library's models are not CA's models under a new path.** `tools/import_systems.py` rebuilds every model that has a heart with the **split heart**: `cardiac_clock`, four `chamber`s and four `valve`s replace the monolithic `heart`. This applies to 17 models:
- 3compartment, 3compartment_extra_ops and 3compartment_nonstiff;
- simple_physiological, physiological, cvs_model_0d, cvs_model_with_arm_0d and neonatal;
- generic_junction_test_closed_loop and generic_junction_test2_closed_loop;
- control_phys, control_phys_asd, control_parasymp and elic;
- FinalModel, FTU_wCVS and new_valve_p_est.

The split heart changes the names that matter to tests:
- outputs: `heart/u_lv` becomes `lv/u`;
- params_for_id rows: `global,E_lv_A` becomes `lv,E_A`, and `q_lv_init` becomes `lv,q_init`.

15 CA test files hard-code such names, for example `test_param_id` (48 hits), `test_solvers` (14), `test_sensitivity_analysis` (11) and `test_fsa_grouped_sensitivities` (8). The library keeps CA's original of each model in `system_models/.../<model>/reference/`, under CA's old module names. This drives decision **D1**.

### 1.2 Resource data files with no library copy (29)

| File(s) | Used by | Proposed home |
|---|---|---|
| `3compartment_q_lv_only_params_for_id.csv` | test_param_id | CA `tests/test_inputs/` (a test fixture) |
| `Lotka_Volterra_multisub_obs_data.json`, `Lotka_Volterra_forced_multi_trace_obs_data.json` | test_param_id, test_solvers | CA `tests/test_inputs/` (next to the existing `Lotka_Volterra_forced_*` fixtures) |
| `control_phys_`, `physiological_`, `new_valve_p_est_prediction_variables.csv` | none in tests (physiological is named in `.gitignore`) | library, next to each model |
| `Goodwin.cellml` (+obs_data, params_for_id), `Teusink.cellml` (+obs_data, params_for_id), `modifications.txt` | `benchmarks/benchmark_specs.py`, test_benchmarks, **#547 `test_reporting`** | **D4**: library `system_models/benchmarks/` as ready-CellML models (needs library #6's `equivalence.reference_cellml`), or CA `benchmarks/data/` |
| `aortic_bif_1d`, `aortic_bif_hybrid_V1/V2`, `cvs_model_with_arm_hybrid` (module array + parameters), and #520's `aortic_bif_0d` | test_autogeneration (cpp), test_module_instances, test_cpp_template_generator, test_vessel_configurations_1d, test_lumped_vessels; `coupler/`, `solver1d/main1D.py`, the tutorial notebook | library `system_models/coupled_1d/` once library cam_testing has a C++ path (library #5 adds the C++ check). Until then, CA `tests/test_inputs/models/` (**D5**) |
| `data_references.txt` | test_issue_fixes (FTU_wCVS parsing) | CA `tests/test_inputs/` |
| `pulmonary_lobe_imped-{pre,post}PEA_obs_data.json` | nothing | library data, or delete (**D6**) |
| `cvs_model_with_arm_0d_*_BACKUP.csv`, `cvs_model_with_arm_extra_parameters.csv` | nothing | delete |
| `excel_to_obs_data_json.py`, `parameters_cellml_2_csv.py` | nothing | move to `src/libcuflynx/scripts/`, or delete |

### 1.3 Built-in modules (`src/libcuflynx/generators/resources/`)

There are 47 files: 22 `*_modules.cellml`, 23 configs, `base_script.cellml` and `units.cellml`. They are package data in the wheel (`pyproject` `"libcuflynx.generators" = ["resources/*.cellml", "resources/*.json", "*.cpp"]`). `resources/` itself is **not** in the wheel.

Each `(vessel_type, BC_type)` was mapped through `cam_testing.library.legacy_renames()` (176 pairs) and #538's `junction_migration.JUNCTION_TWINS`. Built-in entries **without** a library counterpart:

| Built-in | Notes |
|---|---|
| `FV1D_solver/nn` (#520, `coupling_modules_config.json`, `module_format: external_api`) | **Missing.** `convert_0d_to_1d.py` inserts it, and `FV1D_vessel`'s `api.process` names it |
| `volume_sum_2/nn` (coupling) | Missing. No CA model or test uses it. Drop, or add as a `volume_sum` version |
| `tissue_diffusion/nn` (diffusion_volume) | Missing. The library has `transport/tissue_diffusion_volume` and `_face` instead. Unused. Drop, or record the mapping |
| `test_units`, `vessel_type_test1`, `vessel_type_test2`, `init_states_test` (test_modules) | Already in CA `tests/test_inputs/test_modules` (#538). They stay in CA |
| `base_script.cellml`, `units.cellml` | Not modules: the generator's skeleton and shared units. They stay in libcuflynx |

Coupling modules whose names the libcuflynx code hard-codes:
- `FV1D_vessel` (`CVSCellMLGenerator` and `ModelParsers` test for it);
- `FV1D_volume_sum` (`ModelParsers` appends the `volume_sum_1D` rows itself);
- `FV1D_solver` (`convert_0d_to_1d.py`).

The library has `FV1D_vessel` and `FV1D_volume_sum`. However, its `FV1D_vessel` is `module_format: cpp` on `devel_module_tests` and `external_api` with #520's `api` block only on `fix/consistent-junctions` (library #4). See **D3**.

### 1.4 CA tests that need resource models

The tip of the stack is `feat/tests-on-lumped-vessels` (#540, which includes #538).

Shared fixtures in `tests/conftest.py`:
- `resources_dir` returns `<repo>/resources`;
- `_MODEL_INPUT_FILES` names 3compartment, SN_simple and test_init_states;
- `test_model_configs` names ports_test, 3compartment, simple_physiological, parasympathetic_model, test_fft, neonatal, the four generic_junction models, SN_simple, physiological, control_phys and aortic_bif_1d (cpp).

How many test files use each model:

| Model | Test files | Model | Test files |
|---|---|---|---|
| 3compartment | 28 | Lotka_Volterra | 12 |
| SN_simple | 7 | Simple_ODE_Benchmark | 7 |
| 3compartment_nonstiff / _extra_ops | 5 / 5 | simple_physiological | 5 |
| test_fft, aortic_bif_1d | 4 each | aortic_bif_hybrid_V1 | 3 |
| physiological, control_phys, generic_junction_test2_* (both), test_init_states, VanDerPol | 3 each | ports_test, pid_control, neonatal, parasympathetic_model, control_phys_asd, generic_junction_test_* | 2 each |
| 1 each | control_parasymp, delay_test, FTU_wCVS, FinalModel, aortic_bif_hybrid_V2, cerebral, cerebral_elic, cvs_model_with_arm_hybrid, elic, new_valve_p_est, 3compartment_q_lv_only, FitzHugh_Nagumo, Lotka_Volterra_multisub, NKE_pump, TwoMinima, Lotka_Volterra_forced_multi_trace, Goodwin, Teusink | | |

The 43 test files that read inputs from `resources_dir`:
conftest, interactive_tests/test_interactive_tutorial_notebooks, interactive_tests/test_notebook_paths_resolve, test_UQ, test_aadc_conditionals, test_aadc_coupled_pipeline, test_aadc_solvers, test_aadc_vs_casadi_3compartment, test_autogeneration, test_benchmarks, test_casadi_conditionals, test_cpp_template_generator, test_cross_experiment_references, test_emulator_solver_helper, test_emulator_training, test_external_module_library, test_external_user_funcs, test_generation_fixes, test_issue_fixes, test_lumped_vessels, test_module_instances, test_no_writes_into_package, test_obs_data_conventions, test_obs_data_vocabulary, test_operation_kwargs, test_param_id, test_prediction_features, test_protocol_shapes, test_sensitivity_analysis, test_sn_firing_threshold_recovery, test_solvers, test_unit_conversion, test_uq_on_emulator, test_user_inputs_archive, test_weight_zero_series_is_inert, test_vessel_configurations_1d, and others.

Run `git grep -l resources_dir upstream/feat/tests-on-lumped-vessels -- tests` for the exact list. Tests that only use `importlib.resources` or `generators/resources` (test_package_resources, test_package_syntax, test_config_schemas) are unaffected by removing `resources/`. They are affected by removing the built-ins (§2.3).

Also: `tests/data/external_coupling/resources/` (#520) is CA test data and stays.

### 1.5 CA code, docs and packaging that default to `resources/`

- **`utilities/paths.default_resources_dir()`** is `user_data_root()/resources`. It is the default `resources_dir` for:
  - `PrimitiveParsers` (`params_for_id_path`, `vessels_csv_abs_path`, `parameters_csv_abs_path`, the 0D/1D split file names);
  - `CVSCellMLGenerator` (where `*_parameters_unfinished.csv` is written);
  - `CVSCppGenerator`;
  - `paramID`;
  - `identifiabilityAnalysis`.

  `save_dated_user_inputs` **archives every run's config into `resources_dir`**. With no `resources/` in the checkout, that default must still work or be retargeted (§2.5).
- **`user_run_files/user_inputs.yaml`:** `file_prefix: 3compartment`, `input_param_file: 3compartment_parameters.csv`, and `param_id_obs_path: /home/farg967/.../circulatory_autogen/resources/3compartment_obs_data.json` (an absolute path to the developer's own checkout). `module_library_dirs` is only a commented example. `user_run_files/run_module_generator.sh` also uses the default.
- **Tutorial:** 8 docs (getting-started, design-model, example-3compartment-model, model-generation-simulation, parameter-identification, sensitivity-analysis, external-python-solvers, api/index) and the interactive notebooks (`generation_and_calibration*.ipynb`, `image_to_hemodynamics_model.ipynb`, `generate_param_array.py`, `image_proc_funcs.py`). These are tested by `interactive_tests/`.
- **benchmarks/:** README, `benchmark_specs.py`, `profile_benchmark.py` and `run_benchmarks.py` use FitzHugh_Nagumo, 3compartment, Goodwin, Teusink and `modifications.txt`.
- **`src/libcuflynx/scripts`:** 15 files. Most only document the default. `convert_0d_to_1d.py` has a hard-coded HPC path; `_cli.py` has help text.
- **`external_testing/full_pipeline_run.py`:** takes `resources_dir` as an argument.
- **README.md, CLAUDE.md:** describe `resources/` as the place for inputs ("Never commit `resources/user_inputs_<yymmdd>.yaml`").
- **Packaging:** `resources/` was never package data, so the wheel already lacks it. Removing `resources/` changes nothing for pip users. Removing the **built-ins** does (§2.3).
- **CUFLynx** (`/home/farg967/Documents/git_projects/CUFLynx`) doesn't rely on CA's `resources/` or its built-in modules:
  - its examples live in CUFLynx's own `resources/`;
  - every runner passes `resources_dir=os.path.dirname(params_path)`;
  - it never calls the generator;
  - the CI/release smoke tests check out CA only as the `ca_dir` source, identified by `pyproject.toml` + `user_run_files`, not `resources/`.

  Low risk. Keep `user_run_files/` (a repo marker in `paths._REPO_MARKERS`).

---

## 2. Migration plan

### 2.1 How CA finds the library

**What exists.**
- **#521** (merged): the config keys `module_library_dirs` and `use_builtin_modules`.
- **#538:** the environment variable `CUFLYNX_MODULE_LIBRARY` (`os.pathsep` list). It is the default for configs that set neither key, and it turns the built-ins off.
- **#538's `tests/conftest.py`:** sets the variable to `$CUFLYNX_MODULE_LIBRARY`, else to `../circulatory-autogen-modules/modules`, and appends `tests/test_inputs/test_modules`.
- **#538's CI** (`.github/workflows/tests.yml`): every job does `actions/checkout` of `physiomelinks/circulatory-autogen-modules` at `env.MODULE_LIBRARY_REF` (currently the **branch** `fix/consistent-junctions`) into `circulatory-autogen-modules/`.

The options:

| Option | For | Against |
|---|---|---|
| (a) Sibling checkout plus env var (#538 today) | Already written; a developer edits library and CA side by side | A branch ref drifts. Local runs use whatever is checked out. Nothing pins it |
| (b) Git submodule | Pinned by SHA in the tree | Submodule friction for every contributor. A 36 MB `.git` clone. PRs that bump it conflict |
| (c) Test fixture that clones the library at a pinned commit into a cache | Pinned; works with no setup | Network in the test session; another copy beside the developer's checkout |
| (d) Pinned git dependency: a data-only distribution of the library installed by pip | One mechanism for tests, CI **and** runtime (§2.3); pinned by SHA in `pyproject` | Needs a small packaging change in the library; needs care to avoid a dependency cycle |

**Recommendation: (a) with a pinned SHA now, and (d) when the built-ins go.**

1. **Pin now (in #538, or the next CA PR):**
   - `MODULE_LIBRARY_REF` becomes a full commit SHA, kept in one place: `tests/module_library_ref.txt`, or `[tool.libcuflynx] module_library_ref` in `pyproject.toml`.
   - CI reads it from there.
   - conftest prints the checked-out library's SHA. It warns (or fails, behind `--strict-library-ref`) when the SHA differs from the pin.
   - Developers keep the sibling checkout and the env override.
2. **When the built-ins are removed (§2.3):**
   - The library publishes a data-only distribution (working name `cuflynx-modules`) containing `modules/` and `system_models/` as package data. It has **no dependencies**, so there is no cycle: the library's `cam-testing` depends on libcuflynx, and libcuflynx's extra depends only on the data package.
   - It registers an entry point `libcuflynx.module_libraries`, `default = cuflynx_modules:modules_dir`.
   - libcuflynx's `ModuleSources` resolves modules in this order: the config, then `CUFLYNX_MODULE_LIBRARY`, then the entry point, then an error that says how to install.
   - CI then runs `pip install "cuflynx-modules @ git+https://github.com/physiomelinks/circulatory-autogen-modules@<sha>"` instead of the extra checkout steps. conftest uses the installed package unless the env var says otherwise.

### 2.2 Names: CA's old names vs the library's restructured names

There are three layers of rename:
1. **The junction removal:** `junction_migration.JUNCTION_TWINS` (CA #538), matched by library #4.
2. **Versions, the heart move and the 2026-10-06 restructure** (`cam_testing.library.legacy_renames()`: 176 pairs).
   - Sources: `tools/restructure_map.yaml`, `HEART_VERSION_RENAMES` and `tools/restructure_modules_renames.json` (115 entries).
   - The restructure drops the `nn` prefixes, giving names such as `cardiac_clock/Liang2009_v01` and `heart_effector/lv_ArgusUNPUBLISHED_v01`. Those names need **CA #548** (BC prefix only for vessels).
3. **The lumped twins:** `lumped_migration` (#538/#540). It changes only the `BC_type` and a few parameters; outputs keep their names.

#538 renamed CA's `resources/` against `fix/consistent-junctions`, which **doesn't have the 2026-10-06 restructure**. So as soon as CI moves to a library commit that has it, CA's arrays name module pairs that no longer exist.

**Options.**
- (i) **A scripted migration** of CA's test models and user_inputs to the new names, using a rename map the library publishes.
- (ii) **A legacy alias layer** in libcuflynx: on load, unknown `(type, version)` pairs are looked up in an alias table.

**Recommendation: (i), plus a user-facing migration CLI. No silent alias layer.**
- The library publishes **one machine-readable rename map**: `modules/renames.json`, which is `legacy_renames()` + `JUNCTION_TWINS`, each entry `{from, to, why, since}`. It is generated by a `cam_testing` command and checked by a structure test. The library is the single source of truth.
- libcuflynx gets `python -m libcuflynx.utilities.module_migration <array> [--parameters ...] [--check obs/params_for_id...]` (`cuflynx-migrate-model`). It generalises `junction_migration`: it reads `renames.json` from the active library, rewrites arrays in place, applies the junction collapse, and reports changed output names.
- Optionally, libcuflynx turns an unknown pair that appears in `renames.json` into an **error message** that names the new pair and the CLI. This gives users with old models a clear way forward, without hidden aliasing that would let old names live on in tests.
- Why not an alias layer: CA's tests would keep exercising names nobody else uses, and the alias table would have to stay in step with every library rename. A split-heart rebuild can't be expressed as an alias anyway.

### 2.3 Built-in modules: which stay in libcuflynx

| Group | Proposal |
|---|---|
| `base_script.cellml`, `units.cellml` | **Stay.** Generator infrastructure, not modules |
| `FV1D_vessel`, `FV1D_volume_sum`, `FV1D_solver` (coupling, `external_api` in #520) | **D3.** Recommended: **stay in libcuflynx as "core" modules**, always loaded even with `use_builtin_modules: false`, because the generator and parsers hard-code their names and the `api` contract (`process: FV1D_solver`) belongs to the C++ template generator. The library keeps copies for V&V/PhLynx, plus a structure test that they equal libcuflynx's core copies. (The alternative is for libcuflynx to keep nothing and require the library, so #520 changes to FV1D need a library PR in lockstep.) |
| `coupler/pv_0D_1D`, `imposter_1D/nn` | Move to the library (it already has them) |
| `test_modules` (4 entries), `Simple_ODE_Benchmark` | Stay in CA `tests/test_inputs/test_modules` (#538 did this) |
| Everything else (≈245 entries in 21 files) | **Remove** from libcuflynx; the library is the source. Once the tests run on the library (#538), the built-ins are untested duplicates |
| `volume_sum_2`, `tissue_diffusion` | Drop (unused); note them in the CHANGELOG |

`use_builtin_modules` then means "also load the core modules". It defaults to true and keeps working for `false` (the core modules are always there). Deprecate the key's old meaning in the CHANGELOG.

### 2.4 Gaps to close in the library first

1. **`reference/` can't depend on CA's built-ins.** `cam_testing/system.py::generate(reference=True)` generates `reference/` with CA's built-in modules, which happens whenever `use_builtin_modules` isn't false. Before CA removes the built-ins, choose one:
   - **(recommended)** generate each original's CellML once, with a pinned libcuflynx that still has the built-ins, and commit it as `reference/<model>.cellml`. Library #6's `equivalence.reference_cellml` then simulates it as is.
   - Or migrate `reference/` to library names (§2.2) and generate it on library modules. The monolithic `heart/vp*` versions still exist in the library. This no longer compares against CA's real modules, though.
2. **Re-import `reference/`** from CA master. #537 is merged, so SN_simple's xfail may be stale (also SN_simple_flat and SN_simple_supermodules).
3. **The 11 models that don't generate or are skipped:**
   - `generic_junction_test2_closed_loop` and `_open_loop`: add parameters, copied from `generic_junction_test_*_parameters.csv` as CA's tests already do, or point `input_param_file` at them. That leaves **9**.
   - `ports_test`, `test_init_states` and `Simple_ODE_Benchmark`: they test libcuflynx features with libcuflynx test modules. Recommended: **remove them from the library** and move them to CA `tests/test_inputs/models/` beside `test_modules` (**D2**). That leaves **6**.
   - `FinalModel`, `FTU_wCVS`, `new_valve_p_est`, `elic`, `cerebral` and `cerebral_elic` fail in CA too. CA's tests already list them as "no reference" (`test_lumped_vessels`). Either fix the model data in the library, or keep them as documented xfails (**D6**). No CA test needs them to generate.
4. **The C++/1D models** (`aortic_bif_1d`, `aortic_bif_0d`, `aortic_bif_hybrid_V1/V2`, `cvs_model_with_arm_hybrid`):
   - add a `coupled_1d` category with a run-only system test via library #5's C++ path;
   - needs `FV1D_solver` (or core in libcuflynx, D3) and CA #520 merged.
5. **Library #4** (`fix/consistent-junctions`) needs rebasing onto the restructured `devel_module_tests`. It predates `9f9985a`. CA #538's `MODULE_LIBRARY_REF` then moves to the rebased commit.
6. **Library #5** (external_api / FEniCS): needs CA #520 plus the #534–#536 stack. Independent of removing `resources/`, but it is where `external_api` support and the `FV1D_vessel` `api` blocks land on the restructured branch.
7. **`modules/renames.json`** (§2.2) and the data-only distribution with its entry point (§2.1, step 2). Tag a library release that CA can pin.
8. **The 29 data files** (§1.2) that are model data rather than test fixtures (the prediction_variables files, and Goodwin/Teusink if D4 says library).

### 2.5 Sequencing against the open PR stack

The current stack:
- `master` ← #535 ← #536 ← #541 ← #547;
- `master` ← #520 (on hold for Beatrice's 1D tutorial work) ← #538 ← #540;
- #548, #522, #550, #552 and #519 go straight to `master`;
- library: `main` ← #3 ← {#4, #5, #6}.

None of #522, #535, #536, #541, #547, #548, #550, #551 or #552 touches `resources/`, apart from #547's `test_reporting` reading `resources/Goodwin.cellml`. They need no changes beyond rebasing.

**The blocker.** #538 (the switch to the library) is stacked on #520, which is on hold. Everything below waits on #520 unless #538's library switch is split off. #538's library part needs `port_nodes` (because library #4 removes the junction types). It doesn't need #520's C++ generator, apart from the FV1D checks at a node. See **D7**.

#### CA PRs

**CA #548** (BC prefix only for vessels) to master. *Existing.*
- [ ] Merge before any CA test runs on the restructured library (`lv_*`, `Up*` non-vessel versions).

**CA #535 → #536** to master. *Existing.*
- [ ] Merge. Library #3 pins #536's branch.

**CA #538** (node multiports, tests on the library). *Existing; amend.*
- [ ] Pin `MODULE_LIBRARY_REF` to a **SHA**, not a branch, in a single file that CI and conftest read (§2.1).
- [ ] Re-migrate `resources/` to the post-restructure names once library #4 is rebased (§2.4 item 5). The migration is scripted with `legacy_renames` and checked by generating every model.
- [ ] Optional: generalise `junction_migration` into `module_migration`, reading the library's `renames.json` (it could also go in its own PR).
- [ ] If D7 says split: move the library switch (conftest, CI, env var, `test_modules`, `resources/` renames, `junction_migration`, `port_nodes`) into a PR onto master. Leave the C++-specific node checks on #538.

**CA #540** (tests on lumped vessels). *Existing; no change to its intent.*
- [ ] Rebase after #538.
- [ ] Its `lumped_migration` becomes the tool for building lumped test models from library models (see the next PR, under D1 option B).

**NEW "CA: tests read their models from the module library".** Base: after #540, or after the split-off library switch.
- [ ] conftest:
  - [ ] Replace the `resources_dir` fixture with `model_dir(prefix)`. It resolves `<library>/system_models/*/<prefix>/` (or `.../reference/`, D1), with `tests/test_inputs/models/<prefix>/` searched first for CA-only models.
  - [ ] `base_user_inputs` sets `resources_dir` per model.
  - [ ] `_MODEL_INPUT_FILES` and `test_model_configs` change accordingly.
- [ ] Move the test fixtures (§1.2) into `tests/test_inputs/`:
  - [ ] `3compartment_q_lv_only_params_for_id.csv`;
  - [ ] the two Lotka_Volterra obs_data variants;
  - [ ] `data_references.txt`;
  - [ ] ports_test, test_init_states and Simple_ODE_Benchmark (D2);
  - [ ] the C++/1D models, until the library has them (D5).
- [ ] D1 option B: build the lumped variants of #540's 14 models in a session fixture, using `lumped_migration --out-dir` from the library models. Cache them in `tests/test_outputs/`.
- [ ] Update the 43 test files that join `resources_dir` (mostly a fixture swap). Under D1 option A, also update the 15 files with hard-coded heart names and any calibrated numbers.
- [ ] Add a guard test: no file under `tests/` joins `<repo>/resources`.
- [ ] CI: run on the pinned library SHA. Make sure the full suite's results match #540's baseline table (2513 passed and so on).
- [ ] `benchmarks/` and #547's `test_reporting`: Goodwin/Teusink from their new home (D4).

**NEW "CA: remove `resources/`".** Base: after the previous PR.
- [ ] `git rm -r resources/`. Add `resources/` to `.gitignore` (users' own inputs and the dated `user_inputs_*.yaml` archives still default there).
- [ ] `paths.default_resources_dir()`:
  - [ ] keep `user_data_root()/resources` as the default **output** location for archives and `*_parameters_unfinished.csv`, and create it on demand;
  - [ ] give a clear error, naming `resources_dir` / `module_library_dirs` / the library, when a model's input files aren't found there.
- [ ] `user_run_files/user_inputs.yaml`:
  - [ ] `resources_dir` becomes the library's `system_models/closed_loop_cvs/3compartment` (relative to a sibling checkout, with a comment);
  - [ ] drop the absolute `param_id_obs_path`;
  - [ ] document `module_library_dirs` and `CUFLYNX_MODULE_LIBRARY`;
  - [ ] update `run_module_generator.sh`.
- [ ] Tutorial: the 8 docs and the notebooks point at the library models. Getting-started gets "clone circulatory-autogen-modules" (or `pip install libcuflynx[modules]` after the built-ins PR). Make sure `interactive_tests` pass.
- [ ] README, CLAUDE.md and the funcs_user READMEs: no `resources/` as the input home.
- [ ] `src/libcuflynx/scripts`: help text and defaults. Delete `convert_0d_to_1d.py`'s HPC path. Move or delete the two helper scripts.
- [ ] CHANGELOG: a migration note. It lists where each model went (the §1.1 table), `cuflynx-migrate-model`, and the dropped files.

**NEW "CA: built-in modules come from the module library".** Base: after "remove resources/". Can share a release with it.
- [ ] Remove the ≈245 built-in entries and their cellml. Keep the core modules (D3) and `base_script.cellml` / `units.cellml`.
- [ ] `ModuleSources`: core modules always load, and the default library resolves through the `libcuflynx.module_libraries` entry point (§2.1, step 2).
- [ ] `pyproject`:
  - [ ] `[project.optional-dependencies] modules = ["cuflynx-modules @ git+...@<sha>"]` (or a PyPI version);
  - [ ] the `test` extra includes it;
  - [ ] package-data shrinks to the core files.
- [ ] CI: install the data package instead of the checkout steps, while keeping the env override for library-developer jobs.
- [ ] Tests:
  - [ ] update `test_package_resources`, `test_external_module_library` (the duplicate / built-in cases) and `test_no_writes_into_package`;
  - [ ] add a test that a wheel install with the `modules` extra generates 3compartment.
- [ ] Errors for an unknown `(type, version)` that appears in `renames.json` name the new pair (§2.2).

**CA #520** (cpp template generator). *Existing; rebase when it comes off hold.*
- [ ] Its `resources/aortic_bif_0d_*` additions and its edits to the hybrid models go to the library's `coupled_1d` models (or to `tests/test_inputs/models/`, D5), not to `resources/`.
- [ ] Its `coupling_modules_config.json` `FV1D_*` entries stay in libcuflynx core (D3).
- [ ] `tests/data/external_coupling/resources/` stays.

**CA #541 / #547** (calibration workflow, methods LaTeX). *Existing.*
- [ ] #541 already takes `--module-library-dir`. No change.
- [ ] #547: `test_reporting` reads `resources/Goodwin.cellml`. Repoint it per D4 when rebasing onto the "tests read their models" PR.

#### Library PRs

**Library #3** (`devel_module_tests` → `main`). *Existing.* After CA #536 merges.
- [ ] Repoint `requirements.txt` from `feat/prediction-item-features` to master.

**Library #4** (consistent junctions). *Existing; rebase.*
- [ ] Rebase onto the restructured #3 (it predates `9f9985a`).
- [ ] Merge together with CA #538. CA's `MODULE_LIBRARY_REF` moves to its merge SHA.

**NEW "library: CA test models and references".** Base: after #4 (or #3).
- [ ] Generate each `reference/` as committed CellML with a pinned built-ins libcuflynx, using #6's `reference_cellml` (§2.4 item 1). Needs #6 first, or ports that part of it.
- [ ] Re-import `reference/` from CA master, and re-check SN_simple's xfails.
- [ ] generic_junction_test2_*: add parameters.
- [ ] ports_test, test_init_states and Simple_ODE_Benchmark: remove them (D2) or keep them as skips.
- [ ] Add the prediction_variables files (and Goodwin/Teusink if D4 says library).
- [ ] D1 option B: add CA's test form beside each model (e.g. `<model>/libcuflynx/`, or `reference/` migrated to library names), with an identity-map equivalence test against `reference/`.
- [ ] Add `modules/renames.json`, generated and checked (§2.2).

**Library #5** (external_api / FEniCS) and **NEW "library: FV1D_solver and coupled_1d system models".** After CA #520.
- [ ] FV1D_solver module (if D3 says the library), or the core-equality test (if D3 says core).
- [ ] The five C++/1D models in `system_models/coupled_1d/`, run-only.

**NEW "library: data-only distribution + release tag".** Before CA's built-ins PR.
- [ ] A `cuflynx-modules` package: `modules/` and `system_models/` as package data, no dependencies, and the `libcuflynx.module_libraries` entry point.
- [ ] Tag a release. CA pins its SHA or tag.

**Library #6** (cam_testing from another repo). *Existing.* Independent. Its `reference_cellml` is used above.

**Order:**
1. #548, then #535 → #536 (CA).
2. Library #3, then library #4 (rebased) together with CA #538 (pinned SHA, re-migrated names), then CA #540.
3. Library "CA test models and references", then CA "tests read their models from the library", then CA "remove resources/".
4. Library data distribution + tag, then CA "built-in modules come from the library".
5. C++/1D: CA #520, then library #5 + coupled_1d, then the C++ models move from CA's `tests/test_inputs` into the library.

### 2.6 Risks

- **Users who rely on `resources/` paths:**
  - `file_prefix` with the default `resources_dir`, or absolute paths like the one in `user_inputs.yaml`;
  - their own models in `resources/` that git now sees as deleted files on pull. Mitigations:
  - the CHANGELOG migration section;
  - the clear error when inputs are missing;
  - `resources/` added to `.gitignore` so users' own files survive.
  - Consider one release that keeps `resources/` with a README stub ("moved to circulatory-autogen-modules") before deleting it.
- **Old module names** in users' models: `cuflynx-migrate-model`, and errors that name the new pair. The split-heart rebuilt models change output names, so users moving to the library's models (rather than renaming their own) must update obs_data and params_for_id.
- **Tutorials and notebooks:** paths break silently. `interactive_tests` must run in the "remove resources" PR's CI.
- **The wheel without built-ins:** `pip install libcuflynx` alone can't generate physiological models. This needs the `modules` extra and the entry point, and a clear error when neither the env var nor the extra is present. A git-URL extra can't go to PyPI as is, so either publish `cuflynx-modules` to PyPI or document the git install (**D8**).
- **Version skew:** libcuflynx and the library change together (port semantics, `api` blocks, the instance layout). Pin by SHA in both directions:
  - CA pins the library SHA;
  - the library pins the libcuflynx ref in `requirements.txt`;
  - each side's CI uses the other's pin;
  - bump both pins in paired PRs.
- **CUFLynx:** low risk (§1.5). Its release smoke test checks out CA master and its own resources. Re-run CUFLynx `integration.yml` after the built-ins PR. PhLynx (#595) loads the library itself.
- **Library equivalence tests** break the day the built-ins go, unless §2.4 item 1 is done first.
- **Generation cost:** the lumped models are about 7× slower to generate on libcellml 0.6.3 (#540). Building them in a session fixture adds to the suite's time unless they are cached.
- **Other sessions' open work:** "git-projects-2b" has #536/#538/#549, and "CVS-ANS" has #520/library #5. The amendments to #538 and #520 above should go through those sessions.

---

## 3. Decisions needed

- **D1: Which form of each model CA's tests use.**
  - **A.** The library's rebuilt models (split heart). This is "use the library's models" literally, but needs changes to 15 test files and their numbers.
  - **B (recommended for the first step).** CA's own models, renamed to library names (the current `resources/` after #538/#540), hosted in the library beside each model and equivalence-checked. CA's tests change only how they find files.

  Moving the tests to A could come later.
- **D2:** Move ports_test, test_init_states and Simple_ODE_Benchmark (libcuflynx feature tests on CA test modules) into CA `tests/test_inputs/models/`, rather than keeping them as library skips. Recommended: yes.
- **D3:** `FV1D_vessel`, `FV1D_volume_sum` and `FV1D_solver` stay as libcuflynx core modules (recommended), or live only in the library.
- **D4:** Goodwin/Teusink (external CellML benchmarks): library `system_models/benchmarks` via `reference_cellml`, or CA `benchmarks/data/`.
- **D5:** The C++/1D models: wait for #520 and the library's C++ path, keeping them in CA `tests/test_inputs/models/` meanwhile (recommended), or block removing `resources/` on #520.
- **D6:** The six models that don't generate (FinalModel, FTU_wCVS, new_valve_p_est, elic, cerebral, cerebral_elic): fix them in the library, or keep them as documented xfails. Also: keep or drop the unused `pulmonary_lobe_imped` obs_data.
- **D7:** Split #538's switch to the library off #520 onto master, so it doesn't wait for Beatrice's 1D work.
- **D8:** How the library reaches pip users: a data-only distribution on PyPI, a git-URL extra, or a fetch-on-demand cache. Plus its package name.
- **D9:** Whether to ship one release with a `resources/` stub before deleting it.
- **D10:** Whether to remove the built-in modules in the same release as `resources/`, or one release later.

## 4. Decisions taken (user, 2026-10-07)

- **D1: a mixture, for coverage of both.** Some CA tests use the library's rebuilt split-heart system models (option A: their output names and expected numbers change), and others use CA's original unsplit-heart models, hosted in the library under library module names (option B). Split the 43 test files so that every libcuflynx feature area (param_id, solvers, sensitivity analysis, emulators, junctions, supermodules) is covered by at least one test on each form. A table in the "tests read their models from the library" PR records which test uses which form.
- **D2:** yes. ports_test, test_init_states and Simple_ODE_Benchmark move to CA `tests/test_inputs/models/`.
- **D3:** FV1D_vessel, FV1D_volume_sum and FV1D_solver stay as libcuflynx core modules.
- **D4:** Goodwin and Teusink go to the library (`system_models/benchmarks`, via `reference_cellml`).
- **D5:** the C++/1D models stay in CA `tests/test_inputs/models/` until #520 and the library's C++ path land.
- **D6:** the six models that don't generate stay as documented xfails; the unused `pulmonary_lobe_imped` obs_data is dropped.
- **D7:** split #538's switch to the library off #520 onto master.
- **D8: open.** It is under discussion; the options and measurements are below.
- **D9:** one release ships a `resources/` README stub before the directory is deleted.
- **D10:** the built-in modules are removed one release after `resources/`.

### D8 options (measured 2026-10-07)

**Size:** the distributable library (`modules/` + `system_models/`, excluding the generated results, plots, HTML and .omex) is 3,702 files: 39 MB uncompressed, 9.9 MB as .tar.gz. The core model files alone (CellML, configs, parameters, module arrays) are 11 MB uncompressed, roughly 2–3 MB compressed.
