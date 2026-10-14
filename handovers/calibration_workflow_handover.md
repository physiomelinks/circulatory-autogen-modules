# Handover: calibration workflows for supermodules (libcuflynx + CUFLynx)

Written 2026-10-05 from the circulatory-autogen-modules walkthrough session. Hand it to an agent that will implement the feature in **libcuflynx** (repo `circulatory_autogen`) and **CUFLynx**, and wire it into **circulatory-autogen-modules**.

## 1. What the user asked for (verbatim)

> "In the supermodule instances we need to have a calibration_workflow.json file that details the order of calibrations that are done for submodules in the supermodule. Then libcuflynx should be able to repeat the workflow of calibrations with one command. This calibration_workflow.json should also be stored in a CUFLynx run, where you can switch tabs from submodules to supermodules. Implement this in CUFLynx as well. Take note that in the future I want to be able to store parameter distributions at each submodule calibration to use as priors for the supermodule calibration."

There is also a standing rule from the same session:

> "I want the calibrations to data to be done properly through an obs_data file and to use libcuflynx so a user could change the obs_data and redo. Make this the case for all recalibrations."

## 2. Why it is needed

The module library (`circulatory-autogen-modules`) has submodules, such as ion channels. Each is calibrated alone to its own data, in an instance with `<instance>_obs_data.json` and `<instance>_params_for_id.csv`. Supermodules, such as the sympathetic-neuron soma, assemble them.

Some quantities can only be calibrated on the assembled supermodule, after the submodule calibrations. The case that triggered this is the soma's resting Na/K balance:
- On review, the Na/K pump was given realistic values: density 750 /µm², turnover 9 /s, Nai_mid 15 mM (commit 9492c41).
- That cut the pump current about 125×, so the soma no longer holds Nai ≈ 13 mM / Ki ≈ 200 mM / V ≈ −70 mV at rest.
- The leak conductances, g_leak_Na and g_leak_K, must be recalibrated **on the soma**, with every channel's own calibrated values in place.

That is a workflow: run the channel calibrations, apply them, then calibrate the soma.

## 3. Repos, branches and rules

**circulatory_autogen (libcuflynx)**
- Location: `/home/farg967/Documents/git_projects/circulatory_autogen`. Remote `upstream` = physiomelinks, default branch `master`.
- **Rule:** never work in the main checkout. Create `git worktree add .claude/worktrees/<name> -b <branch> upstream/master`, with a venv inside the worktree. Write tests, push to upstream, and open a PR with `gh pr create --base master`.
- Relevant open PRs:
  - **#534** (JSON vessel arrays and supermodules). The supermodule machinery is here; the new feature probably needs to stack on it, or wait for it.
  - **#535** (module versions/instances; obs_data_name).
  - **#536** (prediction items as SA/emulator features; stacked on #535).
  - **#538** (multiports).
  - Decide the base branch after reading them; stacking on #534/#535 is likely.
- Existing **`src/libcuflynx/scripts/sequential_param_id_run_script.py`** / `cuflynx-sequential-param-id`: a *different*, unimplemented idea (staged identifiability fitting within one model; `SequentialParamID` is missing). Don't confuse the two. Either keep them separate, or say in the PR how they relate.
- Calibration entry points: `src/libcuflynx/scripts/param_id_run_script.py` and `src/libcuflynx/param_id/paramID.py`. Also `train_emulator_run_script.py`, and the yaml-driven pipeline used by the sympathetic_neuron repo.

**CUFLynx**
- Location: `/home/farg967/Documents/git_projects/CUFLynx`, remote physiomelinks/CUFLynx. Vue 3 + Vite + PrimeVue (`apps/web`), FastAPI (`apps/api`), pywebview desktop shell.
- **Read `CLAUDE.md` first.** Key rules:
  - CUFLynx uses libcuflynx for *all* backend work.
  - Never spell a CA module in an import; use `ca_from(...)` / `ca_import(...)` from `apps/api/ca_imports.py`.
  - No OpenCOR, ever.
  - The CA dir is selectable at runtime.
- Open PRs from this session: #364 (OMEX import picks obs_data/params_for_id by name), #365 (validation in the Analysis tab), #366 (prediction items as SA/emulator features).
- Same rule: work in a worktree, add tests, open a PR.

**circulatory-autogen-modules**
- Location: `/home/farg967/Documents/git_projects/circulatory-autogen-modules`, branch `devel_module_tests` (draft PR #3).
- Don't stage the root `README.md`, `modules/system/` or `modules/poiseuille_transport/`; they belong to another session.
- Use the standard libcellml; no `PYTHONPATH` fast build.
- Tests: `venv/bin/python -m pytest tests/test_structure.py`, and `tests/test_modules.py --module <name> --include-unreviewed`.
- The repo venv has an editable install of libcuflynx from the worktree `circulatory_autogen/.claude/worktrees/prediction-features` (PR #536). `requirements.txt` pins #536.

**Commits and PRs everywhere**
- Commit messages end with:
  ```
  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01VVEqV1kahn5N7ekiZc7hHd
  ```
- PR bodies end with `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.
- **Never `pkill -f` / `pgrep -f` with a pattern that matches your own command line.**

## 4. What already exists in the module library (build on it)

**Instance layout:** `modules/<category>/<module_type>/versions/<version>/instances/<instance>/` holds:
- `<instance>_parameters.csv`;
- optionally `<instance>_obs_data.json`, whose `obs_data_name` equals the directory name;
- `<instance>_params_for_id.csv`;
- `<instance>_calibrated_parameters.csv` and `<instance>_calibration.json`, both written by calibration;
- `source_figures/`.

The rules are in `modules/README.md` ("Instances", "Source figures") and enforced by `modules/directory_schema.json` with `tests/test_structure.py`.

**Calibrating one instance:** `cam_testing/calibrate.py` drives libcuflynx. Its multi-segment protocols were fixed in commit e9f5910. The calibration test is `validation_test_calibrate` in `tests/test_modules.py`. Redo it with:
`venv/bin/python -m pytest tests/test_modules.py --component <mt>/<version> -k 'calibrate and <instance>' --include-unreviewed`

**Applying a calibration:** `venv/bin/python -m cam_testing.calibrate apply <mt>/<version> <instance> [--dry-run]` (commit eb09858, test `tests/test_calibrate_apply.py`). It copies the calibrated values into:
- the version's default instance;
- every supermodule using it, recursively (`<var>_<submodule path>`);
- the monolithic counterparts, found through `supermodule.equivalent` output maps;
- system models (`<var>_<vessel>[_<path>]`).

Its references read `Calibrated (libcuflynx) to …obs_data.json`, using the optional obs_data keys `calibration_note` and `reference_key`. **This name mapping is exactly what a workflow needs between steps.** Reuse it, or move it into libcuflynx.

**Supermodule naming:**
- A supermodule generated alone is vessel `mod`. libcuflynx expands it into `mod_<submodule>`, and a submodule parameter `<var>` becomes `<var>_mod_<submodule>`.
- In the supermodule's own instance CSV it is `<var>_<submodule>`.
- See `cam_testing/harness.py` (`supermodule_model_name`, `model_parameter_name`).

**Example of a real submodule calibration:** `modules/cell/ion_channels/i_M/versions/Argus2026_v01/instances/davis2020_wistar/`, which calibrates rho_M to Davis 2020's XE-991-sensitive relaxation. A second instance, `davis2020_wistar_kinetics`, refits tau_w and is **not applied**.

**The target supermodule:** `modules/cell/neuron/soma/versions/sympathetic/`. Its submodules are listed in `soma_sympathetic_modules_config.json`: membrane, Na_K, Ca, 13 channels, IP3R, RyR. There is also a monolithic equivalent, `soma/sympathetic_monolithic_v01`, and the neuron supermodule `modules/cell/neuron/versions/sympathetic/` contains the soma.

## 5. Proposed design (refine as needed, and record the decisions in the PRs)

### 5.1 `calibration_workflow.json`

It lives in a **supermodule instance** directory, e.g. `modules/cell/neuron/soma/versions/sympathetic/instances/rest_balance/calibration_workflow.json`. That instance also holds the supermodule-level obs_data and params_for_id for its own step.

```json
{
  "workflow_name": "soma_sympathetic_rest_balance",
  "description": "Channel calibrations, then the soma's resting Na/K balance",
  "steps": [
    {
      "id": "i_M_g",
      "target": {"module_type": "i_M", "version": "Argus2026_v01", "instance": "davis2020_wistar"},
      "submodule": "i_M",
      "depends_on": [],
      "store_distribution": false
    },
    {
      "id": "rest_balance",
      "target": {"module_type": "soma", "version": "sympathetic", "instance": "rest_balance"},
      "depends_on": ["i_M_g"],
      "fixed_from": ["i_M_g"],
      "priors_from": [],
      "store_distribution": false
    }
  ]
}
```

- **A step** calibrates one instance, using that instance's own obs_data, params_for_id and calibration settings.
  - `submodule` names where a submodule step's results land in the supermodule; it can be inferred.
  - A supermodule-level step targets the supermodule itself.
- **Order:** `depends_on` defines it; run steps topologically.
- **`fixed_from`:** before a step runs, the calibrated values of the listed earlier steps are applied as fixed parameter values in that step's model. The submodule-to-supermodule naming is as in `calibrate apply`.
- **Future, design for it now:** `store_distribution` saves each step's parameter distribution (posterior samples or a fitted summary), and `priors_from` uses those distributions as priors in a later step instead of fixing the values.
  - Reserve the keys and validate them.
  - If they are set before priors are supported, fail with a clear "not yet supported" message.
  - Check what libcuflynx's `param_id` supports for priors (MCMC prior specification, `pymc_backend.py`); if it is little work, implement a Gaussian or uniform prior from the stored samples, but say so.
- **Outputs**, per run: `workflow_output/<step_id>/` holds the step's normal param_id outputs, plus a top-level `workflow_result.json` with the final merged parameter set, the steps' provenance (obs_data hashes, times) and the order run.
  - Optionally write the merged set back as `<instance>_calibrated_parameters.csv`, for the module library's `calibrate apply`.

### 5.2 libcuflynx

- **Parser and validator** for `calibration_workflow.json`: a JSON schema, discoverable like CA's other schemas, which CUFLynx consumes.
- **One command**, e.g. `cuflynx-calibration-workflow <calibration_workflow.json> [--from-step X] [--only X] [--dry-run]`, plus a Python API (`run_calibration_workflow(path, ...)`) that CUFLynx can call.
  - It must work under MPI the way `param_id_run_script` does.
  - It resolves each step's model from the module library (`module_library_dirs`, as `#535`/`#534` use).
- **Each step reuses the normal param_id path**, the same as a single instance calibration. Don't fork the optimiser logic.
- **Tests:** a small two-step workflow on a test module library fixture. For example, calibrate a submodule parameter to data, then a supermodule parameter with the first fixed, checking order, value propagation and outputs.
- **Docs** in `tutorial/docs/parameter-identification.md`.
- **Rough edges found in this session** (from the agent that built the i_M calibrations; worth fixing while here, or as separate issues):
  - `sp_minimize` fails with "Gradient not available" for cellml models without forward sensitivity analysis. It should fall back to finite differences, as the emulator path does.
  - The forward-sensitivity gradient refuses cross-segment references.
  - A data item with weight 0 still needs a value and std; allow `value: null`.
  - There is no series "normalise" / "minus initial" operation.
  - There is no public API to evaluate features at given parameters across segments; the agent used `get_cost_obs_and_pred_from_params`, `accumulating_temp_results`, `evaluating_segment` and `_align_series_to_ground_truth`.

### 5.3 CUFLynx

- **Load** a workflow, from a file or from an OMEX/run that contains one, and **store** `calibration_workflow.json` in the CUFLynx run.
- **UI: tabs to switch from submodules to supermodules.** One tab per step, plus the supermodule. Each shows that step's model, obs_data, params_for_id, sliders and calibration results; per the user, they switch between submodule and supermodule views.
  - The supermodule tab shows which values came fixed from earlier steps, and later which priors.
- **Run the whole workflow or one step** through libcuflynx's API (backend in `apps/api`, using `ca_from`/`ca_import`). Show progress per step.
- **Tests:** API tests for workflow load/store/run with a small fixture, and web component tests if the repo has them.
- Keep the existing OMEX conventions (#364): CUFLynx picks members by name. If the module library's per-instance OMEX (`tools/build_instance_omex.py` in circulatory-autogen-modules) should carry the workflow and its step instances, coordinate that.

### 5.4 circulatory-autogen-modules

- Allow `calibration_workflow.json` in supermodule instances: `modules/directory_schema.json`, `tests/test_structure.py` (instance file rules), and a `modules/README.md` section documenting the format and the one-command redo.
- **A runner** in `cam_testing`, e.g. `python -m cam_testing.calibrate workflow <mt>/<version> <instance>`, that calls libcuflynx's workflow command and then `calibrate apply` for the merged result.
- **The report** (`cam_testing/report.py`, `templates/version.html`) shows a supermodule instance's workflow: the step list with links to each step's instance and its calibration result.
- **The first real workflow:** `soma/sympathetic` instance `rest_balance`.
  - Steps: the existing channel calibrations (currently only i_M `davis2020_wistar`), then the soma step.
  - The soma step calibrates **g_leak_Na and g_leak_K** (rho/gamma names in the leak modules; check `i_leak_Na` / `i_leak_K`) to resting targets:
    - Nai ≈ 13 mM (Boron & Boulpaep 2017);
    - Ki ≈ 200 mM (Brown & Scholfield 1974 / Galvan 1984; see the SN_Na_K_concentrations_soma references);
    - V_rest ≈ −70 mV (Davis 2020: Wistar −72.6 ± 1.0 mV).
    - Long run at 0 pA; stds from the sources.
  - Put the targets in the instance's obs_data with source_figures (the source-figures rule applies), and check whether the rest-check run in `cam_testing/checks.py` (`rest_check`) helps.
  - **Ask the user** before adopting the calibrated leak values as defaults: compare first, as was done for i_M. The leak modules (i_leak_Na, i_leak_K) are the next modules in the user's walkthrough, so coordinate with whoever runs that review.

## 6. Acceptance

- One command reproduces the soma `rest_balance` workflow from the committed files, end to end, in the module library. Editing an obs_data file and rerunning changes the result.
- libcuflynx PR with tests and docs; CUFLynx PR with tests; module-library commit(s) with structure tests passing.
- The priors extension is designed in (keys reserved, validated, documented), even if not implemented.
- Report back the PR links, the command, how the soma `rest_balance` result compares with the current defaults (rest V, Nai, Ki and the 30 s rest-check spikes), and any decisions the user must make.

## 7. Current state when this was written

- Soma defaults since 9492c41 (realistic pump): the 30 s rest check fires at 0.2 Hz. Over 600 s at 0 pA, Nai climbs from 13 to about 44 mM; it is not in balance. Before the pump change it was the reverse (Nai 13 → 7.3 mM; the pump was too strong).
- A separate SN_full recalibration (sympathetic_neuron repo, branch `review/sn-model-changes`, PR FinbarArgus/sympathetic_neuron#10) is running and expected to finish around 22:30–23:00 on 2026-10-05. It doesn't include the pump/NCX changes. Don't touch that worktree's model files while it runs.
