# Handover: getting circulatory-autogen-modules ready to share

Written on 2026-10-08 for a fresh agent. **The goal:** this repository is the tested module library for libcuflynx (circulatory_autogen), PhLynx and CUFLynx. Get it to the point where it can be shared with users: they can find a module, trust it (tests, citations, provenance), add or change one the right way, and test it. This handover says where things stand, what rules already exist, and what is left.

Read first, in this order:
1. `modules/README.md`. It is the authoritative contributor guide: layout, placement and naming, versions, instances, supermodules, external models, what runs, the spec files, the archives, the schema, and "Adding things".
2. `README.md`: reports, tests, parameter ranges, failure risk, system models, the PhLynx → CUFLynx pipeline, using cam_testing from another repo, running locally, CI.
3. The memory index, `~/.claude/projects/-home-farg967-Documents-git-projects-circulatory-autogen-modules/memory/MEMORY.md`. It loads automatically in this project and holds the user's standing rules, summarised in §3.

## 1. Where things stand

**Branch `devel_module_tests`** is at 12ffe5f or later; an agent may have pushed two more commits, see §6. **PR #3** (`devel_module_tests` → `main`) is the **single consolidated PR, marked ready for review**. The user intends to merge it, then build on `main`. It contains:
- **Library PRs:** #5 (external FEniCS coupling), #6 (cam_testing usable from another repository; the test suite now lives in `cam_testing/suite`, and the pytest options come from the `cam_testing` plugin entry point), #7 (vessel arrays → module arrays) and #8 (GE_capillary fix).
- **The 2026-10-06 restructure:**
  - mechanism-named module types; source-named versions; no `nn` versions;
  - bib keys rekeyed to lowercase `<surname><year><firstword>`;
  - `required_citations` in every config;
  - `reference/` copies checked against circulatory_autogen's originals.
- **The sympathetic-neuron (SN) review** of 22 modules; the PR #3 description summarises it.

**Not in #3:** library **PR #4** (junction removal, lumped vessels; still open). It needs libcuflynx `feat/node-multiports` (CA #538) and predates the restructure, so it is redone later as its own PR on top of `main`. Leave #4's branch untouched.

**Numbers:**
- 372 versions; **39 reviewed** (`reviewed: true` in the spec), almost all cell (SN), coupling and benchmarks.
- Versions per category: vessels 197, cell 69, control 35, boundary_conditions 19, heart 18, coupling 11, respiratory 10, organs 5, benchmarks 5, transport 3.
- `test_structure`: about 911 passed, 9 skipped, 35 xfailed.
- Pre-existing failures:
  - heart: `vp_simple` ×6, the `vp*` baselines ×5, the `Argus2026_v01` timestep test, `3compartment_supermodules`;
  - vessels: the `capillary/pp_micro` baseline, the BC sweeps of `flow_merge/vp` and `flow_split/pv`;
  - cell: the BC sweep of `soma/sympathetic_monolithic_v01`;
  - `test_report_features::test_supermodule_structure_of_the_sympathetic_neuron`, which expects the axon's old `C_axon`/`R_axon` and needs updating to the geometry parameters.
- CI (`.github/workflows/module-tests.yml`): structure, plus one job per top-level `modules/` directory. Failures of unreviewed versions are warnings (`tools/ci_reviewed_failures.py`); on `main` the Pages site is deployed. The category jobs currently fail on the blocking errors above.

## 2. What "ready to share" needs (the work)

Do these roughly in this order. Confirm scope with the user before big changes, and ask, don't assume (see §3).

**A. A user-facing front door.**
- The README files are complete but contributor-oriented. Add a short **user guide**: what the library is; how to browse the reports (the Pages site); how to use a module or system model from libcuflynx/CUFLynx (`module_library_dirs`, `CUFLYNX_MODULE_LIBRARY`); how to cite (`required_citations`, the report's citation block); the licence (CC0-1.0 per module, `licence` key).
- Keep `modules/README.md` as the contributor guide. Make sure "Adding things" walks through a full new module: the files, config keys (licence, creator, required_citations), parameters with references and the sourced column, spec files, instances, source figures, tests to run, and the report to check. Add a worked example, e.g. a small ion channel.
- Check that every rule in §3 is in the docs, not only in memory.

**B. Green CI and honest status.**
- Fix or document every blocking failure listed in §1 so `main` CI passes. Each must be fixed, or recorded as a known issue / expected failure with a reason. Never weaken a tolerance to make something pass.
- Unreviewed versions show as unreviewed. The report's uniform 11 columns must not imply passes that didn't happen.

**C. Review coverage.** 333 versions are unreviewed. Reviewing all of them before sharing is unrealistic. **Agree with the user** which set must be reviewed for a first public release (e.g. everything the shipped system models use), and how unreviewed ones are labelled. The walkthrough process is §3. Earlier reviews of system-model groups have notes in `reviews/*.yaml`.

**D. Citations and provenance.**
- 274 versions carry best-guess `required_citations_uncertain` (mostly vessel, junction and microvasculature families). Two bib entries need their details checked: `park2020investigating` and `westerhof2009arterial`.
- Instances listed in `MISSING_SOURCE_FIGURES` (in `cam_testing/suite/test_structure.py`) need source screenshots. The list may only shrink. The Wanaverbecq 2003 PMCA calibration instance was added to it because the PDF isn't on disk.
- Owner-original models cite the placeholders `argus2026sympathetic`, `argus2026unpublished`, `gee2026unpublished` and `argus2026circulatory`. They are replaced when the papers exist.

**E. Distribution and versioning.**
- **Decision D8, how pip users get the library, is still open.** The options are measured in `handovers/ca_resources_to_library_plan.md`:
  - a PyPI data package (recommended);
  - a git-URL dependency (PyPI forbids it in libcuflynx's own dependencies);
  - download on demand (9.9 MB archive; ~1.6 s at 50 Mbit/s, ~8 s at 10 Mbit/s; awkward on offline HPC).
- The library needs **release tags**: circulatory_autogen's tests will pin a library commit or tag.

**F. circulatory_autogen moves onto this library.** The user decided that circulatory_autogen's `resources/` will be removed, its models taken from here, and its tests run against this repo. The plan is `handovers/ca_resources_to_library_plan.md`. D1–D7, D9 and D10 are decided, see §4 there; D8 is open. The library-side work is to:
- freeze each `reference/` as committed CellML;
- host CA's original test models under library names, with a mixture of original and rebuilt models in CA's tests;
- add Goodwin/Teusink to `system_models/benchmarks`;
- tag a release.

The CA-side PRs belong to other sessions (§5).

**G. Unfinished SN items** (the user's sympathetic-neuron work):
- **The soma fires at rest** (about 2.8–3 Hz) after the PMCA was added. The leaks need rebalancing through a **calibration workflow** (`calibration_workflow.json`: ordered submodule → supermodule calibrations with one command, CUFLynx tabs, later per-step priors). It is specified in `handovers/calibration_workflow_handover.md` and owned by another agent or session; check with the user.
- The SN decisions review page, an HTML of every accepted SN decision, which the user asked for once the SN model is final.
- PR #4 redo (above).

## 3. Standing rules (from the user; mostly in memory and the READMEs)

**Placement and naming** (`modules/README.md` "Placement and naming"):
- A module type is named by mechanism and placed where the mechanism exists (`cell/Ca_handling`); only anatomical parts nest.
- Different maths makes a **version** named `<system>_<compartment>_<Source><Year>_vXX` (or `<Source><Year>_vXX`; `ArgusUNPUBLISHED_vXX` for owner-modified; `_SI`/`_OLD`). Same maths with different values makes an **instance**.
- No surnames in type names, no repeated directory names, no plain `nn` versions (except libcuflynx's FV1D allowlist).
- Versions of a type may differ in ports, documented.

**Per version:**
- **Config:** `licence` CC0-1.0; `creator` (ask the user per module, never guess); `required_citations` (+ `_uncertain`).
- **Parameters:** each CSV row has a `data_reference` (`<bibkey>; …`, `definitional; …` or `unreferenced_placeholder; …`) and a `sourced` yes/no. Write large departures from the literature in **bold** (`**…**`). Prefer primary sources over secondary model papers; Tao 2011 values are not trusted for now in SN modules.
- **Reports:** only `python -m cam_testing.report` writes them, from repo files.
- **Calibration:**
  - Any value fitted or derived from data comes from a committed obs_data instance through libcuflynx, never a hand calculation.
  - `python -m cam_testing.calibrate apply <mt>/<v> <instance>` copies it into the defaults, supermodules and systems.
  - Data extracted from a paper or book needs a cropped screenshot in `instances/<i>/source_figures/` (with `source_figures.json`).

**The walkthrough with the user** (how modules get reviewed):
- **One module at a time:** regenerate its report, then ask the questions, with the report's `file://` link as the **last line before the questions**.
- **Ask about:** references (trace values to the source, including the owner's git history), validation and calibration data, tests, the creator, and the known issues.
- **A decision for one module is not applied to others without asking.**
- **Model changes in SN modules are mirrored into the sister repo** (`sympathetic_neuron`, SN_full, PR #10, one step per change with an I_in sweep and a report). SN_full uses its own soma, axon and varicosity, never circulatory_autogen built-ins.

**Tests and tools:**
- Use the standard libcellml. The "fast" PYTHONPATH build breaks generators and units.
- The time-step stability refinement rule is in memory (`stability-refine-rule.md`).
- **Never `pkill -f`/`pgrep -f`** with patterns that match your own command line.

**Git:**
- Never touch the main circulatory_autogen checkout; CA changes go in a new worktree off `upstream/master`, then a PR.
- Merge commits only on shared PR branches (no rebase or force-push).
- Commit trailers: `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` and `Claude-Session: <session URL>`. PR bodies end with `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.
- **Ask before filing issues on PhLynx.**
- **Don't stage other sessions' files.** The root `README.md` has another session's uncommitted line (`system/poiseuille/`); `modules/system/`, `modules/poiseuille_transport/` and editor swap files also belong to others. To commit a README change without that line, stage a blob with `git hash-object -w` + `git update-index --cacheinfo`.

## 4. Environment

- **Venv:** the repo venv (`venv/`) has an editable libcuflynx from `circulatory_autogen/.claude/worktrees/prediction-features` (branch `feat/prediction-item-features`, PR #536), which belongs to another session. **Don't change that worktree's checkout**; ask that session first.
- **The libcuflynx pin:** `requirements.txt` pins `integration/modules-external-coupling` (#536 + #520 + #548), as PR #5 needs.
- **After pulling:** run `venv/bin/pip install -e . --no-deps` if the pytest options (`--include-unreviewed`, …) are missing.
- **Disk is very tight:** `/` often has under 0.5 GB free; `/home` has tens of GB. Put `TMPDIR` and pytest `--basetemp` under `$HOME` for anything big. CUFLynx desktop instances leave 2 GB `/tmp/_MEI*` folders while running.
- **The coupled tests** (`tests/test_coupled_systems.py`) need dolfinx, cmake, a C++ compiler and SUNDIALS with `libcuflynx.coupling`. They pass in the CVS-ANS session's fenicsx environment.

## 5. Other sessions working nearby (coordinate by message; never edit their worktrees)

- **"git-projects-2b":** circulatory_autogen #536/#538/#549 and library #4, #7 and #8. It will repoint CA's tests at a pinned library commit.
- **"CVS-ANS libcellml C code generation with API integration":** circulatory_autogen #520 and library #5 (external coupling, FV1D_solver). It noted that `varicosity/NEexchange_v01` must mirror any change to `varicosity/sympathetic_monolithic_v01`.
- **"calibration-workflow-implementation":** probably the calibration-workflow feature (libcuflynx #541 is `feat/calibration-workflow`).

## 6. In flight when this was written

An agent was changing the library on `devel_module_tests`:
1. NE release rates `k_NE` 0.1 /s and `k_NET` 0.03 /s (SN_full's fitted values; were 200 and 100 /s);
2. the Hernández-Cruz 1997 RyR's starting occupancies computed in the model at the starting Ca (Python/C-generation known issue added).

Check `git log` for those two commits. If they're missing, the work didn't land; redo it from this description.
