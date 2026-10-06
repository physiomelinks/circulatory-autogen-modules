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
    <module_type>/                  a directory with versions/ is a module_type (it may also sit directly in modules/)
      <module_type>.html            report across the versions (generated)
      <nested module_type>/         a module_type that only exists within this one, laid out the same way
                                    (its own versions/, and possibly module_types nested in it)
      versions/
        <version>/                  version == the config entry's module_subtype
          <module_type>_<version>_modules.cellml         its CellML component
          <module_type>_<version>_modules_config.json    exactly one config entry, with "default_instance"
          <module_type>_<version>_units.cellml           the units it uses (and what they depend on)
          <module_type>_<version>_verification_config.json   what the checks run (see "The version spec")
          <module_type>_<version>_tests.yaml             review and record-keeping (see "The version spec")
          <module_type>_<version>_references.bib         BibTeX entries its parameters and its required_citations cite
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
              source_figures/                            screenshots of the publication figures/tables the data came from
                source_figures.json                      [{"file", "source", "caption"}]; required for data from a paper or book
```

The top level of `modules/` holds categories and one module_type, `heart`:

| Directory | Holds |
|---|---|
| `vessels/compartments`, `vessels/junctions`, `vessels/terminals`, `vessels/microvasculature`, `vessels/properties` | 0D vessel segments, junctions, terminal beds, microvascular networks, wall material laws |
| `heart` (a module_type) | the heart's versions (monolithic: `vp`, `vp_Ca`, `vp_wCont`, `vp_devel` and the former alternative hearts, see below; supermodule: `Argus2026_v01`), with the nested module_types `cardiac_clock`, `chamber` and `valve` |
| `cell`, `cell/ion_channels` | cell mechanisms, one module_type each (`Ca_handling`, `membrane_potential`, `ion_concentrations`, `reversal_potentials`, `neurotransmitter_release`; every channel, pump and exchanger in `cell/ion_channels`); `cell/neuron` is the module_type neuron, with its anatomical parts (soma, axon, varicosity) nested in it |
| `respiratory` | lungs, gas exchange and transport |
| `control` | baroreflex, autonomic control, effectors, observers, PID control |
| `boundary_conditions` | inlet/outlet pressures and flows, generators, stimuli |
| `coupling` | 0D-1D coupling, volume sums |
| `transport` | tissue diffusion |
| `organs` | organ-level volume control |
| `benchmarks` | test problems (Lotka-Volterra, Van der Pol, FitzHugh-Nagumo, two minima, delay) |

`modules/poiseuille_transport/` and `modules/system/` are another session's work in the old layout and
are excluded from the schema (`directory_schema.json`, `excluded`) until they move.

## Nested module_types

**A module_type may contain the module_types that only exist within it.** Any directory with a
`versions/` directory is a module_type. Its other subdirectories that have `versions/` are **nested
module_types**, which may nest again. A nested module_type is a full module_type: its own versions,
instances, tests and reports, and it is named in configs and module arrays by its name alone
(`"module_type": "soma"`), never by its path. Only where it sits on disk says that it belongs to its
parent.

```
modules/
  heart/                                  module_type heart (directly in modules/)
    versions/vp/ vp_wCont/ vp_wCont_ASD/ ... Argus2026_v01/
    cardiac_clock/versions/Liang2009_v01/ controlled/   nested: used only in heart
    chamber/versions/vv/                         nested: used only in heart
    valve/versions/pp/ pp_linear/ pp_rmod/       nested: used only in heart
  cell/                                   category (there is no module_type "cell")
    neuron/                               module_type neuron
      versions/sympathetic/
      axon/versions/sympathetic_monolithic_v01/      nested: an anatomical part of the neuron
      soma/versions/sympathetic/ sympathetic_monolithic_v01/        nested in neuron
      varicosity/versions/sympathetic/ sympathetic_monolithic_v01/  nested in neuron
    Ca_handling/versions/SN_soma_Argus2026_v01/ SN_varicosity_Argus2026_v01/ cardiomyocyte_Paci2013_v01/
    membrane_potential/  ion_concentrations/  reversal_potentials/  neurotransmitter_release/
    ion_channels/                         category: channels, pumps and exchangers
      i_Na/  i_CaL/  i_NaK/  ...
```

The soma and varicosity supermodules are built from the mechanism types (`Ca_handling`,
`membrane_potential`, ...) and the channels; those are not nested in them, because every cell has
them (see "Placement and naming").

**Categories remain only for grouping where no module_type of that meaning exists**: `cell`,
`cell/ion_channels`, `vessels/...`, `control`, ... A category is never put
beside a module_type of the same meaning: there is no category `cell/neurons` holding a module_type
`neuron`, and no category `cardiac` holding a module_type `heart`. The module_type holds its parts
itself.

**Used only within its parent.** A nested module_type is used only inside its parent: every
supermodule submodule and every test-harness record (`harness.module_array` in a verification config)
that names it belongs to a version of the parent, of a module_type nested (at any depth) in the
parent, or of the nested module_type itself. `tests/test_structure.py` checks this. System models are
exempt: they may wire a parent's parts explicitly, e.g. `system_models/cellular/SN_simple_flat` (the
soma and varicosity parts, wired by hand) and the closed-loop models (the heart as cardiac clock,
chamber and valve records).

The category of a nested module_type is the category path above its outermost enclosing module_type
(`cell` for neuron, soma and varicosity; none for heart and its nested module_types, which form
their own group, `heart`, in the site index).

Selecting a module_type by name (`--module neuron`, `MODULE=heart`) selects it **and the module_types
nested in it**. A path under `modules/` selects everything below it (`--module cell/neuron`: neuron,
axon, soma and varicosity).

## Placement and naming

The rules of 2026-10-06 (approved on `reviews/library_restructure_proposal.html`, applied by
`tools/restructure_modules.py`):

1. **Name by mechanism.** A module type is named after the mechanism it models (`Ca_handling`,
   `membrane_potential`, `ion_concentrations`, `neurotransmitter_release`, `i_Na`, ...), never after
   the cell, compartment, system or author it was built for. Type names never carry surnames or years.
2. **Placement.** A type sits in the most general category where its mechanism is meaningful: Ca
   handling exists in every cell, so `cell/Ca_handling`; ion channels, pumps and exchangers are
   `cell/ion_channels/<channel>`.
3. **Nesting.** Nest a type only when it is an anatomical part of its parent:
   `cell/neuron/{soma,axon,varicosity}`, `heart/{chamber,valve,cardiac_clock}`. Parts are
   supermodules whose versions are cell or organ types (`soma/sympathetic`). No directory name
   appears twice in the tree, and no two sibling categories or types mean the same thing (one
   `heart/`, no `cardiac/`).
4. **Version or instance: compare the equations.**
   - **Different maths is a version**, named `<system>_<compartment>_<Source><Year>_vXX`
     (`Ca_handling/versions/SN_soma_Argus2026_v01`), or `<system>_<Source><Year>_vXX` without a
     compartment (`membrane_potential/versions/cardiomyocyte_Paci2013_v01`), or `<Source><Year>_vXX`
     where the system adds nothing (ion channels: `i_Na/versions/Tao2011_v01`). `Source` is the first
     author's surname in ASCII CamelCase (`HernandezCruz1997`); the owner's unpublished or modified
     models use `ArgusUNPUBLISHED`; an SI-unit duplicate adds `_SI`; a superseded formulation ends
     `_OLD`.
   - **Same maths with different parameter values is an instance** of the existing version
     (`versions/<v>/instances/<instance>`).
   - Vessel-network modules keep their port prefix first (`vp_wCont_ASD`), because libcuflynx reads
     a vessel's module subtype's first two letters as its inlet/outlet port kinds (see "Versions").
5. **One mechanism across cell types.** Modules of the same mechanism built for different cells or
   compartments join one type as versions or instances (rule 4). Versions of one type expose
   compatible ports where the mechanism allows; where they cannot (e.g. a bond-graph pump beside
   current-law pumps), the difference is allowed and documented (below).
6. **required_citations.** Every version's `modules_config.json` entry has
   `"required_citations": [bib keys]`: the papers whose *model* the version's equations come from
   (not the per-parameter value sources, which stay in the CSV `data_reference`). Each key is in the
   version's `references.bib` and follows `<surname><year><first significant title word>`, all
   lowercase (`belluzzi1986quantitative`). A supermodule lists its own; its report and its OMEX
   archives add its submodules'. See "Citations" below.

**Version or instance, in short.** Before adding a version, compare its equations with the type's
existing versions. Identical equations (same variables, same relations) with other parameter values
are an instance of that version, even when they were built for another cell or organ; any change to
an equation, a state or a port's variables is a version.

**Ports may differ between versions of one type.** Versions of a type usually share their ports, so
one can replace another in a model. When the mechanism does not allow it (the Paci and Tao membrane
equations take one port per named current, the SN versions one summed membrane-current port; the
bond-graph Na/K-ATPase cycle has chemostat ports, the current-law pumps a current port), the versions
keep their own ports. This is allowed, and documented: each version's notes (tests.yaml `notes`) end
with the type's port table, "Ports of the versions of <type>", listing every (direction, port_type)
and the variables each version gives it (written by `tools/restructure_modules.py`; update it when a
version's ports change). No structure check compares ports across versions.

**nn.** No version name starts with `nn` (no boundary condition): a former `nn_<x>` is `<x>`
(`inlet_flow/constant`, `cardiac_clock/controlled`) and a former plain `nn` is named by its source.
The exceptions are `coupling/FV1D_vessel/nn` and `coupling/FV1D_volume_sum/nn`, whose pairs
libcuflynx writes itself in 1D-coupled models. `tests/test_structure.py` checks this.

**Superseded.** The earlier "most specific place covering its uses" rule (which nested the SN parts
in soma and varicosity) and the `cell/cardiomyocytes` category. `tools/restructure_modules_renames.json`
lists every old (module_type, version) and its new pair; `cam_testing.library.legacy_renames()` reads
it, so circulatory_autogen's own models and older module arrays still resolve.

## Versions

A version is one piece of math for a module_type: one CellML component and one config entry. **The
version's name is the config entry's `module_subtype`**, so a model picks a version with
`"module_subtype"`.
- **Names** follow "Placement and naming", rule 4: `<Source><Year>_vXX` (`Paci2013_v01`,
  `Tao2011_v01`, `ArgusUNPUBLISHED_v01`), with `<system>_` and `<compartment>_` in front where they
  tell versions apart (`SN_soma_Argus2026_v01`, `cardiomyocyte_Paci2013_v01`), `_SI` and `_OLD`
  suffixes. Vessel versions keep their descriptive subtypes (`vp`, `pp_nonlinear`, `vv_noI_SI`).
- **The port-letter constraint.** libcuflynx reads a vessel's module_subtype's first two letters as
  the inlet and outlet port kinds (`pp`, `pv`, `vp`, `vv`), and only a vessel's (a module with vessel
  ports whose subtype starts with one of the four: `libcuflynx.utilities.vessel_bc`). So a vessel
  version name starts with its port letters, and any other version may have any name (`lv_...`,
  `SN_soma_...`); it should still not start with `pp`, `pv`, `vp` or `vv` unless it follows that
  convention.
- **A module_type folded into another keeps its port letters first.** When a former module_type
  becomes a version of another, the version is named `<old module_subtype>_<old module_type suffix>`,
  so its first two letters (and so libcuflynx's connection rule) are unchanged. The alternative hearts
  became versions of `heart` this way:

  | Former (module_type, module_subtype) | Now (`heart`, version) |
  |---|---|
  | `heart_ASD`, `vp_wCont` | `vp_wCont_ASD` |
  | `heart_LVprop`, `vp` | `vp_LVprop` |
  | `heart_new_valve`, `vp` | `vp_new_valve` |
  | `heart_nonstiff`, `vp_wCont` | `vp_wCont_nonstiff` |
  | `heart_simple`, `vp` | `vp_simple` |
  | `heart_simple_2`, `vp` | `vp_simple_2` (now `vp_simple_SI`) |
  | `heart_simple_OLD`, `vp` | `vp_simple_OLD` |

  Their CellML components keep their names (`heart_new_valve`, `heart_nonstiff`, ...), so generated
  models are unchanged. `cam_testing.library.legacy_renames()` maps circulatory_autogen's old pairs
  to these (`HEART_VERSION_RENAMES`).
- **Shared components.** libcuflynx finds CellML components by name across the whole library, so two
  versions can't define components with the same name. A component that several old entries shared
  was copied into each version and renamed `<component>__<module_type>_<version>` (in the CellML and
  in the config's `component_type`). The generated model names every component after its vessel, so
  results are unchanged.

**Licence and creator.** Every config entry has two record-keeping keys, after `default_instance`:
- `"licence"`: the SPDX id of the licence the version is published under. It is one of `CC0-1.0`
  (the default, chosen as the most open), `CC-BY-4.0`, `MIT`, `Apache-2.0` and `0BSD`
  (`cam_testing.library.LICENCES`). The version page links it to its SPDX page.
- `"creator"`: a list of the names of the people who created the version, e.g. `["A. Author"]`. It is `[]`
  until the creators are given in review, and the version page shows "creator: not set (to be given in
  review)" until then.

`tests/test_structure.py` checks both keys in every entry. libcuflynx's config schema allows extra
keys (these, and the citation keys below), and both config formats carry them unchanged (`library.normalise_config_entry`,
`library.to_phlynx_entry`). The tools that write new config entries (`tools/restructure_to_versions.py`,
`tools/import_from_libcuflynx.py`, `tools/split_module.py`) add them with `library.with_record_keys`.

**Citations.** Every config entry names the papers its model is built on, after `creator`:
- `"required_citations"`: a non-empty list of BibTeX keys, each with its entry in the version's
  `<module_type>_<version>_references.bib`. These are the papers whose *model* (equations) the
  version implements, to be cited by anyone who uses it; the papers its parameter values come from
  are cited per parameter in the instance CSVs' `data_reference` (`<key>; <note>`) and need not be
  listed here.
- `"required_citations_uncertain"`: the keys among them that are a best guess still to confirm
  (e.g. a primary source not checked equation by equation). Left out when there are none. The
  version page shows "(uncertain)" next to them.
- Models with no publication behind them cite placeholders, `@unpublished` entries whose note says
  they are placeholders: `argus2026sympathetic` (the planned sympathetic-neuron paper),
  `argus2026unpublished` (the owner's unpublished or modified models) and `gee2026unpublished` (the
  Gee and Argus respiratory-control work). Generic mechanisms with no model paper (constant
  boundary conditions, sums, observers, generic vessels) cite the software, `argus2026circulatory`
  (circulatory_autogen / libcuflynx and this library). Replace a placeholder with the real key when
  the paper exists.
- A supermodule lists the citations of its own (often a placeholder); its version page adds "From
  its submodules:" (each linked to the submodule's page), and its per-instance OMEX archives carry
  the submodules' configs too.

The version page shows the required citations under the description, each linked to its entry in
the page's References. The initial lists come from the proposal approved on 2026-10-06
(`tools/required_citations.py`, which applies its rules and adds the missing .bib entries; re-running
it fills only entries without citations). A new version needs its own: the tools that write new config
entries leave `required_citations` out, and the structure test fails until it is given.

**BibTeX keys.** Every key is `<first author's surname><year><first significant title word>`, all
lowercase ASCII: `belluzzi1986quantitative`, `vanderpol1926relaxation` (particles join the surname),
`hernandezcruz1997ca` (accents and hyphens dropped); stopwords (a, an, the, on, of, in, for, to, and,
at, by, with, from, is, are) are skipped (`cam_testing.bib.convention_key`). One key per paper,
library-wide: the same paper has the same key in every version's .bib. Two different papers with the
same key take further title words (`<surname><year><first word><second word>`). Undated web pages have no year in
their key and are listed in `UNDATED_BIB_KEYS` in `tests/test_structure.py`. The library was rekeyed
to this convention on 2026-10-06 by `tools/rekey_bib.py` (old -> new in `tools/bib_rekey_map.json`;
four papers that had two keys each now have one). `tests/test_structure.py` checks the
convention, one key per DOI, and every version's required citations.

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
  `benchmarks/Lotka_Volterra/versions/Lotka1925_v01/instances/hudson_bay_lynx_hare/`. Calibration writes
  `<instance>_calibrated_parameters.csv` and `<instance>_calibration.json` beside them.
- **Recalibrating, and applying a calibration.** A calibration to data is always a libcuflynx run on
  a committed obs_data, never a hand calculation, so anyone can change the data and redo it:
  1. edit `<instance>_obs_data.json` (values, std, protocol) or `<instance>_params_for_id.csv`;
  2. run the instance's calibration test,
     `pytest tests/test_modules.py --component <module_type>/<version> -k 'calibrate and <instance>' --include-unreviewed`
     (it rewrites `<instance>_calibrated_parameters.csv` and `<instance>_calibration.json`);
  3. to make the calibrated values the defaults, run
     `python -m cam_testing.calibrate apply <module_type>/<version> <instance>` (`--dry-run` lists the changes first).

  `apply` copies the parameters of `<instance>_params_for_id.csv`, at their calibrated values, into:
  - the version's default instance;
  - every supermodule version that uses the version's default instance, recursively (row
    `<var>_<submodule path>`, e.g. soma/sympathetic `rho_M_i_M`);
  - the system models that use the version or such a supermodule at its default instance (row
    `<var>_<vessel>[_<submodule path>]`, e.g. `rho_M_SN_soma_i_M`);
  - the monolithic counterparts named by a supermodule's `supermodule.equivalent` output_map (row
    `<var>_<vessel>` in the `reproduces` system model, e.g. `rho_M_soma_SN`, and `<var>` in that
    vessel's version's default instance).

  Only rows that already exist are rewritten. Each reference becomes `Calibrated (libcuflynx) to
  instances/<instance>/<instance>_obs_data.json: <note> (was <old value>; applied <date>)`. The
  `<note>` is the obs_data's optional `"calibration_note"`; an optional `"reference_key"` (a BibTeX key)
  prefixes the reference and marks the row sourced.

### Source figures

**Rule: data extracted from a paper or a book is shown beside its source.** When an instance's
validation or calibration data were read from a publication (digitised from a figure, copied from
a table or the text), the instance keeps a screenshot of that figure, table or passage, and the
report shows it next to the instance's validation and calibration plots, so a reader can judge the
extraction and the fit against the original.

- **Where:** `instances/<instance>/source_figures/`: the images (PNG, cropped to the figure or
  table, readable, e.g. rendered with `pdftoppm -r 200` and cropped) and `source_figures.json`:
  ```json
  [{"file": "davis2020_fig3c.png", "source": "davis2020downregulation; Fig. 3C",
    "caption": "Wistar XE-991-sensitive I_M, hold -25 mV, step to -55 mV (digitised points in the obs_data)"}]
  ```
  `source` starts with the BibTeX key (in the version's references.bib) and names the figure,
  table or page; `caption` says what was extracted from it.
- **When it is required:** for every instance whose validation entry (verification_config
  `validation.<instance>.baseline` / `.calibrate`) has a `source` and a `source_kind` of
  `"publication"`, the default when a source is given. Data that is not extracted from a
  publication declares `"source_kind": "dataset"` (a data file or database used as is) or
  `"synthetic"` (generated by a model, e.g. the FitzHugh-Nagumo and TwoMinima benchmarks) and needs
  no screenshot.
- **Checked by** `tests/test_structure.py::test_publication_data_has_source_figures`. Instances
  that predate the rule and still lack screenshots are listed in `MISSING_SOURCE_FIGURES` there;
  the list only shrinks. The report marks such an instance "No source screenshot".
- **Copyright:** a cropped figure or table with its citation, kept for checking the extracted
  values, is the intent; don't include whole pages or more than the extraction needs.

## How system models use versions and instances

System models live in `../system_models/<category>/<model>/` (not in `modules/`; they are not
modules or supermodules). A module-array record names a module_type, a version and an instance:

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
| `heart` `Argus2026_v01` | cardiac clock, four chambers, four valves (`heart/cardiac_clock`, `heart/chamber`, `heart/valve`) |
| `cell/neuron` `sympathetic` | soma `sympathetic`, axon `sympathetic_monolithic_v01`, varicosity `sympathetic` |
| `cell/neuron/soma` `sympathetic` | membrane (`membrane_potential`), Na/K (`ion_concentrations`), Ca handling (`Ca_handling`), each version `SN_soma_Argus2026_v01`, and 15 `cell/ion_channels` versions |
| `cell/neuron/varicosity` `sympathetic` | membrane, Ca handling, NE release (`neurotransmitter_release`), each version `SN_varicosity_Argus2026_v01`, and channels |

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
parameters. Every version runs them, component or supermodule:

| Test | Checks |
|---|---|
| `run_test` | generates the version alone (or in its spec's `harness` network) and simulates it; every output finite |
| `verification_test_invariants` | the spec's invariants hold at the run parameters |
| `verification_test_BC` | each boundary condition / constant swept over its tested range; runs finite, invariants hold. Points outside a valid range are skipped; failures are classified and a solver setting recommended (below) |
| `verification_test_timestep` | fixed-step integration converges at the scheme's order and matches CVODE |
| `stability_test` | the solver/tolerance/step matrix; the declared-supported configurations work |
| `phlynx_export_test`, `cuflynx_simulate_test`, `phlynx_equivalence_test` | the PhLynx → CUFLynx pipeline (`tests/test_phlynx.py`) |
| `supermodule_structure_test`, `supermodule_equivalence_test` | supermodule versions, in addition to the above |

**Supermodule versions** run the same verification and validation tests as a component. The
supermodule is generated alone as the vessel `mod`: libcuflynx expands it into `mod_<submodule>`
(nested supermodules: `mod_<submodule>_<subsubmodule>`), and its parameters are those of the
flattened model: the supermodule's own instance, its submodules' instances, and the boundary
conditions no internal connection closes. In the spec (`run_parameters`, `bc_sweep`, `validation`)
a submodule's parameter is named as in the supermodule's instance, `<var>_<submodule path>` (e.g.
`I_in_membrane`, or `I_in_soma_membrane` in `neuron/sympathetic`); globals keep their names. The
default `outputs` are every state of the model (`mod_<submodule>/<var>`). The default BC sweep
(`bc_sweep.sweep: globals_and_bcs`) varies only what the supermodule adds: its global constants and
those open boundary conditions, plus any `bc_sweep.extra_parameters`. Each submodule's own parameters
are swept in that submodule version's own tests (soma `sympathetic` has 133 parameters; sweeping them
all again would repeat the submodules' sweeps at many times the run time). `bc_sweep.sweep: all`
sweeps every parameter of the flattened model. PhLynx has no supermodule support, so a supermodule
version's `phlynx_export_test` fails (an `expected_failures` known issue in its tests.yaml) and the
other two pipeline tests are skipped.

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

**The report's test overview has the same columns for every version**, component or supermodule,
in this order: Run, Invariants, BC sweep, Timestep, Stability, Calibration (the version's), PhLynx
export, CUFLynx simulate, PhLynx equivalence, Supermodule structure, Supermodule reproduces. The
module_type page and the status counts use the same columns. A cell shows:
- **N/A** when the test doesn't apply to the version, with the reason (hover the cell; the reasons
  are also listed under the overview). For example "not a supermodule" for the last two columns of a
  component, a C++ version's simulation tests, or a supermodule with no `supermodule.equivalent`.
- **Not run** (grey) when the test applies but has no result yet.
- Passed only when the test ran and passed. A test is never shown as passing because it didn't run.

Status rules (decided in the Lotka_Volterra review, 2026-10-01):
- An instance without obs_data has calibration **not applicable** ("no calibration data in this
  instance"); likewise an instance without baseline data has the baseline not applicable. The instance
  stays usable in models.
- The version's **Calibration** column in the report's test overview fails when **no instance** has
  calibration data, fails when any instance's calibration fails, and passes when they all pass.
- **Calibration in a supermodule.** A version with no calibration data of its own that is a submodule
  of one or more supermodule versions (e.g. `membrane_potential/SN_soma_Argus2026_v01`, which can't be validated alone) is
  calibrated as part of them. Its Calibration column shows **Pass in super** when any of those
  supermodules' calibration passes, and **Fail in super** when none does; a supermodule with no
  calibration data counts as failed. The rule applies transitively: an ion channel takes its status
  from `soma/sympathetic`, which (with no data of its own) takes it from `neuron/sympathetic`. The cell
  and the test entry link to the supermodule versions. A version with calibration data of its own keeps
  its own result. `calibration_in_supermodule: false` in tests.yaml opts a version out (the plain rule).
  In the status counts, Pass in super counts as passed and Fail in super as failed.
- A check blocked by an external release (e.g. the released CUFLynx is too old for a feature) is a
  **failure with a known issue** (`expected_failures`), not not-applicable.
- In the stability matrix a configuration that is too inaccurate at a step (e.g. explicit Euler at a
  coarse dt) is shown as not working at that step; the test fails only if a declared-supported
  configuration fails.

What the version page shows besides the tests:
- **Contents.** A line of counts at the top: states (variables with a d/dt equation in the CellML),
  algebraic variables (defined by an algebraic equation), parameters / constants (config kinds `constant`
  and `global_constant`), boundary conditions, ports, equations and instances. A supermodule adds its
  parts and nested levels, and its other counts are summed over all its parts. Each count links to its
  section and explains itself on hover; the module_type page shows them compactly per version.
- **Unit consistency** (informational, not a test column). `checks.unit_consistency` runs libcellml's unit
  checks on the version's CellML with its units file, without simulating: the Validator's units issues
  (undefined or invalid units) and the Analyser's `ANALYSER_UNITS` issues (equation sides or terms in
  units that are not equivalent, or an argument that should be dimensionless). Inputs are treated as
  constants for the check. Each issue is mapped to the equation (and the variable it defines). The result
  is `results/unit_consistency.json`; the report computes it when it is missing or older than the CellML,
  units or config. The page shows a "Unit consistency" block only when an equation fails (or libcellml
  could not check the equations, e.g. an undeclared variable), and the header gives the count.
- **Structure** (supermodule versions). A diagram of the parts (each linked to its version page) and their
  internal connections with the port types the configs connect, the modules outside that couple to a
  part (from the supermodule versions and system models that use this one, through
  `per_submodule_inputs` / `per_submodule_outputs`), and a table of the parts: version, instance, ports,
  connections, the supermodule instance's rows `<var>_<part>` that override the part, and the globals it
  shares. A collapsible tree lists every nested level. The diagram is drawn by mermaid (loaded from
  cdn.jsdelivr.net like MathJax); offline, its source text and the table remain.

The version's spec is described in the next section.

```bash
pytest tests/test_modules.py --module Lotka_Volterra          # a module_type
pytest tests/test_modules.py --module neuron                  # neuron and the module_types nested in it
pytest tests/test_modules.py --module cell                    # every module_type in a category
pytest tests/test_modules.py --component Lotka_Volterra/Lotka1925_v01    # one version (and its instances)
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
| `outputs` | verification_config | variables to log and check (`var`, or `vessel/var` for a harness neighbour); default: every `variable` of the config (a supermodule: every state, `mod_<submodule>/<var>`) |
| `run_parameters` | verification_config | `{parameter: value}`: the operating point of the run and invariant tests |
| `rest_check` | verification_config | `{sim_time, window, voltage, parameters, threshold, dt}`: a long run (e.g. 0 pA injected) that run_test reports, not a pass/fail gate: the spikes (upward crossings of `threshold` mV, default 0) of `voltage` and their rate in the last `window` s, in the message, metrics.rest_check and plots/rest.png |
| `harness` | verification_config | the test network: `module_array` rows `[name, module_subtype, module_type, inp, out, instance]` (the version under test is `mod`) and `parameters` rows `[name, units, value, source]` for the neighbours; an instance may have its own, `validation.<instance>.harness` (same form), used for that instance's model, calibration and archive instead (e.g. an experiment's set-up: `PMCA/Colegrove2000_v01` instance `wanaverbecq2003_scg` puts the pump on the soma Ca handling with a Ca leak and a load, while the version's own tests run the pump alone) |
| `invariants` | verification_config | numpy expressions (`expr`, with `description` and `applies`) that must hold |
| `bc_sweep` | verification_config | the parameter sweep: `sweep` (`all`, a component's default: every boundary condition and constant; `bcs`: the boundary conditions; `globals_and_bcs`, a supermodule's default: its globals and open boundary conditions), `factors`, `points`, `ranges` (`[min, max]` or `{min, max, points, scale}` or `{values}`), `bounds` (`{param: [min, max]}`, `null` for no limit: the valid range; values outside it are skipped), `constraints` (numpy expressions of the parameters by variable name, e.g. `B_Ca_init < B_Ca_total`; a point that breaks one is skipped), `exclude`, `extra_parameters`, `plot_output`, `rationale` |
| `timestep` | verification_config | the convergence test: `scheme`, `dts`, `t_end`, `min_order`, `tol`, `cvode_tol`, `roundoff`, `wrapped` |
| `stability` | verification_config | the solver matrix: `supported` (must work), `cvode`, `solve_ivp`, `fixed_step`, `max_step_start`, `min_step`, `time_budget`, `t_end`, `tol` |
| `parameter_ranges` | verification_config | published parameter intervals, `{parameter: [lo, hi]}` (reported as validated values; usually inside a baseline block) |
| `validation` | verification_config | per instance, `validation.<instance>.harness` (the instance's own test network, see `harness`) and `validation.<instance>.baseline` / `.calibrate`: `status`, `source`, `source_kind` (`publication` default / `dataset` / `synthetic`; see Source figures), data references relative to the version directory (`data: instances/<i>/<file>.csv`, `obs_data`, `params_for_id`), variables, targets, metric, threshold, optimiser settings |
| `supermodule` | verification_config | supermodule versions: `globals` (parameters that name no submodule) and `equivalent` (system models it must reproduce: `model`, `reproduces`, `instance`, `tol`, `solver_info`, `output_map`, `ignore`) |
| `reviewed` | tests | whether the version has been reviewed (unreviewed versions run in CI without blocking it) |
| `description` | tests | what the version is (supermodule versions) |
| `notes` | tests | notes on the version, shown in its report |
| `skip` | tests | a reason not to run the version's tests |
| `known_issues` | tests | substrings of structure problems that are known (xfail) |
| `expected_failures` | tests | `{test: reason}`: known failures, recorded as failed and xfail in pytest |
| `calibration_in_supermodule` | tests | `false`: the version's Calibration column uses the plain rule even when it has no calibration data and is a submodule of supermodule versions (default: derived from those supermodules, "Pass in super" / "Fail in super") |
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
| `<module_type>_<version>_<instance>_module_array.json`, `_model_parameters.csv` | The test network and the parameters file the model was generated from. |

CUFLynx picks members by name: the first `.json` with "obs" in its name is the obs_data, and the first `.csv` with "param" in its name is the params_for_id. So the study members come first.

An instance with no obs_data or params_for_id still loads and simulates, but CUFLynx then tries other members in those roles and shows two harmless "could not read" notes. That limitation is noted for a future change.

`tests/test_instance_omex.py` (`make omex-test`) loads every archive into a released CUFLynx and checks four things:
- CUFLynx takes the master model, the obs_data and the params_for_id;
- the run is finite;
- the outputs match libcuflynx's model of the same instance within 1e-6;
- the version report shows this as the instance's CUFLynx tick.

## Directory schema

`directory_schema.json` describes every level: `modules_root` (categories and top-level module_types),
`category`, `module_type` (its `versions` and nested module_types), `versions`,
`version`, `instances`, `instance`, `risk`, and `system_models_root` / `system_category` /
`system_model`. For each level it gives the name pattern, the required and optional files (with
`{module_type}`, `{version}`, `{instance}`, `{model}` placeholders) and the allowed subdirectories.
`spec_keys` lists which spec keys go in which file. Its `rules` list the cross-level rules:
- the version is the module_subtype, with one config entry;
- the spec keys each sit in their own file (tests.yaml or verification_config.json);
- the default instance exists;
- an instance is its obs_data_name;
- no directory name repeats (nested module_types included);
- no category is named like a module_type;
- a nested module_type is used only within its parent;
- where a module_type goes (placement and nesting);
- how versions are named (by source; vessels' port letters first; no `nn`);
- versions of one type may differ in their ports (documented in their notes);
- every system-model record names an existing version and instance.

`tests/test_structure.py` walks `modules/` and `system_models/` and checks all of this. It also checks
each version's CellML, config and units, and validates the JSON files against libcuflynx's schemas.

## Adding things

- **A new version** of an existing module_type: add
  `versions/<source>_vXX/` with the six files, one config entry (`module_subtype` = the directory name,
  `"default_instance": "default"`, `"licence"` and `"creator"`), a component name no other version uses, and
  `instances/default/default_parameters.csv`. Then run `make structure` and
  `pytest tests/test_modules.py --component <module_type>/<source>_vXX`.
- **A data instance**: add `instances/<name>/` with `<name>_parameters.csv`, the obs_data files
  (`"obs_data_name": "<name>"`), `<name>_params_for_id.csv`, the raw data and `SOURCES.md`. Then add
  `validation.<name>` to the version's verification_config.json.
- **A new module_type**: first check that the mechanism has no type yet (a new cell's Ca handling is
  a version of `Ca_handling`, not a new type). Name it by its mechanism and put it in the most general
  category where the mechanism is meaningful; nest it only when it is an anatomical part of its
  parent (see "Placement and naming"). Then run `make manifests`.

`tools/restructure_to_versions.py` moved the old per-module layout into this one, driven by
`tools/restructure_map.yaml` (made by `tools/restructure_map.py`). The later move to nested
module_types (the category `cell/neurons` and `cardiac` replaced by the module_types `cell/neuron` and
`heart`, the alternative hearts made versions of `heart`) is not in that map. The 2026-10-06
restructure (mechanism types, source-named versions, no `nn`) is `tools/restructure_modules.py`, with
its old -> new map in `tools/restructure_modules_renames.json`.
