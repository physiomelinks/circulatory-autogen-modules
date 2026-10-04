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

The top level of `modules/` holds categories and one module_type, `heart`:

| Directory | Holds |
|---|---|
| `vessels/compartments`, `vessels/junctions`, `vessels/terminals`, `vessels/microvasculature`, `vessels/properties` | 0D vessel segments, algebraic flow nodes, terminal beds, microvascular networks, wall material laws |
| `haemodynamics/lumped_constitutive`, `haemodynamics/lumped_vessel` | the compliance, resistance, inertance and vessel_volume modules the vessels' `_lumped` supermodules are built from, and the empty lumped vessels (templates); see "Lumped vessels" |
| `heart` (a module_type) | the heart's versions (monolithic: `vp`, `vp_Ca`, `vp_wCont`, `vp_devel` and the former alternative hearts, see below; supermodule: `Argus2026_v01`), with the nested module_types `cardiac_clock`, `chamber` and `valve` |
| `cell`, `cell/ion_channels`, `cell/cardiomyocytes` | cell models and their parts; `cell/neuron` is the module_type neuron, with its parts nested in it |
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
instances, tests and reports, and it is named in configs and vessel arrays by its name alone
(`"module_type": "soma"`), never by its path. Only where it sits on disk says that it belongs to its
parent.

```
modules/
  heart/                                  module_type heart (directly in modules/)
    versions/vp/ vp_wCont/ vp_wCont_ASD/ ... Argus2026_v01/
    cardiac_clock/versions/nn/ nn_controlled/    nested: used only in heart
    chamber/versions/vv/                         nested: used only in heart
    valve/versions/pp/ pp_linear/ pp_rmod/       nested: used only in heart
  cell/                                   category (there is no module_type "cell")
    neuron/                               module_type neuron
      versions/sympathetic/
      NE_release_Tao_2011/                nested: a part of the Tao 2011 neuron
      sympathetic_neuron_membrane_voltage/    nested: a part of the Tao 2011 neuron
      axon/versions/sympathetic_monolithic_v01/
      soma/                               nested in neuron
        versions/sympathetic/ sympathetic_monolithic_v01/
        SN_membrane_soma/  SN_Na_K_concentrations_soma/  SN_Ca_handling_soma/   nested in soma
      varicosity/                         nested in neuron
        versions/sympathetic/ sympathetic_monolithic_v01/
        SN_varicosity_membrane/  SN_Ca_handling_varicosity/  SN_NE_release/   nested in varicosity
    ion_channels/                         category: channels shared by neurons and cardiomyocytes
      i_Na/  i_CaL/  ...
    cardiomyocytes/                       category: there is no module_type "cardiomyocyte"
      myocyte_membrane_voltage/  Ca_dynamics_Paci_2013/  ...
```

**Categories remain only for grouping where no module_type of that meaning exists**: `cell`,
`cell/ion_channels`, `cell/cardiomyocytes`, `vessels/...`, `control`, ... A category is never put
beside a module_type of the same meaning: there is no category `cell/neurons` holding a module_type
`neuron`, and no category `cardiac` holding a module_type `heart`. The module_type holds its parts
itself.

**Used only within its parent.** A nested module_type is used only inside its parent: every
supermodule submodule and every test-harness record (`harness.vessel_array` in a verification config)
that names it belongs to a version of the parent, of a module_type nested (at any depth) in the
parent, or of the nested module_type itself. `tests/test_structure.py` checks this. System models are
exempt: they may wire a parent's parts explicitly, e.g. `system_models/cellular/SN_simple_flat` (the
soma and varicosity parts, wired by hand) and the closed-loop models (the heart as cardiac clock,
chamber and valve records).

The category of a nested module_type is the category path above its outermost enclosing module_type
(`cell` for neuron, soma and `SN_membrane_soma`; none for heart and its nested module_types, which form
their own group, `heart`, in the site index).

Selecting a module_type by name (`--module neuron`, `MODULE=heart`) selects it **and the module_types
nested in it**. A path under `modules/` selects everything below it (`--module cell/neuron/soma`:
soma and its three parts).

## Placement rule

A module_type goes in **the most specific place that covers every context it is used in**, not as
shallow as possible:
- **Inside the module_type it is only used within**, as a nested module_type, when there is one:
  - `soma`, `axon` and `varicosity` are only parts of the neuron, so they are in `cell/neuron/`.
  - The membrane, Na/K and Ca-handling parts of the soma (`SN_membrane_soma`,
    `SN_Na_K_concentrations_soma`, `SN_Ca_handling_soma`) are only parts of the soma, so they are in
    `cell/neuron/soma/`, not in `cell/neuron/`.
  - `SN_varicosity_membrane`, `SN_Ca_handling_varicosity` and `SN_NE_release` are only parts of the
    varicosity, so they are in `cell/neuron/varicosity/`.
  - `NE_release_Tao_2011` and `sympathetic_neuron_membrane_voltage` are parts of the Tao 2011 neuron,
    neither soma- nor varicosity-specific, so they are directly in `cell/neuron/`.
  - `cardiac_clock`, `chamber` and `valve` are only used to build hearts, so they are in `heart/`.
- **Otherwise in the most specific category covering all its uses**:
  - Every ion channel is in `cell/ion_channels/`, one module_type per channel kind (`i_Na`, `i_CaL`,
    ...), because channels of one kind are used across cell types (the soma's channels by the neuron,
    the Paci channels by the cardiomyocyte). They are not nested in `soma` even though the soma
    supermodule uses them. The cell type a version was built for is in its config entry's `"notes"`,
    e.g. `"Built for: sympathetic neuron (Tao et al. 2011)"`.
  - A sarcoplasmic reticulum (or a generic Ca store) would belong to several cell types, so it would
    go in `cell/`, as `NKE_pump` does today.
  - `myocyte_membrane_voltage` is only used by the Paci cardiomyocyte, which is not a module_type, so it
    is in the category `cell/cardiomyocytes/`.
- **One module_type, not several.** A variant of a module_type with the same role and ports is a
  version of it, not a module_type beside it: the alternative hearts `heart_ASD`, `heart_LVprop`, ...
  are versions of `heart` (see "Versions").

**No directory name appears twice** anywhere in `modules/`, among the categories, module_types and
nested module_types, so a category is never named like a module_type. The generic names (`versions`,
`instances`, `risk`, `plots`, `results`), version names (which repeat across module_types: `nn`, `vp`,
...) and instance names (which name data sets, and `default`) are exempt.

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
  | `heart_simple_2`, `vp` | `vp_simple_2` |
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
keys, and both config formats carry them unchanged (`library.normalise_config_entry`,
`library.to_phlynx_entry`). The tools that write new config entries (`tools/restructure_modules.py`,
`tools/import_from_libcuflynx.py`, `tools/split_module.py`) add them with `library.with_record_keys`.

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

## Connecting vessels

There are no junction module types. Any number of vessels can meet at a point (a **node**), on
either side of it, because every vessel's **compliant end** sums over its node: a `vessel_port`
whose flow is an input and whose pressure is an output has `"multi_port": ["sum", "True"]`. The
BC letters say which ends those are: `v` is a compliant end (it takes the summed flow and sets
the pressure), `p` is not (it takes the node's pressure and gives its own flow).

| BC type | inlet | outlet |
|---|---|---|
| `vp` | sums | plain |
| `pv` | plain | sums |
| `vv` | sums | sums |
| `pp` | plain | plain |

Each node has exactly one summing end, its **owner**. libcuflynx (`generators/port_nodes.py`)
gives the owner the signed sum of every other end's flow (flows into the node add, flows out of it
subtract) and maps its pressure to every end. So:

- **many vessels into one** (the old `Min`): the downstream vessel's inlet owns the node (`vp` or
  `vv`), and the upstream vessels end in `p` (`vp` or `pp`);
- **one vessel into many** (the old `Nout`): the upstream vessel's outlet owns it (`pv` or `vv`),
  and the downstream vessels start with `p` (`pv` or `pp`);
- **inflows and other outflows at one node** (which the junction types could not express): one
  owner on either side, everything else `p` at that node.

Two summing ends at one node (for example three `pv` vessels merging) would be two compliances
for one pressure, and generation stops naming the node; so does a node of three or more ends
where none sums. Where only two modules meet the connection is one-to-one, so a vessel whose port
can sum costs nothing there. A summing port left unconnected is a boundary condition set by the
parameters file, like any other open port. A 1D vessel can only meet one 0D module at a node.

`python -m libcuflynx.utilities.junction_migration <vessel_array>...` converts a vessel array that
still names the removed junction types (and the microvasculature `<vessel>_Min/_Nout/_Minlet/
_Noutlet/_MinNout` and `artery_inlet/_outlet`) to these vessels, with its parameters file.

## Supermodule versions

A supermodule version is a version whose config entry has `"module_format": "supermodule"`. It is
built from other versions (its `submodules`, each with its own `module_type`, `module_subtype` and
`instance`), with internal connections only. It replaces the monolithic version of the same
module_type:

| Supermodule version | Built from |
|---|---|
| `heart` `Argus2026_v01` | cardiac clock, four chambers, four valves (`heart/cardiac_clock`, `heart/chamber`, `heart/valve`) |
| `cell/neuron` `sympathetic` | soma `sympathetic`, axon `sympathetic_monolithic_v01`, varicosity `sympathetic` |
| `cell/neuron/soma` `sympathetic` | membrane, Na/K, Ca handling (nested in soma) and 15 `cell/ion_channels` versions |
| `cell/neuron/varicosity` `sympathetic` | membrane, Ca handling, NE release (nested in varicosity) and channels |

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

## Lumped vessels

Every 0D vessel version that splits into a compliance, a resistance and an inertance has a
supermodule twin, `<version>_lumped` beside it (`arterial` `vp` -> `arterial` `vp_lumped`), built
from the modules in `haemodynamics/lumped_constitutive`:

| module_type | Versions | What it is |
|---|---|---|
| `compliance` | `vv_linear`, `vv_linear_novisco`, `vv_linear_geometric`(`_novisco`), `vv_nonlinear`, `vv_nonlinear_geometric`(`_visco`), `vv_linear_controlled`(`_novisco`, `_stressed`), `vv_linear_volume_flux` | the pressure-volume law P(A): a compliant node (`vv`, both ports sum) |
| `resistance` | `<form>_<law>`: form `pv` (u_out = u_in - f R v, before an inertance), `vp` (u_in = u_out + f R v, after one), `pp` (v = (u_in - u_out)/(f R), no inertance); law `linear`, `geometric` (Poiseuille from r_0), `nonlinear` (Poiseuille from the vessel's volume), `volume_scaled`, `controlled`, `controlled_share`, `controlled_additive` | the pressure drop |
| `inertance` | `pp_linear`, `pp_geometric` (I from r_0, l, rho, with the hydrostatic term) | the flow's inertia |
| `vessel_volume` | `nn`, `nn_geometric` | the whole vessel's volume (summed over its compliances) and, `nn_geometric`, its r_0, l and radius r: for hosts that read a vessel's volume or radius (volume sums, material-property modules) and resistance laws that need the whole volume |

Every compliance, resistance and inertance has a dimensionless `fraction`, the share of the vessel it
stands for: a `vv` vessel is two half compliances (fraction 0.5) around R and I; a `pp` vessel is a
compliance between two halves of R and I. The order follows the BC letters: `vp` is C -> R -> I,
`pv` is I -> R -> C, `vv` is C -> R -> I -> C, `pp` is I -> R -> C -> R -> I, and the vessels without
inertance put a flow-form resistance where R and I were.

A `_lumped` version's config entry has the keys that let a model use it exactly as it used the
monolithic version, with the same names:

- **`routes`** say which submodule a host connects to, by port type: `{"inputs": {"vessel_port": "C"},
  "outputs": {"vessel_port": "I", "volume_port": "V"}}`. A model keeps its one record for the vessel
  (`{"name": "aortic_root", "module_type": "arterial", "module_subtype": "vp_lumped", ...}`) with its
  usual `inp_instances`/`out_instances`, and libcuflynx links each neighbour to the routed submodule.
  Ports of the monolithic version that no submodule offers are listed in the version's `notes`.
- **`shared_parameters`** keep the monolithic version's parameter names. An entry is a variable
  (`"l"`: the vessel's `l`, in every submodule that has one) or `{"name": "C_T", "variable": "C",
  "submodules": ["C"]}` (the terminal's `C_T` is its compliance's `C`). Each is one parameter of the
  model, `<name>_<vessel>` (`C_T_systemic_T`, as before), mapped to every submodule that takes it,
  so parameters files, params_for_id and calibration keep their names. What the vessel splits
  (`fraction_R_p`, a terminal's `q_C_init_C = q_init - q_us`) is `<var>_<vessel>_<submodule>`.
- **`outputs`** keep the monolithic version's output names: `{"u": "C_p/u", "v": "I/v", "q": "V/q"}`
  gives the model `aortic_root/u` (computed from the submodule's variable in a component named after
  the vessel), so obs_data, prediction variables and plots keep working.
- **`replaces`** says which monolithic version it is the twin of and how to compute the split
  parameters from its parameters (`{"module_subtype": "pp_RICRI", "parameters": {"fraction_R_d":
  "(1.0 - frac_R_T_1_of_R_T)"}}`). `python -m libcuflynx.utilities.lumped_migration --library
  modules <vessel_array> --parameters <parameters.csv>` moves a model onto the twins with it: the
  records' `module_subtype`, and only the computed rows of the parameters file, change.
- its default instance carries the monolithic default instance's values (TODOs filled with
  representative values), and its spec's `supermodule.equivalent_version` checks that it reproduces
  the monolithic version (`supermodule_version_equivalence_test`). Each output's difference is
  relative to that output's own scale (its largest magnitude, or its range), so a microvascular flow
  of 1e-11 m^3/s is judged against 1e-11, not against the pressures. A state whose scale the
  solver's `atol` does not resolve (atol above 1% of it) fails the check rather than passing on
  noise; the microvascular equivalences run at atol 1e-20 for this.

The monolithic versions stay in the library, and models that name them still build, but their config
entries say `"available_to_phlynx": false`, so the PhLynx manifests (`make manifests`) offer the
`_lumped` versions instead. `tools/lumped/` generates all of it: `python -m tools.lumped.write` rewrites
the constitutive modules, the `_lumped` versions and the templates (edit `tools/lumped/constitutive.py`
and `tools/lumped/vessels.py`, not the generated files).

### Empty lumped vessels (templates)

`haemodynamics/lumped_vessel` holds one empty supermodule per layout: `vp_empty`, `pv_empty`,
`vv_empty`, `pp_empty`, `vp_noI_empty`, `pv_noI_empty`, `pp_noI_empty`. An empty supermodule has
`"template": true`, and each submodule is a **slot**: a `name`, the `module_type` that fits it and the
versions that fit its place (`choices`; for a resistance, the ones of the right form), with no
`module_subtype`:

```json
{"name": "R", "module_type": "resistance", "choices": ["pv_controlled", "pv_geometric", "pv_linear", ...],
 "inp_instances": ["C"], "out_instances": ["I"]}
```

A template is dragged into a model, and a version is chosen for each slot; the result is a lumped
vessel like the `_lumped` versions. libcuflynx refuses to generate a template, naming the slots that
still need a version. Only `supermodule_structure_test` runs on a template (its slots name library
module_types and versions); every other test is skipped (`skip` in its tests.yaml).

### What PhLynx needs for this

PhLynx builds models from CellML module configs only, so it has no supermodules yet (each
`_lumped` version's `phlynx_export_test` is an expected failure). To offer lumped vessels it needs to:

1. **Read supermodule configs from the manifests.** A `"module_format": "supermodule"` entry has no
   CellML of its own; its `submodules` name library versions (`module_type`, `module_subtype`) whose
   CellML, units and configs are in the same manifest (`haemodynamics/` is in the curated set).
2. **Show a supermodule as one node** with the ports its `routes` give: a port type under
   `routes.inputs` is an inlet of the node, under `routes.outputs` an outlet, each with the port's
   variables taken from the routed submodule's config.
3. **Templates: let the user choose each slot.** For a `"template": true` entry, offer each slot's
   `choices` (versions of its `module_type`), and store the choice as the slot's `module_subtype`.
   The filled entry (`template` dropped, each slot with its `module_subtype` and `"instance":
   "default"`) is a supermodule config like the `_lumped` ones, saved with the model (its
   `module_library_dirs` or external modules) under a new `module_subtype`. libcuflynx has no way
   yet to choose the slots inline in a vessel-array record.
4. **Edit parameters per vessel.** Show `shared_parameters` once for the vessel, and the others per
   submodule as `<var>_<submodule>`.
5. **Export the model as libcuflynx reads it:** one record per supermodule
   (`module_type`, `module_subtype`, `instance`, `inp_instances`, `out_instances`), not the expanded
   submodules; libcuflynx expands it and routes the connections. The parameters file names a
   submodule parameter `<var>_<vessel>_<submodule>`, or a shared one `<var>_<vessel>`.
6. **Honour `"available_to_phlynx": false`** for anything that reads configs outside the manifests.

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
| `supermodule_version_equivalence_test` | a supermodule version with `supermodule.equivalent_version`: it reproduces that monolithic version |

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
  of one or more supermodule versions (e.g. `SN_membrane_soma/nn`, which can't be validated alone) is
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
| `outputs` | verification_config | variables to log and check (`var`, or `vessel/var` for a harness neighbour); default: every `variable` of the config (a supermodule: every state, `mod_<submodule>/<var>`) |
| `run_parameters` | verification_config | `{parameter: value}`: the operating point of the run and invariant tests |
| `rest_check` | verification_config | `{sim_time, window, voltage, parameters, threshold, dt}`: a long run (e.g. 0 pA injected) that run_test reports, not a pass/fail gate: the spikes (upward crossings of `threshold` mV, default 0) of `voltage` and their rate in the last `window` s, in the message, metrics.rest_check and plots/rest.png |
| `harness` | verification_config | the test network: `vessel_array` rows `[name, module_subtype, module_type, inp, out, instance]` (the version under test is `mod`) and `parameters` rows `[name, units, value, source]` for the neighbours |
| `invariants` | verification_config | numpy expressions (`expr`, with `description` and `applies`) that must hold |
| `bc_sweep` | verification_config | the parameter sweep: `sweep` (`all`, a component's default: every boundary condition and constant; `bcs`: the boundary conditions; `globals_and_bcs`, a supermodule's default: its globals and open boundary conditions), `factors`, `points`, `ranges` (`[min, max]` or `{min, max, points, scale}` or `{values}`), `bounds` (`{param: [min, max]}`, `null` for no limit: the valid range; values outside it are skipped), `constraints` (numpy expressions of the parameters by variable name, e.g. `B_Ca_init < B_Ca_total`; a point that breaks one is skipped), `exclude`, `extra_parameters`, `plot_output`, `rationale` |
| `timestep` | verification_config | the convergence test: `scheme`, `dts`, `t_end`, `min_order`, `tol`, `cvode_tol`, `roundoff`, `wrapped` |
| `stability` | verification_config | the solver matrix: `supported` (must work), `cvode`, `solve_ivp`, `fixed_step`, `max_step_start`, `min_step`, `time_budget`, `t_end`, `tol` |
| `parameter_ranges` | verification_config | published parameter intervals, `{parameter: [lo, hi]}` (reported as validated values; usually inside a baseline block) |
| `validation` | verification_config | per instance, `validation.<instance>.baseline` / `.calibrate`: `status`, `source`, data references relative to the version directory (`data: instances/<i>/<file>.csv`, `obs_data`, `params_for_id`), variables, targets, metric, threshold, optimiser settings |
| `supermodule` | verification_config | supermodule versions: `globals` (parameters that name no submodule), `equivalent` (system models it must reproduce: `model`, `reproduces`, `instance`, `tol`, `solver_info`, `output_map`, `ignore`) and `equivalent_version` (the monolithic version it reproduces: `version`, `monolithic_parameters`, `output_map` from monolithic outputs to `<submodule>/<var>`, `tol`, `solver_info`; both run alone, or in the monolithic version's `harness` network) |
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
| `<module_type>_<version>_<instance>_vessel_array.json`, `_model_parameters.csv` | The test network and the parameters file the model was generated from. |

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
- where a module_type goes (placement);
- how versions are named (the port letters first);
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
- **A new module_type**: put it in the most specific place covering all its uses: inside the
  module_type it is only used within (a nested module_type), otherwise in the most specific category
  (create the category if none fits and no module_type of that meaning exists, with a name no
  module_type has). Then run `make manifests`.

`tools/restructure_modules.py` moved the old per-module layout into this one, driven by
`tools/restructure_map.yaml` (made by `tools/restructure_map.py`). The later move to nested
module_types (the category `cell/neurons` and `cardiac` replaced by the module_types `cell/neuron` and
`heart`, the alternative hearts made versions of `heart`) is not in that map.
