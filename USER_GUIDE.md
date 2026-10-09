# Using the module library

circulatory-autogen-modules is a library of tested model components (**modules**) for building
physiological models with [libcuflynx](https://github.com/physiomelinks/circulatory_autogen)
(circulatory_autogen), [PhLynx](https://github.com/physiomelinks/PhLynx) and
[CUFLynx](https://github.com/physiomelinks/CUFLynx): vessels, heart chambers and valves, ion
channels, neurons, control loops, gas exchange, boundary conditions and more, plus whole
**system models** built from them.

This guide is for people who want to **use** the library. To add or change a module, see
[`modules/README.md`](modules/README.md) ("Adding things").

## The words used here

| Term | Meaning | Where |
|---|---|---|
| **module_type** | A mechanism, e.g. `i_CaL` (an L-type Ca channel) or `arterial` | `modules/<category>/<module_type>/` |
| **version** | One set of equations for it. Different maths makes a different version, e.g. `Argus2026_v01` | `versions/<version>/` |
| **parameterisation** | The same equations with a particular set of values: the version's `default`, or a set fitted to a dataset or experiment, with its data, calibration and source figures | `versions/<version>/parameterisations/<name>/` |
| **instance** | One use of a version in a model: a record of a module array (or a node in PhLynx). It names the parameterisation it uses, and several instances can share one | `*_module_array.json` |
| **system model** | A whole model built from library modules | `system_models/<category>/<model>/` |

(Parameterisations were called *instances* before October 2026. The old names, `instances/`,
`default_instance` and a record's `"instance"`, are still read.)

## Browse the library

The reports are published at **https://physiomelinks.github.io/circulatory-autogen-modules/**.

- The **index** groups module_types by category, with each one's versions, test status and
  PhLynx/CUFLynx compatibility.
- A **module_type page** compares its versions.
- A **version page** shows:
  - the equations and ports;
  - every parameter with its value, units, source (`data_reference`) and whether that source
    was checked (`sourced`);
  - the test results, with plots;
  - the required citations;
  - each parameterisation, with its data, validation plots, source figures and a downloadable
    CUFLynx archive (`.omex`).

## Can I trust a version?

Each version page shows:

- **Reviewed**: whether the version has been through review: its equations, parameter sources,
  tests and known issues checked with its owner. **A version not reviewed yet fails this
  check.** It is still public and usable, but treat its values and behaviour as provisional.
- **The test columns**:
  - *verification*: it runs, conserves what it should, converges with the time step, and stays
    stable across solvers and parameter ranges;
  - *validation*: each parameterisation against its data;
  - *the PhLynx → CUFLynx pipeline*: the version exported by PhLynx and run in CUFLynx
    reproduces libcuflynx.

  A cell shows Pass, Fail, Not applicable, or Not run. Not run is used, for example, when
  continuous integration (CI) only runs the quick tests for an unreviewed version. Known
  issues and expected failures are listed with their reasons.
- **Parameter sources**: each value cites its source; values not taken from a source say so
  (`unreferenced_placeholder`, `definitional`), and large departures from the literature are in
  **bold**.

The tests run on every change (only the versions affected by it) and on all versions weekly.

## Use a module or system model

### In CUFLynx
Download a parameterisation's `.omex` archive from its version page (or build them all with
`make omex`) and open it in CUFLynx. The archive holds the flattened model, the
parameterisation's parameters and, when it has them, its observation data and parameters to
identify, ready to simulate or calibrate.

### In PhLynx
PhLynx loads the versions listed in `manifests/index.json`. Place them as nodes (instances) and
export the model; the exported archive runs in CUFLynx.

### With libcuflynx
Point libcuflynx at the library in your `user_inputs.yaml`:

```yaml
module_library_dirs: [/path/to/circulatory-autogen-modules/modules]
use_builtin_modules: false
```

Then name library versions in your model's module array. Each record is an instance; give it
the parameterisation to use (default: the version's `default`):

```json
[
  {"name": "aorta", "module_type": "arterial_simple", "module_subtype": "vp", "parameterisation": "default",
   "inp_instances": ["heart"], "out_instances": ["artery_2"]},
  {"name": "artery_2", "module_type": "arterial_simple", "module_subtype": "vp", "parameterisation": "default",
   "inp_instances": ["aorta"], "out_instances": ["venous"]}
]
```

An instance's parameters come from its parameterisation and are renamed per instance (the
parameterisation's `C` becomes `C_aorta` and `C_artery_2`). A row in your model's own parameters
file overrides them, so two instances of one parameterisation can still differ.

The system models in `system_models/` are complete examples: each has its module array,
parameters file and the checks it passes. Your libcuflynx needs the `module_library_dirs`
support (on circulatory_autogen master); `requirements.txt` pins the version this library is
tested with.

## Cite what you use

- Every version lists its **required citations**: the papers its model is built on (config key
  `"required_citations"`; BibTeX in the version's `*_references.bib`). Cite them when you use the
  version. The version page shows them with links.
- Individual parameter values cite their own sources in the parameters file (`data_reference`).
  Cite those you rely on.
- Some versions also list `required_citations_uncertain`, a best guess still being checked, and
  owner-original models cite placeholders (`argus2026...`) until their papers are published.

## Licence

Each module version states its licence in its config (`"licence"`); the library's modules are
**CC0-1.0** (public domain dedication). The test and report code (`cam_testing/`, `tools/`) is
under the repository's Apache-2.0 licence (`LICENSE`).

## Test your own modules the same way

`cam_testing` (this repository's test package) runs the same checks on any repository laid out
like this one, or on your own module library next to this one. See README.md, "Using cam_testing
from another repository".

## Contribute

Adding a module, a version or a parameterisation, and the review process, are described in
[`modules/README.md`](modules/README.md). Questions and problems: open an issue on this repository.
