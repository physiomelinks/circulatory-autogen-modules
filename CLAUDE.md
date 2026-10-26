# circulatory-autogen-modules

## Pull requests

A PR is a **module PR**, a **migration PR** or another PR (README "Pull requests"):

- A module PR's title starts with `Module:` and it changes only module content: `modules/` (not
  `modules/README.md` or `modules/directory_schema.json`), `system_models/`, `manifests/` and `reviews/`.
- Anything else (`cam_testing`, `tests/`, `tools/`, CI, docs, the schema, setup) goes in a separate PR
  without the prefix, which changes no module content; merge it first, then merge `main` into the module PR.
  A change to the library's modules made by an infrastructure change (e.g. a rename across modules or
  regenerated manifests) goes in a `Module:` PR of its own.
- A migration PR's title starts with `Migration:` and it may change both: a format conversion that changes
  module content together with the code that reads it (e.g. converting the supermodule configs and `cam_testing`).
- During a module review walkthrough, put infrastructure fixes on their own branch as you go, not on the
  review branch. Check with `python tools/pr_scope.py "<title>" origin/main` before opening the PR; the
  `PR scope` workflow runs the same check.

## Data sources

- Schmidt-Nielsen, *Animal Physiology: Adaptation and Environment*, fifth edition (1997; bib key
  `schmidtnielsen1997animal`), a local copy at
  `/home/farg967/Documents/books/Physiology/Animal physiology adaptation and environment Fifth Edition.pdf`:
  a good source of parameter values and validation data for modules. Cite it in the version's `references.bib`,
  and keep `source_figures/` screenshots of the figures and tables used.
