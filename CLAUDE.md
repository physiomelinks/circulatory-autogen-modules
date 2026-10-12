# circulatory-autogen-modules

## Pull requests

A PR is either a **module PR** or another PR (README "Pull requests"):

- A module PR's title starts with `Module:` and it changes only module content: `modules/` (not
  `modules/README.md` or `modules/directory_schema.json`), `system_models/`, `manifests/` and `reviews/`.
- Anything else (`cam_testing`, `tests/`, `tools/`, CI, docs, the schema, setup) goes in a separate PR
  without the prefix, which changes no module content; merge it first, then merge `main` into the module PR.
  A change to the library's modules made by an infrastructure change (e.g. a rename across modules or
  regenerated manifests) goes in a `Module:` PR of its own.
- During a module review walkthrough, put infrastructure fixes on their own branch as you go, not on the
  review branch. Check with `python tools/pr_scope.py "<title>" origin/main` before opening the PR; the
  `PR scope` workflow runs the same check.
