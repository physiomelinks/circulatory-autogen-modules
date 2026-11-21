#!/bin/bash
# The CI setup cache keys (.github/actions/setup), one per cache; the setup-cache job computes them
# once and the other jobs of the run take them as inputs.
#   debs: the system libraries (.deb files, their dependency closure and the apt lists): the runner
#         OS and ci/apt-packages.txt only, so a new runner image build or a moved libcuflynx never
#         rebuilds them (apt installs what the image lacks from the cached files).
#   venv: the Python packages: runner OS, Python, the requirement files and the commit libcuflynx's
#         git requirement resolves to now (a moved branch is a new key).
set -euo pipefail
case "${1:-}" in
  debs)
    echo "debs-${ImageOS:-local}-$(sha256sum ci/apt-packages.txt | cut -c1-16)"
    ;;
  venv)
    sha=nosha
    line=$(grep -E '^libcuflynx @ git\+' requirements.txt || true)
    if [ -n "$line" ]; then
      url=${line#*git+}; repo=${url%@*}; ref=${url##*@}
      sha=$(git ls-remote "$repo" "$ref" | head -1 | cut -c1-12)
      sha=${sha:-${ref:0:12}}     # the ref is itself a commit
    fi
    files=$(cat requirements.txt pyproject.toml ci/constraints.txt | sha256sum | cut -c1-16)
    py=$(python -c 'import sys; print("%d.%d.%d" % sys.version_info[:3])')
    echo "venv-${ImageOS:-local}-py$py-$files-$sha"
    ;;
  *)
    echo "usage: ci/setup_key.sh debs|venv" >&2
    exit 2
    ;;
esac
