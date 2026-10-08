#!/bin/bash
# The CI setup cache key (.github/actions/setup): runner image, Python, the setup files, and the
# commit libcuflynx's git requirement resolves to now, so a moved branch or a new image is a new key.
set -euo pipefail
sha=nosha
line=$(grep -E '^libcuflynx @ git\+' requirements.txt || true)
if [ -n "$line" ]; then
  url=${line#*git+}; repo=${url%@*}; ref=${url##*@}
  sha=$(git ls-remote "$repo" "$ref" | head -1 | cut -c1-12)
  sha=${sha:-${ref:0:12}}     # the ref is itself a commit
fi
files=$(cat requirements.txt pyproject.toml ci/apt-packages.txt ci/setup_key.sh .github/actions/setup/action.yml | sha256sum | cut -c1-16)
py=$(python -c 'import sys; print("%d.%d.%d" % sys.version_info[:3])')
echo "setup-${ImageOS:-local}-${ImageVersion:-0}-py$py-$files-$sha"
