#!/bin/bash
# CI's setup-check job: a fresh install (apt + pip, no cache) compared with the cached setup
# (~/.cam-setup) for the same key, printing any difference.
# - Python packages: every one is pinned (ci/constraints.txt; libcuflynx's commit is in the key), so
#   a difference fails: a package outside the pins (refresh ci/constraints.txt, which also changes
#   the key and so rebuilds the cache).
# - apt packages: Ubuntu's updates (e.g. security fixes) to the system libraries change them without
#   changing the key, and the cached .deb files never downgrade what the runner image has, so a
#   difference is a warning; a new runner image rebuilds the cache.
set -uo pipefail
CACHE=$HOME/.cam-setup
if [ ! -x "$CACHE/venv/bin/python" ]; then
  echo "::notice::no cached setup for this key yet (the first job to finish saves it); nothing to compare"
  exit 0
fi
APT="-o Acquire::Retries=3 -o Acquire::http::Timeout=20 -o Acquire::https::Timeout=20"
status=0

# apt: the .deb files a fresh install downloads vs the cached ones (name_version_arch)
mkdir -p "$HOME/fresh-debs/partial"
timeout 120 sudo apt-get $APT update -qq
timeout 300 sudo apt-get $APT install -y -qq --download-only -o Dir::Cache::archives="$HOME/fresh-debs" $(cat ci/apt-packages.txt)
diff <(ls "$CACHE/debs" 2>/dev/null | grep '\.deb$' | sort) <(ls "$HOME/fresh-debs" | grep '\.deb$' | sort) > apt.diff \
  && echo "apt: same packages and versions ($(grep -c . <(ls "$CACHE/debs" | grep '\.deb$')) .deb files)" \
  || { echo "::warning::apt packages differ from the cached setup (< cached, > fresh; Ubuntu updates, rebuilt with the next runner image):"; cat apt.diff; }

# pip: a fresh virtualenv vs the cached one (the repo's own editable install left out)
python -m venv "$HOME/fresh-venv"
"$HOME/fresh-venv/bin/python" -m pip install -q --upgrade pip
"$HOME/fresh-venv/bin/pip" install -q -r requirements.txt -c ci/constraints.txt
"$CACHE/venv/bin/pip" freeze --exclude-editable | sort > cached.txt
"$HOME/fresh-venv/bin/pip" freeze --exclude-editable | sort > fresh.txt
diff cached.txt fresh.txt > pip.diff \
  && echo "pip: same $(wc -l < fresh.txt) packages and versions" \
  || { echo "::error::Python packages differ from the cached setup (< cached, > fresh): pin them in ci/constraints.txt"; cat pip.diff; status=1; }
exit $status
