#!/bin/bash
# CI's setup-check job: a fresh install (apt + pip, no cache) compared with the cached setup
# (~/.cam-setup, ~/.cam-setup-apt) for the same keys, printing any difference.
# - Python packages: every one is pinned (ci/constraints.txt; libcuflynx's commit is in the key), so
#   a difference fails: a package outside the pins (refresh ci/constraints.txt, which also changes
#   the key and so rebuilds the cache).
# - apt packages: Ubuntu's updates (e.g. security fixes) to the system libraries change them without
#   changing the key (ci/apt-packages.txt only), so a difference in a top-level package's version is a
#   warning.
set -uo pipefail
export DEBIAN_FRONTEND=noninteractive   # no debconf prompt can block (sudo -E keeps it)
CACHE=$HOME/.cam-setup
APT_CACHE=$HOME/.cam-setup-apt
if [ ! -x "$CACHE/venv/bin/python" ]; then
  echo "::notice::no cached setup for this key yet (the first job to finish saves it); nothing to compare"
  exit 0
fi
APT="-o Acquire::Retries=3 -o Acquire::http::Timeout=20 -o Acquire::https::Timeout=20"
status=0

# apt: the version of each package of ci/apt-packages.txt in the cached .deb files vs the candidate
# a fresh apt-get update offers (the dependency closure is not compared)
timeout 120 sudo apt-get $APT update -qq
for p in $(cat ci/apt-packages.txt); do
  deb=$(ls "$APT_CACHE"/debs/${p}_*.deb 2>/dev/null | head -1)
  cached=${deb:+$(dpkg-deb -f "$deb" Version)}
  fresh=$(apt-cache policy "$p" | awk '/Candidate:/ {print $2}')
  echo "$p ${cached:-none} $fresh"
done > apt.versions
awk '$2!=$3' apt.versions > apt.diff
[ -s apt.diff ] \
  && { echo "::warning::apt package versions differ from the cached setup (package, cached, fresh; Ubuntu updates, rebuilt when ci/apt-packages.txt changes):"; cat apt.diff; } \
  || echo "apt: same versions for $(wc -l < apt.versions) packages"

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
