# Local entry points; CI runs exactly these targets.
PYTHON ?= venv/bin/python
MODULE ?=
MODULE_ARGS = $(if $(MODULE),--module $(MODULE),)
# a local libcuflynx checkout to install instead of the pinned git ref, e.g. LIBCUFLYNX=../circulatory_autogen
LIBCUFLYNX ?=
# extra pytest arguments, e.g. PYTEST_ARGS=--quick-unreviewed (CI)
PYTEST_ARGS ?=
# PhLynx -> CUFLynx pipeline: a PhLynx checkout (with node_modules) and a released CUFLynx binary
PHLYNX_DIR ?= ../phlynx
CUFLYNX_BIN ?= $(HOME)/software/CUFLynx

.PHONY: setup structure test test-all systems pipeline pipeline-setup omex omex-test report site serve clean manifests risk

setup:
	python3 -m venv venv
ifneq ($(LIBCUFLYNX),)
	venv/bin/pip install -e $(LIBCUFLYNX) -e .
else
	venv/bin/pip install -r requirements.txt
endif

structure:
	$(PYTHON) -m pytest tests/test_structure.py

# V&V tests, skipping slow calibration. MODULE=<module_type or path under modules/> runs those versions
# (reviewed or not), e.g. MODULE=Lotka_Volterra, MODULE=heart (with cardiac_clock, chamber, valve) or MODULE=cell.
test:
	$(PYTHON) -m pytest tests/test_modules.py -m "not slow" $(MODULE_ARGS) $(PYTEST_ARGS)

test-all:
	$(PYTHON) -m pytest tests/test_modules.py $(MODULE_ARGS) $(PYTEST_ARGS)

# the system models (system_models/) and the supermodule versions' structure and "reproduces" tests
systems:
	$(PYTHON) -m pytest tests/test_systems.py $(PYTEST_ARGS)
	$(PYTHON) -m pytest tests/test_modules.py -k supermodule --include-unreviewed $(PYTEST_ARGS)

# PhLynx build & export -> CUFLynx import & simulate -> compare with libcuflynx, per version
pipeline-setup:
	cd cam_testing/bridges/phlynx && npm ci --no-audit --no-fund

pipeline:
	PHLYNX_DIR=$(PHLYNX_DIR) CUFLYNX_BIN=$(CUFLYNX_BIN) $(PYTHON) -m pytest tests/test_phlynx.py $(MODULE_ARGS) $(PYTEST_ARGS)

# a COMBINE archive per instance for CUFLynx, generated from the library's files (not committed)
omex:
	$(PYTHON) tools/build_instance_omex.py $(MODULE_ARGS)

# each instance's archive loads and runs in a released CUFLynx and reproduces libcuflynx
omex-test:
	CUFLYNX_BIN=$(CUFLYNX_BIN) $(PYTHON) -m pytest tests/test_instance_omex.py $(MODULE_ARGS) $(PYTEST_ARGS)

report:
	$(PYTHON) -m cam_testing.report $(MODULE_ARGS)

# site/ exactly as GitHub Pages serves it
site:
	$(PYTHON) -m cam_testing.report --site $(MODULE_ARGS)

serve: site
	$(PYTHON) -m http.server -d site 8000

manifests:
	$(PYTHON) -m cam_testing.manifests

clean:
	rm -rf site
	find modules system_models -type d \( -name plots -o -name results \) -prune -exec rm -rf {} +
	find modules -name '*.html' -delete

# joint failure-risk analysis (on request): MODULE=<module_type or category>; RISK_ARGS=--reviewed-only in CI
RISK_ARGS ?=
risk:
	$(PYTHON) -m cam_testing.risk $(MODULE_ARGS) $(RISK_ARGS)
