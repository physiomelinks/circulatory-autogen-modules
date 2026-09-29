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

.PHONY: setup structure test test-all systems pipeline pipeline-setup report site serve clean manifests risk

setup:
	python3 -m venv venv
ifneq ($(LIBCUFLYNX),)
	venv/bin/pip install -e $(LIBCUFLYNX) -e .
else
	venv/bin/pip install -r requirements.txt
endif

structure:
	$(PYTHON) -m pytest tests/test_structure.py

# V&V tests, skipping slow calibration. MODULE=name runs one module (reviewed or not).
test:
	$(PYTHON) -m pytest tests/test_modules.py -m "not slow" $(MODULE_ARGS) $(PYTEST_ARGS)

test-all:
	$(PYTHON) -m pytest tests/test_modules.py $(MODULE_ARGS) $(PYTEST_ARGS)

# the system models (modules/system) and supermodules
systems:
	$(PYTHON) -m pytest tests/test_systems.py tests/test_supermodules.py $(PYTEST_ARGS)

# PhLynx build & export -> CUFLynx import & simulate -> compare with libcuflynx, per component
pipeline-setup:
	cd tools/phlynx_bridge && npm ci --no-audit --no-fund

pipeline:
	PHLYNX_DIR=$(PHLYNX_DIR) CUFLYNX_BIN=$(CUFLYNX_BIN) $(PYTHON) -m pytest tests/test_phlynx.py $(MODULE_ARGS) $(PYTEST_ARGS)

report:
	$(PYTHON) -m cam_testing.report $(MODULE_ARGS)

# site/ exactly as GitHub Pages serves it
site:
	$(PYTHON) -m cam_testing.report --site $(MODULE_ARGS)

serve: site
	$(PYTHON) -m http.server -d site 8000

manifests:
	$(PYTHON) tools/build_manifests.py

clean:
	rm -rf site modules/*/plots modules/*/results modules/*/*.html
