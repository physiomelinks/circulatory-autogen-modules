# Local entry points; CI runs exactly these targets.
PYTHON ?= venv/bin/python
MODULE ?=
MODULE_ARGS = $(if $(MODULE),--module $(MODULE),)
# a local libcuflynx checkout to install instead of the pinned git ref, e.g. LIBCUFLYNX=../circulatory_autogen
LIBCUFLYNX ?=
# extra pytest arguments, e.g. PYTEST_ARGS=--quick-unreviewed (CI)
PYTEST_ARGS ?=

.PHONY: setup structure test test-all report site serve clean manifests ranges risk

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
	$(PYTHON) -m pytest tests $(MODULE_ARGS) $(PYTEST_ARGS)

report:
	$(PYTHON) -m cam_testing.report $(MODULE_ARGS)

# site/ exactly as GitHub Pages serves it
site:
	$(PYTHON) -m cam_testing.report --site $(MODULE_ARGS)

serve: site
	$(PYTHON) -m http.server -d site 8000

# record verified/validated parameter ranges from the latest test results into <name>_parameters.csv
# joint failure-risk analysis (Sobol samples, corner plots), on request
SAMPLES ?= 512
risk:
	$(PYTHON) -m cam_testing.risk $(MODULE_ARGS) --samples $(SAMPLES)

ranges:
	$(PYTHON) -m cam_testing.ranges update $(MODULE_ARGS)

manifests:
	$(PYTHON) tools/build_manifests.py

clean:
	rm -rf site modules/*/plots modules/*/results modules/*/*.html
