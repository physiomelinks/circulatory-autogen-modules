"""
Supermodules: modules/supermodules/<name>/, a module made of other modules.

  <name>_modules_config.json   a config entry with module_format "supermodule": its module_type /
                               module_subtype, the submodules (a vessel array of library modules,
                               internal connections only) and default_parameters
  <name>_parameters.csv        default parameters, {var}_{submodule} (local names) or globals
  <name>_supermodule.yaml      description and tests

A model uses one instance, e.g. {"name": "heart", "module_type": "heart", "module_subtype":
"supermodule", "per_submodule_inputs": {"ra": ["venous_svc"]}, "per_submodule_outputs": {"aov":
["aortic_root"]}}. libcuflynx expands it: each submodule becomes <instance>_<submodule>, the host
modules in per_submodule_inputs / per_submodule_outputs are coupled to that submodule, and the
default parameters are renamed to the prefixed names (the model's own parameter values win).
"""
import glob
import json
import os
from dataclasses import dataclass

import yaml

from cam_testing.library import MODULES_DIR

SUPERMODULE_DIR = os.path.join(MODULES_DIR, 'supermodules')


@dataclass
class Supermodule:
    name: str
    dir: str
    entry: dict       # the module_format "supermodule" config entry
    spec: dict

    @property
    def submodules(self):
        return self.entry.get('submodules') or []

    @property
    def submodule_names(self):
        return [s['name'] for s in self.submodules]

    @property
    def defaults_path(self):
        d = self.entry.get('default_parameters')
        return os.path.join(self.dir, d) if d else None


def supermodules():
    out = []
    for cfg in sorted(glob.glob(os.path.join(SUPERMODULE_DIR, '*', '*_modules_config.json'))):
        d = os.path.dirname(cfg)
        name = os.path.basename(d)
        spec_path = os.path.join(d, f'{name}_supermodule.yaml')
        spec = {}
        if os.path.isfile(spec_path):
            with open(spec_path) as f:
                spec = yaml.safe_load(f) or {}
        with open(cfg) as f:
            for entry in json.load(f):
                if entry.get('module_format') == 'supermodule':
                    out.append(Supermodule(entry.get('module_type', name), d, entry, spec))
    return out


def load_supermodule(name):
    for s in supermodules():
        if s.name == name:
            return s
    raise KeyError(f'no supermodule {name!r} in {SUPERMODULE_DIR}')


def prefixed_output_map(names, instance, submodules):
    '''Output names of a model with the submodules flattened ('ra/u') -> the supermodule model's ('heart_ra/u').'''
    subs = set(submodules)
    out = {}
    for n in names:
        vessel, _, var = n.partition('/')
        out[n] = f'{instance}_{vessel}/{var}' if vessel in subs else n
    return out
