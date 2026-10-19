"""Supermodule assembly for module hierarchies: a cell type extends a base by adding submodules and links, and every
level is written as a flat supermodule (libcuflynx has no nested supermodules).

Parameter rows of a supermodule are generated from its submodules' own default instances: every constant of a submodule
that no port connection supplies becomes <variable>_<submodule> (global constants keep their bare name), with the value,
data_reference and sourced flag of the submodule's CSV. A cell then overrides the rows its calibration changes.
"""
import csv
import glob
import json
import os

from module_writer import REPO


def module_dir(module_type, version):
    hits = glob.glob(os.path.join(REPO, 'modules', '**', module_type, 'versions', version), recursive=True)
    assert len(hits) == 1, f'{module_type}/{version}: {hits}'
    return hits[0]


def module_config(module_type, version):
    d = module_dir(module_type, version)
    return json.load(open(os.path.join(d, f'{module_type}_{version}_modules_config.json')))[0]


def module_rows(module_type, version, instance='default'):
    d = module_dir(module_type, version)
    return list(csv.DictReader(open(os.path.join(d, 'parameterisations', instance, f'{instance}_parameters.csv'))))


def _ports(cfg, kind):
    return cfg.get(kind) or []


def connected_variables(submodules):
    """{submodule name: set of its variables supplied or taken by a port connection} for a submodule list
    (sub dicts with name, module_type, module_subtype, out_instances)."""
    cfgs = {s['name']: module_config(s['module_type'], s['module_subtype']) for s in submodules}
    out = {s['name']: set() for s in submodules}
    for s in submodules:
        a = cfgs[s['name']]
        a_ports = _ports(a, 'exit_ports') + _ports(a, 'general_ports')
        for tname in s['out_instances']:
            b = cfgs[tname]
            b_ports = _ports(b, 'entrance_ports') + _ports(b, 'general_ports')
            b_types = {p['port_type'] for p in b_ports}
            a_types = {p['port_type'] for p in a_ports}
            for p in a_ports:
                if p['port_type'] in b_types:
                    out[s['name']].update(p['variables'])
            for p in b_ports:
                if p['port_type'] in a_types:
                    out[tname].update(p['variables'])
    return out


def supermodule_rows(submodules, instances=None):
    """Parameter rows [name, units, value, data_reference, sourced] of a flat supermodule from its submodules' CSVs.
    instances: {submodule name: instance} where it is not the submodule's own 'instance' entry."""
    conn = connected_variables(submodules)
    rows, seen = [], {}
    for s in submodules:
        cfg = module_config(s['module_type'], s['module_subtype'])
        roles = {v[0]: v[3] for v in cfg['variables_and_units']}
        inst = (instances or {}).get(s['name'], s.get('parameterisation', s.get('instance', 'default')))
        for r in module_rows(s['module_type'], s['module_subtype'], inst):
            n = r['variable_name']
            role = roles.get(n)
            if role == 'global_constant':
                if n not in seen:
                    seen[n] = len(rows)
                    rows.append([n, r['units'], float(r['value']), r['data_reference'], r['sourced']])
                continue
            if n in conn[s['name']] and role in ('boundary_condition', 'variable'):
                continue
            rows.append([f'{n}_{s["name"]}', r['units'], float(r['value']), r['data_reference'], r['sourced']])
    return rows
