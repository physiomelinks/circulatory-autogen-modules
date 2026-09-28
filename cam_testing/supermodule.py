"""
Supermodules: modules/supermodules/<name>/, a vessel array of library modules with named
interface ports, used by a host model as one vessel.

  <name>_vessel_array.csv   the internal vessels; '@<port>' in inp_vessels/out_vessels marks
                            where a host vessel connects
  <name>_parameters.csv     default parameters of the internal vessels (the host's win)
  <name>_supermodule.yaml   the ports ({vessels, side: inp|out}) and test harnesses

A host vessel array uses a supermodule through one row

    heart,supermodule,supermodule:heart,systemic_venous_in:venous_svc pulmonary_venous_in:pvn,pulmonary_arterial_out:par ...

and names the instance (here 'heart') in its own vessels' inp/out lists. expand() replaces the
row by the internal vessels, the '@<port>' markers by the host vessels, and the instance name in
each host vessel's lists by the internal vessels of the port it connects to. Unconnected
optional ports are dropped. Internal vessel names are used as they are, so one host can hold one
instance of a supermodule. libcuflynx can't load supermodules yet; this flattening is the
stand-in until it can.
"""
import csv
import glob
import os
from dataclasses import dataclass

import yaml

from cam_testing.library import MODULES_DIR

SUPERMODULE_DIR = os.path.join(MODULES_DIR, 'supermodules')
VESSEL_FIELDS = ['name', 'BC_type', 'vessel_type', 'inp_vessels', 'out_vessels']


@dataclass
class Supermodule:
    name: str
    dir: str
    spec: dict

    def path(self, suffix):
        return os.path.join(self.dir, f'{self.name}_{suffix}')


def supermodules():
    out = []
    for spec_path in sorted(glob.glob(os.path.join(SUPERMODULE_DIR, '*', '*_supermodule.yaml'))):
        with open(spec_path) as f:
            spec = yaml.safe_load(f) or {}
        out.append(Supermodule(spec.get('supermodule') or os.path.basename(os.path.dirname(spec_path)),
                               os.path.dirname(spec_path), spec))
    return out


def load_supermodule(name):
    for s in supermodules():
        if s.name == name:
            return s
    raise KeyError(f'no supermodule {name!r} in {SUPERMODULE_DIR}')


def read_vessels(path):
    with open(path) as f:
        rows = [{k: (v or '').strip() for k, v in r.items()} for r in csv.DictReader(f)]
    return [r for r in rows if r['name'] and not r['name'].startswith('#')]


def read_parameters(path):
    with open(path) as f:
        reader = csv.reader(f)
        header = next(reader)
        return header, [r for r in reader if r and r[0].strip()]


def _split(s):
    return s.split()


def _port_map(tokens, instance):
    out = {}
    for tok in tokens:
        if ':' not in tok:
            raise ValueError(f'supermodule instance {instance}: {tok!r} should be <port>:<host vessel>')
        port, vessel = tok.split(':', 1)
        out[port] = vessel
    return out


def expand(rows):
    '''Host vessel rows (dicts) -> rows with every supermodule instance flattened.'''
    instances = [r for r in rows if r['vessel_type'].startswith('supermodule:')]
    if not instances:
        return [dict(r) for r in rows]
    out = [dict(r) for r in rows]
    for inst in instances:
        sm = load_supermodule(inst['vessel_type'].split(':', 1)[1])
        ports = sm.spec.get('ports') or {}
        connect = {**_port_map(_split(inst['inp_vessels']), inst['name']),
                   **_port_map(_split(inst['out_vessels']), inst['name'])}
        unknown = set(connect) - set(ports)
        if unknown:
            raise ValueError(f'supermodule {sm.name} has no port(s) {sorted(unknown)} (ports: {sorted(ports)})')
        missing = [p for p, d in ports.items() if p not in connect and not d.get('optional')]
        if missing:
            raise ValueError(f'supermodule instance {inst["name"]}: port(s) {missing} not connected')

        internal = read_vessels(sm.path('vessel_array.csv'))
        clash = {r['name'] for r in internal} & {r['name'] for r in out if not r['vessel_type'].startswith('supermodule:')}
        if clash:
            raise ValueError(f'supermodule {sm.name}: internal vessel(s) {sorted(clash)} clash with host vessels')
        for r in internal:
            for col in ('inp_vessels', 'out_vessels'):
                names = []
                for n in _split(r[col]):
                    if n.startswith('@'):
                        if n[1:] not in ports:
                            raise ValueError(f'supermodule {sm.name}: {n} is not one of its ports')
                        if n[1:] in connect:
                            names.append(connect[n[1:]])
                    else:
                        names.append(n)
                r[col] = ' '.join(names)

        # host vessels: the instance name -> the internal vessels of the port(s) they connect to
        for r in out:
            if r['vessel_type'].startswith('supermodule:'):
                continue
            for col, side in (('out_vessels', 'inp'), ('inp_vessels', 'out')):
                names = []
                for n in _split(r[col]):
                    if n != inst['name']:
                        names.append(n)
                        continue
                    via = [p for p, v in connect.items() if v == r['name'] and ports[p]['side'] == side]
                    if not via:
                        raise ValueError(f'{r["name"]} lists {inst["name"]} in {col} but connects to none of its '
                                         f'{"entrance" if side == "inp" else "exit"} ports')
                    for p in via:
                        names.extend(v for v in ports[p]['vessels'] if v not in names)
                r[col] = ' '.join(names)

        i = next(k for k, r in enumerate(out) if r['name'] == inst['name'] and r['vessel_type'] == inst['vessel_type'])
        out[i:i + 1] = internal
    return out


def expand_parameters(host_rows, host_params, header=None):
    '''The host's parameters plus each instanced supermodule's defaults the host doesn't set.'''
    header = header or ['variable_name', 'units', 'value', 'data_reference']
    params = [list(p) for p in host_params]
    have = {p[0].strip() for p in params}
    for inst in (r for r in host_rows if r['vessel_type'].startswith('supermodule:')):
        sm = load_supermodule(inst['vessel_type'].split(':', 1)[1])
        _, defaults = read_parameters(sm.path('parameters.csv'))
        for p in defaults:
            if p[0].strip() not in have:
                params.append(p)
                have.add(p[0].strip())
    return header, params


def write_model(host_vessel_array, host_parameters, dest_dir, name):
    '''Writes <dest_dir>/<name>_vessel_array.csv and _parameters.csv with supermodules expanded.'''
    rows = read_vessels(host_vessel_array)
    header, params = read_parameters(host_parameters)
    os.makedirs(dest_dir, exist_ok=True)
    with open(os.path.join(dest_dir, f'{name}_vessel_array.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=VESSEL_FIELDS, extrasaction='ignore')
        w.writeheader()
        w.writerows(expand(rows))
    header, params = expand_parameters(rows, params, header)
    with open(os.path.join(dest_dir, f'{name}_parameters.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(params)
    return dest_dir
