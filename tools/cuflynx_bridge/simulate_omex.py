"""
CUFLynx half of the module pipeline test: launches a released CUFLynx binary, imports PhLynx
.omex exports through its HTTP API (as PhLynx's hand-off and the upload box do), and simulates
them. Standard library only, so it runs in any Python.

    CUFLYNX_BIN=~/software/CUFLynx python tools/cuflynx_bridge/simulate_omex.py job.json

job.json: {"models": [{"id", "omex", "sim_time", "extra_outputs", "out"}], "solver_info": {...}}
For each model: POST /api/omex/upload, then POST /api/simulate for the model's states plus
extra_outputs ('component/variable', skipped when CUFLynx can't resolve them), writing <out> with the time, outputs and status. Prints a JSON summary on stdout.

The binary runs as `CUFLynx --port <free> --browser` (serve only) with BROWSER=true, so no
browser window opens, and a throwaway CUFLYNX_CONFIG_DIR, so the user's saved settings (e.g.
a CA dir) neither leak in nor get overwritten: the bundled engine is what simulates.
Exit 3 when the binary is missing or never becomes healthy.
"""
import json
import os
import socket
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request

DEFAULT_BIN = os.path.expanduser('~/software/CUFLynx')


def _req(method, url, data=None, timeout=600):
    body = json.dumps(data).encode() if data is not None else None
    req = urllib.request.Request(url, data=body, method=method,
                                 headers={'Content-Type': 'application/json'} if body else {})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, json.loads(r.read().decode() or 'null')
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode(errors='replace')


def _upload(url, path, timeout=600):
    boundary = '----cam-cuflynx-bridge'
    with open(path, 'rb') as f:
        payload = f.read()
    body = (f'--{boundary}\r\nContent-Disposition: form-data; name="file"; filename="{os.path.basename(path)}"\r\n'
            'Content-Type: application/zip\r\n\r\n').encode() + payload + f'\r\n--{boundary}--\r\n'.encode()
    req = urllib.request.Request(url, data=body, method='POST',
                                 headers={'Content-Type': f'multipart/form-data; boundary={boundary}'})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode(errors='replace')


def _free_port():
    with socket.socket() as s:
        s.bind(('127.0.0.1', 0))
        return s.getsockname()[1]


def launch(binary):
    port = _free_port()
    base = f'http://127.0.0.1:{port}'
    env = dict(os.environ, BROWSER='true', CUFLYNX_CONFIG_DIR=tempfile.mkdtemp(prefix='cam_cuflynx_config_'),
               CUFLYNX_OUTPUT_DIR=tempfile.mkdtemp(prefix='cam_cuflynx_out_'))
    env.pop('CIRCULATORY_AUTOGEN_SRC', None)
    log = open(os.path.join(env['CUFLYNX_OUTPUT_DIR'], 'cuflynx.log'), 'w')
    proc = subprocess.Popen([binary, '--port', str(port), '--browser'], stdout=log, stderr=subprocess.STDOUT,
                            env=env, cwd=env['CUFLYNX_OUTPUT_DIR'])
    deadline = time.time() + 180
    while time.time() < deadline:
        if proc.poll() is not None:
            break
        try:
            status, health = _req('GET', f'{base}/api/health', timeout=2)
            if status == 200:
                return proc, base, health
        except OSError:
            pass
        time.sleep(1)
    proc.kill()
    log.close()
    with open(log.name) as f:
        tail = f.read()[-2000:]
    print(f'BRIDGE_UNSUPPORTED: CUFLynx at {binary} did not become healthy:\n{tail}', file=sys.stderr)
    sys.exit(3)


def run(job, binary):
    proc, base, health = launch(binary)
    summary = []
    try:
        config = {'generated_model_format': 'cellml', 'solver': 'CVODE_myokit',
                  'solver_info': dict(job.get('solver_info') or {}, method='CVODE')}
        config_status, _ = _req('POST', f'{base}/api/config', config)
        for m in job['models']:
            res = {'id': m['id'], 'ok': False, 'stage': 'upload', 'cuflynx_version': health.get('version'),
                   'config_status': config_status}
            status, body = _upload(f'{base}/api/omex/upload', m['omex'])
            if status != 200:
                res['error'] = f'upload {status}: {str(body)[:2000]}'
            else:
                res['model_id'] = body.get('model_id')
                res['import_warnings'] = body.get('warnings') or []
                # what CUFLynx's importer took from the archive
                def _role(x):
                    return {k: x.get(k) for k in ('filename', 'error')} if isinstance(x, dict) else None
                res['imported'] = {'model_filename': body.get('model_filename'),
                                   'obs_data': _role(body.get('obs_data')),
                                   'params_for_id': _role(body.get('params_for_id'))}
                res['stage'] = 'simulate'
                # the states, plus the component's own outputs where CUFLynx can resolve them
                outputs = list(dict.fromkeys((body.get('odes') or []) + list(m.get('extra_outputs') or [])))
                status, out = _req('POST', f'{base}/api/simulate', {
                    'model_id': res['model_id'], 'params': {}, 'sim_time': float(m.get('sim_time', 1.0)),
                    'pre_time': 0.0, 'outputs': outputs or None, 'best_effort_outputs': True})
                if status != 200:
                    res['error'] = f'simulate {status}: {str(out)[:2000]}'
                else:
                    res['time'] = out.get('time') or []
                    res['outputs'] = out.get('outputs') or {}
                    res['n_points'] = len(res['time'])
                    res['ok'] = res['n_points'] > 1 and bool(res['outputs'])
                    res['stage'] = 'done'
            with open(m['out'], 'w') as f:
                json.dump(res, f)
            summary.append({k: res.get(k) for k in ('id', 'ok', 'stage', 'error', 'n_points')})
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=20)
        except subprocess.TimeoutExpired:
            proc.kill()
    return summary


if __name__ == '__main__':
    binary = os.environ.get('CUFLYNX_BIN') or DEFAULT_BIN
    if not os.access(binary, os.X_OK):
        print(f'BRIDGE_UNSUPPORTED: no CUFLynx binary at {binary} (set CUFLYNX_BIN)', file=sys.stderr)
        sys.exit(3)
    with open(sys.argv[1]) as f:
        print(json.dumps(run(json.load(f), binary)))
