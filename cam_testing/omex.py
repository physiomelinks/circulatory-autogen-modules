"""
Per-instance COMBINE archives (.omex) for CUFLynx, generated from the library's own files so
nothing is kept in two places:

    modules/<module_type path>/versions/<version>/instances/<instance>/<module_type>_<version>_<instance>.omex

Built by ``python tools/build_instance_omex.py`` (``make omex``); not committed. Each archive holds

  <module_type>_<version>_<instance>.cellml   master: the flattened model of the version's test
                                              network (its harness, or the version alone) with this
                                              instance's parameters, as libcuflynx generates it
  <instance>_obs_data.json                    the instance's calibration data, if it has any
  <instance>_params_for_id.csv                its parameters to identify, if any
  inputs, unchanged: <instance>_parameters.csv, raw data and
    SOURCES.md; the version's _modules.cellml, _modules_config.json (with its required_citations),
    _units.cellml and _verification_config.json; a supermodule's submodules' _modules_config.json
    (their required_citations); and the test network the model was generated from
    (<module_type>_<version>_<instance>_module_array.json / _model_parameters.csv)
  manifest.xml                                COMBINE manifest (master marked)

CUFLynx classifies members by name: the first ``*.cellml`` the manifest marks master is the
model, the first JSON with "obs" in its name is the obs_data, and the first CSV with "param" in
its name is the params_for_id. So the study members are written first, in that order.
"""
import io
import os
import shutil
import tempfile
import zipfile
from xml.sax.saxutils import quoteattr

FORMATS = {
    '.cellml': 'http://identifiers.org/combine.specifications/cellml',
    '.json': 'http://purl.org/NET/mediatypes/application/json',
    '.csv': 'http://purl.org/NET/mediatypes/text/csv',
    '.md': 'http://purl.org/NET/mediatypes/text/markdown',
    '.txt': 'http://purl.org/NET/mediatypes/text/plain',
    '.yaml': 'http://purl.org/NET/mediatypes/application/x-yaml',
}
OMEX_FORMAT = 'http://identifiers.org/combine.specifications/omex'
VERIFICATION_CONFIG_FORMAT = 'http://purl.org/NET/mediatypes/application/x.vnd.cam-verification-config+json'


def omex_path(version, instance):
    return os.path.join(instance.dir, f'{version.vessel_type}_{version.name}_{instance.name}.omex')


def _manifest(entries):
    lines = ['<?xml version="1.0" encoding="UTF-8"?>',
             '<omexManifest xmlns="http://identifiers.org/combine.specifications/omex-manifest">',
             f'  <content location="." format="{OMEX_FORMAT}"/>']
    for location, fmt, master in entries:
        m = ' master="true"' if master else ''
        lines.append(f'  <content location={quoteattr("./" + location)} format={quoteattr(fmt)}{m}/>')
    lines.append('</omexManifest>')
    return '\n'.join(lines) + '\n'


def _fmt(name):
    if name.endswith('_verification_config.json'):
        return VERIFICATION_CONFIG_FORMAT
    return FORMATS.get(os.path.splitext(name)[1].lower(), 'http://purl.org/NET/mediatypes/application/octet-stream')


def members(version, instance, work_dir):
    '''[(archive name, bytes, is_master)] in archive order.'''
    from cam_testing import harness
    stem = f'{version.vessel_type}_{version.name}_{instance.name}'
    model = harness.generate(version, work_dir, parameters=instance.parameters(), instance=instance)
    flat = model[:-len('.cellml')] + '_flat.cellml'
    if not os.path.isfile(flat):
        raise FileNotFoundError(f'{version.key}: libcuflynx wrote no flattened model ({flat})')
    res = os.path.join(work_dir, 'resources')
    prefix = version.id
    out = [(f'{stem}.cellml', open(flat, 'rb').read(), True)]

    def add(path, name=None):
        if path and os.path.isfile(path):
            out.append((name or os.path.basename(path), open(path, 'rb').read(), False))

    # the study, first: CUFLynx takes the first "obs" JSON and the first "param" CSV
    add(instance.obs_data_path)
    add(instance.params_for_id_path)
    # the inputs, as they are in the library
    add(instance.parameters_path)
    for p in sorted(os.listdir(instance.dir)):
        full = os.path.join(instance.dir, p)
        if os.path.isfile(full) and not p.endswith('.omex') and p not in {n for n, _, _ in out}:
            add(full)
    for p in (version.cellml_path, version.config_path, version.units_path, version.verification_config_path):
        add(p)
    # a supermodule's submodules' configs (recursively), so their required_citations travel with it
    seen = {version.key}
    for sub, _, _ in version.submodule_citations() if version.is_supermodule else []:
        if sub.key not in seen:
            seen.add(sub.key)
            add(sub.config_path)
    # the test network the model was generated from
    add(os.path.join(res, f'{prefix}_module_array.json'), f'{stem}_module_array.json')
    add(os.path.join(res, f'{prefix}_parameters.csv'), f'{stem}_model_parameters.csv')
    return out


def build(version, instance, work_dir=None, path=None):
    '''Writes the instance's .omex; returns its path.'''
    own = work_dir is None
    work_dir = work_dir or tempfile.mkdtemp(prefix=f'cam_omex_{version.id}_{instance.name}_')
    try:
        entries = members(version, instance, work_dir)
        path = path or omex_path(version, instance)
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as z:
            for name, data, _ in entries:
                z.writestr(name, data)
            z.writestr('manifest.xml', _manifest([(n, _fmt(n), m) for n, _, m in entries]))
        with open(path, 'wb') as f:
            f.write(buf.getvalue())
        return path
    finally:
        if own:
            shutil.rmtree(work_dir, ignore_errors=True)
