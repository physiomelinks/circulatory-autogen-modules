"""
Parameters and params_for_id files written by the tools, in the forms libcuflynx reads: nested JSON
parameters (``libcuflynx.utilities.parameters_file``) and flat JSON params_for_id
(``libcuflynx.utilities.params_for_id_file``). The tools build rows; these write them.

* ``write_model_parameters``: a system model's ``<model>_parameters.json`` from rows named
  ``<vessel>/<variable>``, or ``{variable}_{vessel}`` as the older tools write them, which are
  rewritten against the model's module array by libcuflynx's ``parameter-names``;
* ``write_parameterisation``: a parameterisation's ``<p>_parameters.json`` (rooted at its module,
  globals ``global/<variable>``);
* ``write_module_table`` / ``read_module_table``: the per-module parameters table of the
  pre-versions layout (``modules/<name>/<name>_parameters.json``), one row per (vessel_type,
  BC_type, variable): nested ``{vessel_type: {BC_type: {variable: ...}}}`` (``{vessel_type:
  {variable: ...}}`` with no BC_type) with the globals under ``global``;
* ``write_params_for_id``: a params_for_id file from entries, or from older CSV rows.
"""
import contextlib
import io
import os

from libcuflynx.utilities import parameters_file as pf
from libcuflynx.utilities import params_for_id_file as pfi

#: the columns a pre-versions module table row has besides those a parameter keeps
TABLE_KEYS = ('vessel_type', 'BC_type')


def write_model_parameters(path, rows, module_array=None, module_library_dirs=()):
    """Write ``rows`` as the model's nested parameters file ``path``. A name without a ``/`` is
    first read as a global; with ``module_array`` given, a ``{variable}_{vessel}`` name of one of
    its vessels is then rewritten ``<vessel>/<variable>`` (``cuflynx-upgrade parameter-names``,
    checked against the modules of ``module_library_dirs``, without libcuflynx's built-in
    ones)."""
    pf.write_parameter_rows(path, rows, pf.MODEL)
    if module_array is not None:
        from libcuflynx.scripts.migrate_parameter_names import main as parameter_names
        args = ['--module-array', str(module_array), '--parameters', str(path)]
        for directory in module_library_dirs:
            # the library has its own copies of libcuflynx's built-in modules
            args += ['--module-library-dir', str(directory), '--no-builtin-modules']
        with contextlib.redirect_stdout(io.StringIO()):
            parameter_names(args)
    return path


def write_parameterisation(path, rows, global_names=()):
    """Write ``rows`` as the parameterisation file ``path`` (rooted at its module); a row whose
    name is in ``global_names`` is a global (``global/<variable>``)."""
    rows = [dict(r, variable_name=f'{pf.GLOBAL}/{r["variable_name"]}')
            if r.get('variable_name') in set(global_names) else r for r in rows]
    return pf.write_parameter_rows(path, rows, pf.MODULE)


def write_module_table(path, rows):
    """Write the pre-versions per-module table ``rows`` (``vessel_type``, ``BC_type``,
    ``variable_name``, ...; ``vessel_type`` "global" for a global) as nested JSON."""
    named = []
    for row in rows:
        leaf = {k: v for k, v in row.items() if k not in TABLE_KEYS}
        if row.get('vessel_type') == 'global':
            leaf['variable_name'] = f'{pf.GLOBAL}/{row["variable_name"]}'
        else:
            leaf['variable_name'] = '/'.join([row['vessel_type']] + (
                [row['BC_type']] if row.get('BC_type') else []) + [row['variable_name']])
        named.append(leaf)
    return pf.write_parameter_rows(path, named, pf.MODEL)


def read_module_table(path):
    """The rows of a pre-versions per-module table written by ``write_module_table``."""
    rows = []
    for row in pf.read_parameter_rows(path, kind=pf.MODEL):
        parts = row['variable_name'].split('/')
        if len(parts) == 1:
            vessel_type, bc_type, variable = 'global', '', parts[0]
        elif len(parts) == 2:
            (vessel_type, variable), bc_type = parts, ''
        else:
            vessel_type, bc_type, variable = parts
        rows.append(dict(row, vessel_type=vessel_type, BC_type=bc_type, variable_name=variable))
    return rows


def write_params_for_id(path, entries):
    """Write a params_for_id JSON file from ``entries`` (its entries), or from CSV text or rows
    (``[header, row, ...]`` lists of cells) in the older ``vessel_name,param_name,...`` form."""
    if isinstance(entries, list) and entries and isinstance(entries[0], list):
        entries = '\n'.join(','.join(str(c) for c in row) for row in entries) + '\n'
    doc = pfi.csv_to_document(entries) if isinstance(entries, str) else entries
    return pfi.write_params_for_id(str(os.path.splitext(str(path))[0] + pfi.EXTENSION), doc)
