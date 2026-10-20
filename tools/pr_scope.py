'''
A pull request is a module PR, a migration PR or another PR (README "Pull requests").

A module PR's title starts with "Module:" and it changes only module content, so its review is a
review of modules alone: the library's module versions under modules/ (not modules/README.md or
modules/directory_schema.json), the system models in system_models/, the manifests generated from
them in manifests/ and the review records in reviews/. Any other PR (no prefix) changes only the
rest: cam_testing, tests, CI, tools, documentation, the schema, the setup. A migration PR's title
starts with "Migration:" and it may change both, for a format conversion that changes the module
content together with the code that reads it.

    python tools/pr_scope.py "<PR title>" <base ref>     # changed files from git diff <base>...HEAD

exits 1 and lists the files that don't belong: anything but module content in a module PR, module
content in any other PR, nothing in a migration PR.
'''
import subprocess
import sys

MODULE_PREFIX = 'Module:'
MIGRATION_PREFIX = 'Migration:'
# module content: these directories, apart from the repository files under them
MODULE_DIRS = ('modules/', 'system_models/', 'manifests/', 'reviews/')
NOT_MODULE_CONTENT = ('modules/README.md', 'modules/directory_schema.json')


def is_module_pr(title):
    return title.strip().startswith(MODULE_PREFIX)


def is_migration_pr(title):
    return title.strip().startswith(MIGRATION_PREFIX)


def is_module_content(path):
    return path.startswith(MODULE_DIRS) and path not in NOT_MODULE_CONTENT


def out_of_scope(title, paths):
    '''The changed paths a PR titled ``title`` may not change: for a module PR, every path that isn't
    module content; for a migration PR, none; for any other PR, every path that is module content.'''
    if is_migration_pr(title):
        return []
    module_pr = is_module_pr(title)
    return sorted(p for p in paths if is_module_content(p) != module_pr)


def changed_files(base):
    out = subprocess.run(['git', 'diff', '--name-only', f'{base}...HEAD'], check=True,
                         capture_output=True, text=True).stdout
    return [p for p in out.splitlines() if p]


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    title, base = argv[1], argv[2]
    bad = out_of_scope(title, changed_files(base))
    content = f'module content ({", ".join(MODULE_DIRS)}, apart from {", ".join(NOT_MODULE_CONTENT)})'
    if is_migration_pr(title):
        print('Migration PR: module content and the rest may both change.')
        return 0
    if is_module_pr(title):
        rule, ok = f'A module PR changes only {content}.', 'Module PR: every changed file is module content.'
    else:
        rule = (f'A PR whose title doesn\'t start with "{MODULE_PREFIX}" or "{MIGRATION_PREFIX}" changes '
                f'no {content}; module changes go in a module PR.')
        ok = 'Not a module PR: no changed file is module content.'
    if bad:
        print(f'{rule} Move these {len(bad)} file(s) to a PR of their own:')
        for p in bad:
            print(f'  {p}')
        return 1
    print(ok)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
