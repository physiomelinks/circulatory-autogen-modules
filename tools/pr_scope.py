'''
A pull request is either a module PR or another PR (README "Pull requests").

A module PR's title starts with "Module:" and it changes only module content, so its review is a
review of modules alone: the library's module versions under modules/ (not modules/README.md or
modules/directory_schema.json), the system models in system_models/, the manifests generated from
them in manifests/ and the review records in reviews/. Every other change (cam_testing, tests, CI,
tools, documentation, the schema, the setup) goes in a PR of its own, without the prefix.

    python tools/pr_scope.py "<PR title>" <base ref>     # changed files from git diff <base>...HEAD

exits 1 and lists the files that don't belong when a module PR changes anything but module content.
'''
import subprocess
import sys

MODULE_PREFIX = 'Module:'
# module content: these directories, apart from the repository files under them
MODULE_DIRS = ('modules/', 'system_models/', 'manifests/', 'reviews/')
NOT_MODULE_CONTENT = ('modules/README.md', 'modules/directory_schema.json')


def is_module_pr(title):
    return title.strip().startswith(MODULE_PREFIX)


def is_module_content(path):
    return path.startswith(MODULE_DIRS) and path not in NOT_MODULE_CONTENT


def out_of_scope(title, paths):
    '''The changed paths a PR titled ``title`` may not change: for a module PR, every path that isn't
    module content; nothing for any other PR.'''
    if not is_module_pr(title):
        return []
    return sorted(p for p in paths if not is_module_content(p))


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
    if not is_module_pr(title):
        print(f'Not a module PR (the title doesn\'t start with "{MODULE_PREFIX}"): any files may change.')
        return 0
    if bad:
        print(f'A module PR changes only module content ({", ".join(MODULE_DIRS)}, apart from '
              f'{", ".join(NOT_MODULE_CONTENT)}). Move these {len(bad)} file(s) to a PR of their own:')
        for p in bad:
            print(f'  {p}')
        return 1
    print('Module PR: every changed file is module content.')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
