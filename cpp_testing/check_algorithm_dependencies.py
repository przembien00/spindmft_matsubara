#!/usr/bin/env python3
"""Reject source/build dependencies between algorithm directories.

Run from any directory. Shared libraries are also checked for dependencies
back into algorithms. Comparisons of simulation outputs and optional input
data are outside this source/build boundary.
"""
import argparse
import os
from pathlib import Path
import re


SOURCE_SUFFIXES = {'.h', '.hpp', '.hh', '.hxx', '.c', '.cc', '.cpp', '.cxx', '.cmake'}
SKIP = {'Data', 'Plots', 'Notes', 'logs', '.git', '__pycache__', 'CMakeFiles'}


def algorithm_roots(root):
    candidates = list(root.iterdir())
    archive = root / 'deprecated_algorithms'
    if archive.is_dir():
        candidates.extend(archive.iterdir())
    return sorted(p for p in candidates if p.is_dir()
                  and ((p / 'Algorithm').is_dir() or p.name == 'ExactDiagonalization'))


def source_files(folder):
    for directory, subdirs, files in os.walk(folder):
        subdirs[:] = [name for name in subdirs if name not in SKIP
                      and not name.startswith(('build', 'cmake-build'))]
        for name in files:
            path = Path(directory) / name
            if path.suffix in SOURCE_SUFFIXES or name.startswith(('CMakeLists', 'Makefile')):
                yield path
        # Do not silently allow a source directory symlink into another algorithm.
        for name in subdirs:
            path = Path(directory) / name
            if path.is_symlink():
                yield path


def violations(root):
    root = root.resolve()
    algorithms = algorithm_roots(root)
    shared = [root / name for name in ('cpp_libs', 'python_libs')
              if (root / name).is_dir()]
    errors = []
    for owner in algorithms + shared:
        others = [p for p in algorithms if p != owner]
        if not others:
            continue
        foreign_path = re.compile(r'(?<![\w-])(?:'
                                  + '|'.join(re.escape(p.name) for p in others) + r')/')
        for path in source_files(owner):
            resolved = path.resolve()
            if any(resolved == other or other in resolved.parents for other in others):
                errors.append(f'{path.relative_to(root)}: source symlink points into another algorithm')
                continue
            if path.is_dir():
                continue
            build_file = path.suffix == '.cmake' or path.name.startswith(('CMakeLists', 'Makefile'))
            for number, line in enumerate(path.read_text().splitlines(), 1):
                if build_file:
                    line = line.split('#', 1)[0]
                elif not re.match(r'\s*#\s*include\b', line):
                    continue
                if foreign_path.search(line):
                    errors.append(f'{path.relative_to(root)}:{number}: {line.strip()}')
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    errors = violations(root)
    if errors:
        print('Cross-algorithm source/build dependencies:\n' + '\n'.join(errors))
        return 1
    print(f'No cross-algorithm source/build dependencies in {len(algorithm_roots(root))} algorithms.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
