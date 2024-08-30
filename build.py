#!/usr/bin/env python

"""
This is to simplify executing CMake commands for developpers.
Consumers should use CMake directly
"""

import os
import sys

from subprocess import run
from argparse import ArgumentParser, Namespace


def project_root() -> str:
    return os.path.dirname(sys.argv[0])


def build_dir(build_type: str) -> str:
    return f'{project_root()}/out/{build_type}'

# TODO
# def required_options(targets: list[str]) -> list[str]:


def parse_args(args: list[str]) -> Namespace:
    parser = ArgumentParser()

    parser.add_argument('targets',
                        choices=['shockwave', 'shockwave-test',
                                 'shockwave-test-coverage', 'clean', 'all'],
                        default='all',
                        nargs='*',
                        help='The target(s) to build')

    parser.add_argument('--build-type',
                        choices=['Debug', 'Release',
                                 'RelWithDebInfo', 'MinSize'],
                        default='Debug',
                        help='The CMake build type to build')

    parsed = parser.parse_args(args)

    # Cannot set default as list
    if isinstance(parsed.targets, str):
        parsed.targets = [parsed.targets]

    return parsed


def main(argv):
    args = parse_args(argv)

    print("Configuring...")
    run(['cmake',
         '-B', build_dir(args.build_type),
        f'-DCMAKE_BUILD_TYPE={args.build_type}',
         '-DSHOCKWAVE_BUILD_TESTS:BOOL=ON',
         '-DSHOCKWAVE_TEST_COVERAGE:BOOL=ON'
         ])

    print("Building...")
    run(['cmake',
         '--build', build_dir(args.build_type),
         '--config', f'{args.build_type}',
         '--target', *args.targets])


if __name__ == "__main__":
    main(sys.argv[1:])
