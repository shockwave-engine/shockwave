#!/usr/bin/env python

"""
This is to simplify the development flow for developpers.
Consumers do not need this, use CMake directly.
"""

import os
import sys
import subprocess

from argparse import ArgumentParser, Namespace


def project_root() -> str:
    return os.path.dirname(sys.argv[0])


def build_dir(build_type: str) -> str:
    return f'{project_root()}/out/{build_type}'


def banner(title: str, cmd: list[str]):
    print("\n------------------------------------------------------------")
    print(f"-- {title}")
    print(f"-- Executing: {' '.join(cmd)}")
    print("------------------------------------------------------------\n")


def execute(title: str, cmd: list[str], **kwargs):
    banner(title, cmd)
    subprocess.run(cmd, **kwargs)


def contains_any(lst: list[str], elements: list[str]) -> bool:
    return len(set(elements).intersection(lst)) > 0


def option(name: str, value: bool) -> str:
    return f'-D{name}:BOOL={'ON' if value else 'OFF'}'


def needs_coverage(targets: list[str]) -> bool:
    return contains_any(targets, ['shockwave-test-coverage', 'shockwave-test-coverage-viewer'])


def needs_tests(targets: list[str]) -> bool:
    return needs_coverage(targets) or contains_any(targets, ['shockwave-test'])


def required_options(targets: list[str]) -> list[str]:
    return [
        option('SHOCKWAVE_BUILD_TESTS', needs_tests(targets)),
        option('SHOCKWAVE_TEST_COVERAGE', needs_coverage(targets))
    ]


def build(args: Namespace):
    execute('Configuring', ['cmake',
                            '-B', build_dir(args.build_type),
                            f'-DCMAKE_BUILD_TYPE={args.build_type}',
                            *required_options(args.targets)
                            ])

    execute('Building', ['cmake',
                         '--build', build_dir(args.build_type),
                         '--config', f'{args.build_type}',
                         '--target', *args.targets])


def run(args: Namespace):
    build(args)
    for target in args.targets:
        execute(f"Running target {target}", [
                f'{build_dir(args.build_type)}/{target}'
                ])


def lint(args: Namespace):
    # TODO:
    ...


def add_build_args(parser: ArgumentParser, targets=['shockwave', 'shockwave-test', 'shockwave-test-coverage', 'shockwave-test-coverage-viewer', 'clean', 'all']):
    parser.add_argument('targets',
                        choices=targets,
                        default='all',
                        nargs='*',
                        help='The target(s) to build')
    parser.add_argument('--build-type',
                        choices=['Debug', 'Release',
                                 'RelWithDebInfo', 'MinSize'],
                        default='Debug',
                        help='The CMake build type to build')


def parse_args(args: list[str]) -> Namespace:
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(required=True)

    build_parser = subparsers.add_parser('build')
    add_build_args(build_parser)
    build_parser.set_defaults(func=build)

    run_parser = subparsers.add_parser('run')
    add_build_args(run_parser, targets=['shockwave-test'])
    run_parser.set_defaults(func=run)

    parsed = parser.parse_args(args)

    # Cannot set default as list
    if parsed.targets is not None and isinstance(parsed.targets, str):
        parsed.targets = [parsed.targets]

    return parsed


def main(argv):
    args = parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
