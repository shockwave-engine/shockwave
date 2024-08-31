#! /usr/bin/env python

"""
Copies the compile commands, keeping commands the are not in destination.

Without this, building tests produces the correct compile commands for test files, but then building the lib erases them.
"""

import json
import os
from argparse import ArgumentParser


def load_commands(file: str) -> dict[str, dict]:
    with open(file, 'r') as f:
        content: list = json.load(f)
        return {obj['file']: obj for obj in content}


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--from', dest='src', help='source commands')
    parser.add_argument('--to', help='destination commands')
    args = parser.parse_args()

    if not os.path.isfile(args.src):
        print(f"Not copying {args.src} since it doesn't exist")
        exit(0)

    print(f"Copying {args.src} to {args.to}")

    src = load_commands(args.src)
    dst = load_commands(args.to) if os.path.isfile(args.to) else {}

    dst.update(src)

    with open(args.to, 'w') as f:
        json.dump(list(dst.values()), f, indent=4)
