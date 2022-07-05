#!/usr/bin/env python
"""
Extract JSON schemas into files.
"""
import argparse
import json
import os
import sys

#
from idaes.dmf import resource

DMF_SOURCE_DIR = os.path.join("source", "dmf")


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "-d",
        "--dir",
        dest="dir",
        default="schemas",
        help="Destination directory "
        "for .json files. If absolute, used as-is. "
        'If relative, added to "{}". (default=schemas)'.format(DMF_SOURCE_DIR),
    )
    args = p.parse_args()
    if os.path.isabs(args.dir):
        schema_dir = args.dir
    else:
        schema_dir = os.path.join(DMF_SOURCE_DIR, args.dir)
    for schema, name in ((resource.RESOURCE_SCHEMA, "resource"),):
        json_file = os.path.join(schema_dir, name) + ".json"
        json_fp = open(json_file, "w")
        json.dump(schema, json_fp, indent=4)
        print("Wrote {}".format(json_file))
    return 0


if __name__ == "__main__":
    sys.exit(main())
