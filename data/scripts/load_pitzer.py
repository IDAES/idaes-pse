"""
Simple script to load Pitzer data.

This script should only be used if the Pitzer data is not already in the DMF.

You can check this from the command-line by looking at the output
of the `dmf ls -s name -s type` command and looking for a resource with the
name "pitzer-paper" and type "publication". You can then double-check that all
the data is linked to this publication with `dmf related <id>`, where `<id>` is the
identifier associated with the publication resource you just listed.

The function `get_pitzer_data()` in this file is an example of how to fetch
the data tables.
"""
# stdlib
import argparse
import json
import logging
from pathlib import Path
# package
from idaes.dmf.datasets import PublicationDataset


_log = logging.getLogger(__name__)
_hnd = logging.StreamHandler()
_hnd.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(msg)s"))
_log.addHandler(_hnd)


def set_log_level(vb):
    if vb > 1:
        level = logging.DEBUG
    elif vb > 0:
        level = logging.INFO
    else:
        level = logging.WARNING

    _log.setLevel(level)
    logging.getLogger("idaes.dmf").setLevel(level)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--datadir",
        "-d",
        default=None,
        help="Data & configuration directory (default is current working directory)",
    )
    p.add_argument(
        "--workspace",
        "-w",
        default=None,
        help="Workspace directory (default is idaes-pse data workspace)",
    )
    p.add_argument(
        "--verbose", "-v", action="count", dest="vb", default=0,
        help="Increase logging verbosity (repeatable)"
    )

    args = p.parse_args()

    set_log_level(args.vb)

    if args.workspace is None:
        workspace = None
    else:
        workspace = Path(args.workspace).absolute()
    _log.info(f"Workspace directory is {workspace}")

    if args.datadir is None:
        data_directory = Path(".").absolute()
    else:
        data_directory = Path(args.datadir).absolute()
    _log.info(f"Data directory is {data_directory}")

    pd = PublicationDataset(workspace=workspace)
    pd.load(data_directory)

