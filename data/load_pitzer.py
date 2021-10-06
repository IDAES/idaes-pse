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
import glob
import logging
from pathlib import Path
# package
from idaes.dmf import resource as rsrc

paper_meta = {
    "date": "1984",  # DMF parses dates automatically
    "doi": "10.1063/1.555709",
    "language": "english",
    "source": "Kenneth S. Pitzer, J. Christopher Peiper, and R. H. Busey, "
    '"Thermodynamic Properties of Aqueous Sodium Chloride Solutions". '
    "Journal of Physical and Chemical Reference Data 13, 1 (1984)",
}

PITZER_TAG = "pitzer"
PITZER_NAME = "pitzer-paper"

_log = logging.getLogger(__name__)
_hnd = logging.StreamHandler()
_hnd.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(msg)s"))
_log.addHandler(_hnd)


def pitzer_paper():
    r = rsrc.Resource(type_=rsrc.ResourceTypes.publication, name=PITZER_NAME)
    r.sources.append(paper_meta)
    return r


def pitzer_table(path):
    r = rsrc.Resource(type_=rsrc.ResourceTypes.tabular)
    r.add_table(path, do_copy=False)
    r.add_tag(PITZER_TAG)
    return r


def add_pitzer_data(dmf, directory) -> int:
    n = 0
    root = pitzer_paper()
    root.add_data_file(directory / "Pitzer_1984.pdf", do_copy=False)
    dmf.add(root)
    for csv_file in directory.glob("pitzer_*.csv"):
        ptable = pitzer_table(csv_file)
        ptable.add_tag(PITZER_TAG)
        dmf.add(ptable)
        rsrc.create_relation(root, rsrc.Predicates.derived, ptable)
        n += 1
    dmf.update()
    return n


def get_pitzer_data(dmf, include_pub=True):
    pub = dmf.find_one(name=PITZER_NAME)
    if pub:
        if include_pub:
            yield pub
        _log.debug("Finding data files attached to Pitzer publication..")
        n = 0
        for r in dmf.find_related_resources(pub, outgoing=True):
            n += 1
            yield r
        _log.debug(f"Found {n} data files")
    else:
        _log.debug("Pitzer publication not found")
        if include_pub:
            yield None


def remove_pitzer_data(dmf) -> int:
    total = 0
    while True:
        _log.debug("Try to remove old Pitzer data ..")
        n, pub = 0, None
        for r in get_pitzer_data(dmf):
            if n == 0:
                pub = r
            else:
                _log.debug(f"Removing tabular resource ({r.id})")
                dmf.remove(r.id, update_relations=False)
            n += 1
        if not pub:
            break
        if pub:
            _log.debug(f"Removing publication resource ({pub.id})")
            dmf.remove(pub.id, update_relations=False)
        total += n
        dmf.update()
    return total


def main(workspace=None, data_dir=None):
    from idaes.dmf import DMF

    dmf = DMF(workspace)
    _log.debug("Removing previous Pitzer data (if any)..")
    n = remove_pitzer_data(dmf)
    if n == 0:
        _log.debug("Nothing to remove")
    else:
        _log.debug(f"Removed {n} resources from DMF")
    _log.debug("Adding new Pitzer data..")
    n = add_pitzer_data(dmf, data_dir)
    _log.info(f"Added {n} tables, and one publication record, to the DMF")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--datadir",
        "-d",
        default=None,
        help="Data directory (default is current working directory)",
    )
    p.add_argument(
        "--workspace",
        "-w",
        default=None,
        help="Workspace directory (default is current working directory)",
    )
    args, params = p.parse_args(), {}

    _log.setLevel(logging.INFO)

    if args.workspace is None:
        params["workspace"] = Path(".").absolute()
    else:
        params["workspace"] = Path(args.workspace).absolute()
    _log.info(f"Workspace directory is {params['workspace']}")

    if args.datadir is None:
        params["data_dir"] = Path(".").absolute()
    else:
        params["data_dir"] = Path(args.datadir).absolute()
    _log.info(f"Data directory is {params['data_dir']}")

    main(**params)
