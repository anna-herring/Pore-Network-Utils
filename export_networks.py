"""Copy the extraction outputs to their published names, ready for upload.

The working files are named after the segmentation/processing chain that
produced them (e.g. ``segmented_Jangda_04_Dry_Segmented_IMP_RPI_3_G.pickle``).
This copies each one to the name used in the papers and figures
(``clashach_sandstone_diamorse.pickle``), driven by ``network_names.csv``.

COPIES rather than renames on purpose: the analysis code matches these
networks to their persistence-pair CSVs and to the voxel-size spreadsheet by
filename prefix, so renaming in place would silently break it. The originals
stay where they are.

    python export_networks.py --dest /path/to/upload_staging
    python export_networks.py --dest ... --dry-run     # list, copy nothing
    python export_networks.py --dest ... --link        # hardlink (no space)

Writes ``network_names.csv`` alongside the copies so the mapping travels
with the data. ``--link`` needs the destination on the same volume as the
source and uses no extra disk; it is the right choice when staging tens of
GB, but edit/delete of one name affects the other only for deletes.
"""

import argparse
import csv
import os
import shutil


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=os.path.join(
        os.path.dirname(here), "Diamorse_persistence_pairs"),
        help="folder holding the *_G.pickle files")
    ap.add_argument("--dest", required=True, help="staging folder to create")
    ap.add_argument("--names", default=os.path.join(here,
                                                    "network_names.csv"))
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--link", action="store_true",
                    help="hardlink instead of copying (same volume only)")
    args = ap.parse_args()

    with open(args.names, newline="") as fh:
        rows = list(csv.DictReader(fh))

    if not args.dry_run:
        os.makedirs(args.dest, exist_ok=True)

    total = 0.0
    missing = []
    for r in rows:
        src = os.path.join(args.source, r["original_stem"] + "_G.pickle")
        dst = os.path.join(args.dest, r["publication_name"] + ".pickle")
        if not os.path.exists(src):
            missing.append(r["original_stem"])
            print("  MISSING   %s" % r["original_stem"][:66])
            continue
        mb = os.path.getsize(src) / 1048576.0
        total += mb
        flag = "  rescale" if r["rescale_needed"].strip().upper() == "YES" \
            else ""
        print("  %-44s %8.1f MB%s" % (r["publication_name"], mb, flag))
        if args.dry_run:
            continue
        if os.path.exists(dst):
            continue
        if args.link:
            try:
                os.link(src, dst)
                continue
            except OSError:
                pass                       # cross-volume: fall back to copy
        shutil.copy2(src, dst)

    print("\n%d files, %.1f GB total" % (len(rows) - len(missing),
                                         total / 1024.0))
    if missing:
        print("%d MISSING - check the source folder:" % len(missing))
        for m in missing:
            print("   ", m)
    if not args.dry_run:
        shutil.copy2(args.names, os.path.join(args.dest, "network_names.csv"))
        print("copied network_names.csv into", args.dest)
    else:
        print("(dry run - nothing written)")


if __name__ == "__main__":
    main()
