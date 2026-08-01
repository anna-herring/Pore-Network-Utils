"""Write SHA256 checksums for the published network files.

Output is the standard `sha256sum` format -- lowercase hex, two spaces,
filename -- so a downloader can verify with the tool they already have:

    sha256sum -c CHECKSUMS.txt          # Linux / macOS
    python make_checksums.py --check    # anywhere

On Windows without sha256sum:

    Get-FileHash file.pickle -Algorithm SHA256

Usage:
    python make_checksums.py --dir /path/to/staging
    python make_checksums.py --dir ... --check
"""

import argparse
import hashlib
import os

CHUNK = 1 << 22          # 4 MiB


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", default="CHECKSUMS.txt")
    ap.add_argument("--check", action="store_true",
                    help="verify against an existing CHECKSUMS.txt")
    ap.add_argument("--ext", default=".pickle",
                    help="only checksum files with this extension; '' for "
                         "all. Defaults to the data files only -- README and "
                         "scripts get revised, and a checksum failure should "
                         "mean a bad download, not an edited doc.")
    args = ap.parse_args()

    out = os.path.join(args.dir, args.out)

    if args.check:
        bad = missing = ok = 0
        with open(out) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                want, name = line.split("  ", 1)
                p = os.path.join(args.dir, name)
                if not os.path.exists(p):
                    print("MISSING  %s" % name)
                    missing += 1
                elif sha256(p) == want:
                    ok += 1
                else:
                    print("FAILED   %s" % name)
                    bad += 1
        print("\n%d OK, %d failed, %d missing" % (ok, bad, missing))
        return

    names = sorted(f for f in os.listdir(args.dir)
                   if not f.startswith(".") and f != args.out
                   and (not args.ext or f.endswith(args.ext))
                   and os.path.isfile(os.path.join(args.dir, f)))
    total = 0
    lines = []
    for n in names:
        p = os.path.join(args.dir, n)
        mb = os.path.getsize(p) / 1048576.0
        total += mb
        lines.append("%s  %s" % (sha256(p), n))
        print("  %-46s %8.1f MB" % (n[:46], mb), flush=True)

    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n%d files, %.2f GB -> %s" % (len(names), total / 1024.0, out))


if __name__ == "__main__":
    main()
