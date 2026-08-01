"""Bake every loading caveat into the published network files, once.

The working pickles carry four things a downstream user has to know about:
lengths in a placeholder voxel size (5 files), leftover max-flow solver
nodes, a mixed (x,y,z)/(z,y,x) coordinate convention, and a `capacity`
attribute inherited from ISAAC that is built on the EQUIVALENT diameter --
not the inscribed diameter the published min cuts use. This script resolves
all four and writes corrected copies under their publication names, so the
distributed dataset needs no corrections at load time.

Per file:
  1. drop non-integer nodes (`super_source`/`super_sink`) and reindex
  2. rescale lengths to the true voxel size where the stored one is a
     placeholder (lengths f, areas f^2, volumes f^3, *_size_factors f)
  3. reverse global_peak / geometric_centroid / local_peak (z,y,x)->(x,y,z)
     so EVERY vector attribute is (x,y,z)
  4. recompute `capacity` = (inscribed_diameter/2)^4 in m^4 on every throat,
     overwriting the equivalent-diameter values, so all files agree and the
     attribute reproduces the published min cuts
  5. stamp G.graph with the conventions, so the file is self-describing
  6. write with pickle protocol 4

RUN THIS FROM THE PYTHON 3.8 / NUMPY 1.x ENVIRONMENT. Arrays pickled by
numpy 2.x reference `numpy._core` and will not load under numpy 1.x; the
reverse direction is fine. Writing from py3.8 keeps both audiences working.

    python correct_networks.py --dest /path/to/staging
    python correct_networks.py --dest ... --only castlegate   # substring
    python correct_networks.py --dest ... --verify-only       # re-check

Writes `correction_report.csv` (per-file inventory + verification) next to
this script.
"""

import argparse
import csv
import gc
import os
import pickle
import shutil
import time

import numpy as np
import networkx as nx

FORMAT_VERSION = "1.0"
ZYX_FIELDS = ("global_peak", "geometric_centroid", "local_peak")

LENGTH_EXPONENT = {
    "coords": 1, "diameter": 1, "equivalent_diameter": 1,
    "extended_diameter": 1, "inscribed_diameter": 1, "max_size": 1,
    "geometric_centroid": 1, "global_peak": 1, "local_peak": 1,
    "surface_area": 2, "volume": 3, "region_volume": 3,
    "cross_sectional_area": 2, "direct_length": 1, "length": 1,
    "lens_volume": 3, "perimeter": 1, "spacing": 1, "total_length": 1,
    "total_volume": 3, "hydraulic_size_factors": 1,
    "diffusive_size_factors": 1,
}


def implied_voxel_um(G):
    """Voxel size the extraction used, from throat `global_peak`: those sit
    on the voxel grid, so the smallest non-zero spacing is one voxel."""
    gp = [d["global_peak"] for _, _, d in G.edges(data=True)
          if "global_peak" in d]
    if not gp:
        return float("nan")
    vals = np.unique(np.round(np.asarray(gp, dtype=float).ravel() * 1e12))
    diffs = np.diff(vals)
    diffs = diffs[diffs > 0]
    return float(diffs.min() / 1e6) if len(diffs) else float("nan")


def _plain(v):
    """numpy -> exact Python equivalent; item()/tolist() preserve dtype."""
    if isinstance(v, np.generic):
        return v.item()
    if isinstance(v, np.ndarray):
        return v.tolist()
    return v


def _rescale(attrs, f):
    for k, n in LENGTH_EXPONENT.items():
        if k in attrs:
            attrs[k] = np.asarray(attrs[k], dtype=float) * (f ** n)


def correct(path, true_vox_um, sample, extraction):
    """Load, correct, stamp. Returns (graph, info dict)."""
    with open(path, "rb") as fh:
        G = pickle.load(fh)

    info = {"in_pores": G.number_of_nodes(), "in_throats": G.number_of_edges()}

    # 1. solver nodes. convert_node_labels_to_integers() builds a whole
    # second copy of the graph, which is the difference between fitting in
    # RAM and not for the largest networks -- so only relabel when the node
    # set is not already 0..N-1.
    bad = [n for n in G.nodes() if not isinstance(n, (int, np.integer))]
    info["solver_nodes"] = ";".join(map(str, bad)) if bad else ""
    if bad:
        G.remove_nodes_from(bad)
    nodes = G.nodes()
    contiguous = (len(nodes) == 0
                  or (min(nodes) == 0 and max(nodes) == len(nodes) - 1))
    info["relabelled"] = not contiguous
    if not contiguous:
        G = nx.convert_node_labels_to_integers(G, ordering="sorted")
        gc.collect()

    # 2. voxel size
    stored = implied_voxel_um(G)
    info["stored_voxel_um"] = round(stored, 4) if stored == stored else ""
    if stored == stored and abs(stored - true_vox_um) / true_vox_um > 0.01:
        f = true_vox_um / stored
        info["rescale_factor"] = round(f, 6)
        for _, d in G.nodes(data=True):
            _rescale(d, f)
        for _, _, d in G.edges(data=True):
            _rescale(d, f)
    else:
        info["rescale_factor"] = 1.0

    # 3. coordinate order. .copy() matters: a[::-1] is a VIEW that keeps its
    # parent array alive, so storing views would hold two arrays per field
    # across millions of throats.
    n_flip = 0
    for _, d in G.nodes(data=True):
        for k in ZYX_FIELDS:
            if k in d:
                d[k] = np.asarray(d[k], dtype=float)[::-1].copy()
                n_flip += 1
    for _, _, d in G.edges(data=True):
        for k in ZYX_FIELDS:
            if k in d:
                d[k] = np.asarray(d[k], dtype=float)[::-1].copy()
                n_flip += 1
    info["fields_reordered"] = n_flip

    # 4. capacity, on the inscribed basis used by the published min cuts
    had, set_, missing = 0, 0, 0
    for _, _, d in G.edges(data=True):
        if "capacity" in d:
            had += 1
        di = d.get("inscribed_diameter")
        if di is None:
            d.pop("capacity", None)
            missing += 1
        else:
            d["capacity"] = float(di) ** 4 / 16.0     # (d/2)**4, m^4
            set_ += 1
    info["capacity_before"] = had
    info["capacity_after"] = set_
    info["throats_no_inscribed_d"] = missing

    # 5. plain Python attribute types. pickle memoises every numpy scalar
    # (they serialise through __reduce__) but writes Python floats inline
    # with no memo entry -- at 45M attributes that memo is the difference
    # between dumping and a MemoryError. It also drops the dependency on
    # numpy's array ABI, so the files load under any numpy, and shrinks
    # them by roughly a third. .item()/.tolist() are exact.
    for _, d in G.nodes(data=True):
        for k in list(d):
            d[k] = _plain(d[k])
    for _, _, d in G.edges(data=True):
        for k in list(d):
            d[k] = _plain(d[k])
    gc.collect()

    # 6. stamp
    G.graph = {
        "format_version": FORMAT_VERSION,
        "coord_order": "xyz",
        "attribute_types": "plain Python (float / int / list)",
        "units": "SI (metres)",
        "voxel_size_um": float(true_vox_um),
        "capacity_definition": "(inscribed_diameter/2)**4, m^4",
        "sample_name": sample,
        "extraction": extraction,
        "solver_nodes_removed": bool(bad),
    }

    info["out_pores"] = G.number_of_nodes()
    info["out_throats"] = G.number_of_edges()
    return G, info


def verify(path, true_vox_um):
    """Reload a corrected file and re-derive its conventions from scratch."""
    with open(path, "rb") as fh:
        G = pickle.load(fh)
    out = {"v_pores": G.number_of_nodes(), "v_throats": G.number_of_edges(),
           "v_version": G.graph.get("format_version", ""),
           "v_coord_order": G.graph.get("coord_order", "")}
    vox = implied_voxel_um(G)
    out["v_voxel_um"] = round(vox, 4) if vox == vox else ""
    out["v_voxel_ok"] = (vox == vox
                         and abs(vox - true_vox_um) / true_vox_um <= 0.01)
    out["v_nonint_nodes"] = sum(
        1 for n in G.nodes() if not isinstance(n, (int, np.integer)))
    # capacity must equal (inscribed/2)^4 on a sample of throats
    err, n = 0.0, 0
    for _, _, d in G.edges(data=True):
        if "capacity" in d and "inscribed_diameter" in d:
            p = float(d["inscribed_diameter"]) ** 4 / 16.0
            if p > 0:
                err = max(err, abs(d["capacity"] - p) / p)
            n += 1
            if n >= 20000:
                break
    out["v_capacity_max_relerr"] = "%.2e" % err if n else "n/a"
    del G
    gc.collect()
    return out


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=os.path.join(
        os.path.dirname(here), "Diamorse_persistence_pairs"))
    ap.add_argument("--dest", required=True)
    ap.add_argument("--names", default=os.path.join(here,
                                                    "network_names.csv"))
    ap.add_argument("--only", default=None,
                    help="only files whose publication name contains this")
    ap.add_argument("--verify-only", action="store_true")
    ap.add_argument("--skip-existing", action="store_true",
                    help="resume: leave already-corrected outputs alone")
    args = ap.parse_args()

    with open(args.names, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if args.only:
        rows = [r for r in rows if args.only in r["publication_name"]]
    if not args.verify_only:
        os.makedirs(args.dest, exist_ok=True)

    report = []
    for i, r in enumerate(rows, 1):
        name = r["publication_name"]
        src = os.path.join(args.source, r["original_stem"] + "_G.pickle")
        dst = os.path.join(args.dest, name + ".pickle")
        vox = float(r["voxel_size_um"])
        rec = {"publication_name": name, "sample": r["data_name"],
               "extraction": r["extraction"], "true_voxel_um": vox}

        if args.verify_only:
            if os.path.exists(dst):
                rec.update(verify(dst, vox))
                print("[%2d/%d] %-42s verified" % (i, len(rows), name[:42]),
                      flush=True)
            report.append(rec)
            continue

        if not os.path.exists(src):
            print("[%2d/%d] %-42s SOURCE MISSING" % (i, len(rows), name[:42]),
                  flush=True)
            rec["error"] = "source missing"
            report.append(rec)
            continue

        if args.skip_existing and os.path.exists(dst):
            rec.update(verify(dst, vox))
            print("[%2d/%d] %-42s already done (%s)"
                  % (i, len(rows), name[:42],
                     "OK" if rec.get("v_voxel_ok") else "*** CHECK ***"),
                  flush=True)
            report.append(rec)
            continue

        t0 = time.time()
        print("[%2d/%d] %-42s ..." % (i, len(rows), name[:42]),
              end=" ", flush=True)
        G, info = correct(src, vox, r["data_name"], r["extraction"])
        rec.update(info)
        gc.collect()
        try:
            with open(dst, "wb") as fh:
                pickle.dump(G, fh, protocol=4)
        except (MemoryError, OSError) as exc:
            # never leave a truncated file behind for the next run to trust
            del G
            gc.collect()
            if os.path.exists(dst):
                os.remove(dst)
            print("FAILED (%s) -- partial output removed"
                  % type(exc).__name__, flush=True)
            rec["error"] = "%s during write" % type(exc).__name__
            report.append(rec)
            continue
        del G
        gc.collect()

        rec.update(verify(dst, vox))
        rec["mb"] = round(os.path.getsize(dst) / 1048576.0, 1)
        rec["seconds"] = round(time.time() - t0, 1)
        print("%d pores %d throats | vox %s->%.4f | cap %d->%d | %s | %.0fs"
              % (rec["out_pores"], rec["out_throats"],
                 rec.get("stored_voxel_um"), vox,
                 rec.get("capacity_before", 0), rec.get("capacity_after", 0),
                 "OK" if rec.get("v_voxel_ok") else "*** CHECK ***",
                 rec["seconds"]), flush=True)
        report.append(rec)

    if report:
        keys = []
        for rec in report:
            for k in rec:
                if k not in keys:
                    keys.append(k)
        out = os.path.join(here, "correction_report.csv")
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(report)
        print("\nreport:", out)

    if not args.verify_only:
        for f in ("network_names.csv", "network_to_openpnm.py", "README.md"):
            p = os.path.join(here, f)
            if os.path.exists(p):
                shutil.copy2(p, os.path.join(args.dest, f))
        print("copied loader + metadata into", args.dest)

    bad = [r["publication_name"] for r in report
           if "error" in r or r.get("v_voxel_ok") is False
           or r.get("v_nonint_nodes")]
    print("\n%d file(s) processed, %d need attention%s"
          % (len(report), len(bad), (": " + ", ".join(bad)) if bad else ""))


if __name__ == "__main__":
    main()
