"""Load a pore-network pickle from this dataset and convert it to OpenPNM.

The networks are distributed as pickled NetworkX ``Graph`` objects: one node
per pore, one edge per throat, with OpenPNM-style property names carried as
node/edge attributes (the ``pore.``/``throat.`` prefixes are stripped). This
script restores those prefixes and rebuilds an ``openpnm.network.Network``.

    from network_to_openpnm import load_graph, to_openpnm

    G  = load_graph("clashach_sandstone_diamorse.pickle", voxel_size_um=7.0)
    pn = to_openpnm(G)
    print(pn)

Or from the command line::

    python network_to_openpnm.py clashach_sandstone_diamorse.pickle \
        --voxel-um 7.0 --save-vtk

CONVENTIONS
-----------
Published files carry ``G.graph["format_version"] = "1.0"``. In those files
every correction below has already been applied, so you can simply load and
use them -- ``load_graph`` and ``to_openpnm`` detect the stamp and become
no-ops. The list matters only if you are handling a pre-1.0 file.

    G.graph == {
        "format_version": "1.0",
        "coord_order": "xyz",
        "units": "SI (metres)",
        "voxel_size_um": <true voxel size, already applied>,
        "capacity_definition": "(inscribed_diameter/2)**4, m^4",
        "sample_name": ..., "extraction": ...,
    }

1. UNITS ARE METRES. Every length/area/volume attribute is SI, not voxels
   and not microns. Multiply by 1e6 for microns.

2. COORDINATE ORDER. In 1.0 files EVERY vector attribute is ``(x, y, z)``.
   In pre-1.0 files node ``coords`` was ``(x, y, z)`` but
   ``global_peak``/``geometric_centroid``/``local_peak`` were ``(z, y, x)``
   (the porespy region convention); ``to_openpnm`` reverses those. It keys
   off ``coord_order``, so passing a 1.0 file through does NOT double-flip.

3. VOXEL SIZE. Some extractions recorded a placeholder voxel size, leaving
   absolute lengths off by a constant factor -- topology right, permeability
   and capillary pressures wrong. 1.0 files are already corrected to the
   true voxel size recorded in the stamp. For a pre-1.0 file, pass the
   correct ``voxel_size_um`` (see ``network_names.csv``) and ``load_graph``
   rescales: lengths by f, areas by f^2, volumes by f^3, and the OpenPNM
   ``*_size_factors`` (units of length) by f. It is safe to pass every time
   -- a no-op when the stored size is already right.

4. SOLVER NODES. Some pre-1.0 files were saved after a max-flow run and
   still contain ``super_source``/``super_sink`` nodes with no geometry.
   ``load_graph`` drops any non-integer node and reindexes contiguously.
   1.0 files contain none.

5. THROAT ``capacity``. In 1.0 files this is ``(inscribed_diameter/2)**4``
   in m^4 on every throat -- the definition behind the published minimum
   cuts. PRE-1.0 FILES CARRY A DIFFERENT QUANTITY: an inherited
   ``(equivalent_diameter/2)**4``, which yields a different min cut. If you
   are working from a pre-1.0 file, recompute it from
   ``inscribed_diameter`` rather than trusting the stored value.

Requires: networkx, numpy, openpnm >= 3.0 (only for ``to_openpnm``).
"""

import argparse
import pickle

import numpy as np
import networkx as nx

# node/edge attributes stored (z, y, x) rather than (x, y, z)
ZYX_FIELDS = {"global_peak", "geometric_centroid", "local_peak"}

# attributes whose units are length^n, for voxel-size rescaling.
# OpenPNM's *_size_factors have units of LENGTH (they are area/length), so
# they take exponent 1, not 3.
LENGTH_EXPONENT = {
    # pores
    "coords": 1, "diameter": 1, "equivalent_diameter": 1,
    "extended_diameter": 1, "inscribed_diameter": 1, "max_size": 1,
    "geometric_centroid": 1, "global_peak": 1, "local_peak": 1,
    "surface_area": 2, "volume": 3, "region_volume": 3,
    # throats
    "cross_sectional_area": 2, "direct_length": 1, "length": 1,
    "lens_volume": 3, "perimeter": 1, "spacing": 1, "total_length": 1,
    "total_volume": 3, "hydraulic_size_factors": 1,
    "diffusive_size_factors": 1,
    # capacity is (d/2)^4 in this dataset's own min-cut convention
    "capacity": 4,
}


def implied_voxel_size_um(G):
    """Infer the voxel size the extraction used, from the throat
    ``global_peak`` values: those sit on the voxel grid, so the smallest
    non-zero spacing between distinct coordinates is one voxel."""
    gp = np.array([d["global_peak"] for _, _, d in G.edges(data=True)
                   if "global_peak" in d], dtype=float)
    if not len(gp):
        return float("nan")
    vals = np.unique(np.round(gp.ravel() * 1e12))       # pm, kills fp noise
    diffs = np.diff(vals)
    diffs = diffs[diffs > 0]
    return float(diffs.min() / 1e6) if len(diffs) else float("nan")


def load_graph(path, voxel_size_um=None, verbose=True):
    """Load a network pickle, drop solver nodes, and (optionally) correct
    the voxel size. ``voxel_size_um`` is the TRUE voxel size in microns --
    look it up in ``network_names.csv``."""
    with open(path, "rb") as fh:
        G = pickle.load(fh)

    bad = [n for n in G.nodes() if not isinstance(n, (int, np.integer))]
    if bad:
        if verbose:
            print("dropping %d non-pore node(s): %s" % (len(bad), bad[:4]))
        G.remove_nodes_from(bad)
    G = nx.convert_node_labels_to_integers(G, ordering="sorted")

    if voxel_size_um is not None:
        implied = implied_voxel_size_um(G)
        if np.isnan(implied):
            if verbose:
                print("cannot infer the stored voxel size; no rescale applied")
        elif abs(implied - voxel_size_um) / voxel_size_um > 0.01:
            f = voxel_size_um / implied
            if verbose:
                print("voxel size: stored %.4f um, true %.4f um -> "
                      "rescaling lengths by %.4f" % (implied, voxel_size_um, f))
            for _, d in G.nodes(data=True):
                _rescale(d, f)
            for _, _, d in G.edges(data=True):
                _rescale(d, f)
        elif verbose:
            print("voxel size %.4f um already correct; no rescale" % implied)

    if verbose:
        print("%d pores, %d throats" % (G.number_of_nodes(),
                                        G.number_of_edges()))
    return G


def _rescale(attrs, f):
    for k, n in LENGTH_EXPONENT.items():
        if k in attrs:
            attrs[k] = np.asarray(attrs[k], dtype=float) * (f ** n)


def to_openpnm(G, reorder=None):
    """Rebuild an ``openpnm.network.Network`` from the graph. Node
    attributes become ``pore.*`` and edge attributes ``throat.*``.

    ``reorder`` converts the (z,y,x) fields to (x,y,z). Leave it as ``None``
    (the default) and the right thing happens: files stamped
    ``coord_order = "xyz"`` are already normalised and are left alone, older
    files are reordered. Only pass True/False to override that.
    """
    import openpnm as op

    if reorder is None:
        reorder = G.graph.get("coord_order", "mixed") != "xyz"

    nodes = sorted(G.nodes())
    index = {n: i for i, n in enumerate(nodes)}
    Np = len(nodes)
    edges = list(G.edges(data=True))
    conns = np.array([[index[a], index[b]] for a, b, _ in edges], dtype=int)

    pn = op.network.Network(conns=conns,
                            coords=np.array([G.nodes[n]["coords"]
                                             for n in nodes], dtype=float))

    def stack(dicts, key):
        vals = [np.asarray(d[key], dtype=float) for d in dicts]
        arr = np.array(vals)
        if reorder and key in ZYX_FIELDS and arr.ndim == 2 and \
                arr.shape[1] == 3:
            arr = arr[:, ::-1]
        return arr

    ndicts = [G.nodes[n] for n in nodes]
    for key in sorted({k for d in ndicts for k in d}):
        if key == "coords" or not all(key in d for d in ndicts):
            continue
        pn["pore." + key] = stack(ndicts, key)

    edicts = [d for _, _, d in edges]
    for key in sorted({k for d in edicts for k in d}):
        if key == "conns" or not all(key in d for d in edicts):
            continue
        pn["throat." + key] = stack(edicts, key)

    return pn


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("pickle_path")
    ap.add_argument("--voxel-um", type=float, default=None,
                    help="TRUE voxel size in microns (see network_names.csv)")
    ap.add_argument("--no-reorder", action="store_true",
                    help="keep global_peak etc. as stored (z,y,x)")
    ap.add_argument("--save-vtk", action="store_true",
                    help="also write <stem>.vtk for ParaView")
    args = ap.parse_args()

    G = load_graph(args.pickle_path, args.voxel_um)
    pn = to_openpnm(G, reorder=not args.no_reorder)
    print(pn)

    d = np.asarray(pn["throat.inscribed_diameter"], dtype=float) * 1e6
    print("\nthroat inscribed diameter: median %.2f um, range %.2f-%.2f um"
          % (np.median(d), d.min(), d.max()))
    ext = pn["pore.coords"].max(axis=0) - pn["pore.coords"].min(axis=0)
    print("domain extent: %.2f x %.2f x %.2f mm" % tuple(ext * 1e3))
    print("mean coordination Z = %.2f" % (2.0 * pn.Nt / pn.Np))

    if args.save_vtk:
        import openpnm as op
        out = args.pickle_path.rsplit(".", 1)[0]
        op.io.project_to_vtk(pn.project, filename=out)
        print("wrote:", out + ".vtk")


if __name__ == "__main__":
    main()
