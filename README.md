# Heterogeneous pore networks

Pore-network models extracted from micro-CT images of 15 porous media,
spanning bead packs and sand packs, synthetic micromodels, sandstones, and
carbonates. Where possible each sample is provided under three independent
extraction methods, so that method-dependence can be separated from rock
properties.

**The networks themselves are hosted on the Digital Porous Media Portal**
(see [Data](#data)). This repository holds the loader/converter code and the
metadata needed to use them correctly.

## Contents

| file | purpose |
| --- | --- |
| `network_to_openpnm.py` | load a network pickle, correct its voxel size, convert to OpenPNM |
| `network_names.csv` | maps each file to its publication name, extraction method and **true voxel size** |
| `export_networks.py` | renames the original extraction outputs to the published names (provenance/reproducibility) |

## Quick start

```python
from network_to_openpnm import load_graph, to_openpnm

G  = load_graph("clashach_sandstone_diamorse.pickle", voxel_size_um=7.0)
pn = to_openpnm(G)
print(pn)
```

or from the command line:

```bash
python network_to_openpnm.py clashach_sandstone_diamorse.pickle --voxel-um 7.0
```

Requires `numpy`, `networkx`, and `openpnm >= 3.0` (the latter only for the
OpenPNM conversion; `load_graph` works without it).

## Format

Each file is a pickled NetworkX `Graph`: one node per pore, one edge per
throat. Node and edge attributes use OpenPNM property names with the
`pore.`/`throat.` prefixes stripped — `to_openpnm` puts them back. Pores
carry `coords`, `volume`, `inscribed_diameter`, `equivalent_diameter`,
`surface_area`, `region_label` and others; throats carry `conns`,
`inscribed_diameter`, `equivalent_diameter`, `length`, `total_length`,
`cross_sectional_area`, `hydraulic_size_factors`, `diffusive_size_factors`
and others.

## Four things to know before using these files

**1. Units are metres.** Every length, area and volume attribute is SI —
not voxels, not microns. Multiply by `1e6` for microns.

**2. Coordinate order is not uniform.** Pore `coords` is `(x, y, z)`.
But `global_peak`, `geometric_centroid` and `local_peak` are `(z, y, x)`,
following the porespy region convention. `to_openpnm` reverses those fields
so that everything comes out `(x, y, z)`; pass `reorder=False` to keep them
as stored.

**3. Some networks were extracted with a placeholder voxel size.** In those
files, absolute lengths are off by a constant factor — the topology is
correct but permeability and capillary pressures computed from them would
not be. `network_names.csv` flags these (`rescale_needed = YES`) and gives
the correct voxel size. Passing `voxel_size_um` to `load_graph` applies the
correction: lengths scale by *f*, areas by *f²*, volumes by *f³*, and the
OpenPNM `*_size_factors` by *f* (they have units of length, being
area/length). **Always pass the voxel size from `network_names.csv`** — it
is a no-op for files that are already correct, so it is safe to pass every
time.

**4. Some files contain leftover solver nodes.** A few networks were saved
after a max-flow calculation and still contain `super_source`/`super_sink`
nodes carrying no geometry. `load_graph` removes any non-integer node and
reindexes the remainder contiguously.

## Extraction methods

| method | family | notes |
| --- | --- | --- |
| Diamorse | discrete Morse theory | topology from the signed Euclidean distance transform |
| SNOW2 | marker-based watershed | `porespy.networks.snow2` |
| MAGNET | medial axis / skeleton | `porespy.networks.magnet` |

These are **not** interchangeable. Pore and throat counts differ by up to an
order of magnitude for the same image, and mean coordination *Z* = 2*N*t/*N*p
differs systematically: Diamorse gives 3.8–8.1, SNOW2 3.9–6.9, and MAGNET
2.1–2.8. MAGNET's *Z* ≈ 2 reflects its chain topology — it inserts pores
along skeleton branches, so its "throats" are skeleton segments rather than
pore necks. It is not reporting a more poorly connected rock. Compare
counts and coordination *within* an extraction method, not across.

Quantities that do transfer across all three methods include the critical
percolation (breakthrough) diameter and the location of the network's
minimum cut.

## Data

The network files are available from the Digital Porous Media Portal:

> *(add the portal DOI / dataset link here once the upload is complete)*

## Citation

> *(add citation once published)*

## License

> *(choose a license — CC-BY-4.0 is common for data, MIT for the code)*
