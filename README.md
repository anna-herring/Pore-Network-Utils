# Pore-Network-Utils

Utilities for working with pore-network models extracted from micro-CT
images of porous geologic media — 15 samples spanning bead and sand packs,
a synthetic micromodel, sandstones and carbonates. Where possible each
sample is provided under three independent extraction methods, so that
method-dependence can be separated from rock properties.

**The networks themselves are on the Digital Porous Media Portal** (see
[Data](#data)); this repository holds the code to load them and the
conventions you need to use them correctly.

> **Licence pending.** No licence has been chosen yet, which means the code
> is under default copyright and cannot yet be reused. See
> [Licence](#licence).

## Contents

Distributed with the networks:

| file | purpose |
| --- | --- |
| `network_to_openpnm.py` | load a network pickle and convert it to OpenPNM |
| `network_names.csv` | maps each file to its publication name, extraction method and voxel size |

In this repository only, documenting how the published files were produced:

| file | purpose |
| --- | --- |
| `correct_networks.py` | the one-off pass that produced the v1.0 files from the raw extraction output |
| `export_networks.py` | maps the original extraction filenames to the published names |
| `make_checksums.py` | writes / verifies `CHECKSUMS.txt` for the deposit |
| `correction_report.csv` | per-file record of what `correct_networks.py` changed, and its verification results |

Verify a download against the published checksums with `sha256sum -c
CHECKSUMS.txt`, or `python make_checksums.py --dir . --check`.

## Quick start

```python
from network_to_openpnm import load_graph, to_openpnm

G  = load_graph("clashach_sandstone_diamorse.pickle")
pn = to_openpnm(G)
print(pn)
```

or from the command line:

```bash
python network_to_openpnm.py clashach_sandstone_diamorse.pickle
```

No voxel size or coordinate handling is needed — see [Conventions](#conventions).

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

## Conventions

Every file carries a stamp describing itself:

```python
G.graph
# {'format_version': '1.0',
#  'coord_order': 'xyz',
#  'attribute_types': 'plain Python (float / int / list)',
#  'units': 'SI (metres)',
#  'voxel_size_um': 7.0,
#  'capacity_definition': '(inscribed_diameter/2)**4, m^4',
#  'sample_name': 'Clashach Sandstone',
#  'extraction': 'Diamorse'}
```

**Units are metres.** Every length, area and volume attribute is SI — not
voxels, not microns. Multiply by `1e6` for microns.

**Every vector attribute is `(x, y, z)`.** That includes `global_peak`,
`geometric_centroid` and `local_peak`, which the extraction tools emit as
`(z, y, x)`. They were reordered when these files were built, so no
transposition is needed on load.

**Voxel sizes are already correct.** A few extractions recorded a
placeholder voxel size, leaving absolute lengths off by a constant factor.
That was corrected before publication — `voxel_size_um` in the stamp is the
true value, and it has already been applied. Topology, permeability and
capillary pressures are all on the correct length scale as distributed.

**Throat `capacity` is `(inscribed_diameter/2)⁴`, in m⁴.** This is the
definition behind the published minimum cuts: nonwetting-phase invasion is
controlled by the constriction, so the inscribed radius is the relevant
one. If you want the equivalent-area version instead, it is one line —
`equivalent_diameter` is stored on every throat. The two give substantially
different minimum cuts, so be explicit about which you are using.

**No solver nodes.** Nothing but integer pore ids appears in the node set.

**Attributes are plain Python types** — `float`, `int`, and `list` — not
numpy scalars or arrays. Loading therefore needs only `networkx`, with no
dependence on numpy's array ABI, so these files are not tied to the numpy
version they were written with. Wrap anything you want to compute on in
`np.asarray()` as usual; `to_openpnm` already does.

Earlier, uncorrected copies of these networks exist in the authors' working
directories and carry none of the above guarantees. Files with
`format_version` absent are pre-1.0 and should not be used; in particular
their `capacity` attribute is built on the equivalent diameter and will not
reproduce the published cuts.

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

## Licence

Not yet chosen. Until a `LICENSE` file is added this code is under default
copyright — publicly visible, but not licensed for reuse, modification or
redistribution.

Intended: **MIT** for the code here, **CC-BY-4.0** for the network data on
the portal. Pending a check of the DE-SC0025400 award terms and a software
disclosure to the UT Research Foundation, as the work was produced under a
federal award. Once settled, add the `LICENSE` file and uncomment the
`license:` field in `CITATION.cff`.
