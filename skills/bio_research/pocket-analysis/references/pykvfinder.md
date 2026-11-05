# pyKVFinder Reference

## What it does

pyKVFinder is a Python package for detecting protein cavities by scanning a 3D
grid around and within the protein. It provides a Python API as an alternative
to invoking the external fpocket executable.

## Installation

```bash
pip install pyKVFinder
```

Python dependencies include numpy and toml. Provision pyKVFinder in the
registry-selected scientific environment; it is separate from the fpocket CLI.

## Basic Usage

The example follows the pyKVFinder 0.9.x API used by
[`scripts/pykvfinder_utils.py`](../scripts/pykvfinder_utils.py):

```python
import pyKVFinder

# 1. Read the PDB
atomic = pyKVFinder.read_pdb("protein.pdb")

# 2. Obtain the grid vertices
vertices = pyKVFinder.get_vertices(atomic, probe_out=4.0, step=0.6)

# 3. Detect cavities
ncavities, cavities = pyKVFinder.detect(
    atomic,
    vertices,
    step=0.6,          # Grid spacing (Å): smaller gives finer resolution but is slower
    probe_in=1.4,      # Inner probe radius: defines the protein–cavity boundary
    probe_out=4.0,     # Outer probe radius: defines the cavity–solvent boundary
    removal_distance=2.4,  # Remove cavity regions too close to solvent
    volume_cutoff=5.0,     # Minimum cavity volume (Å³)
)

print(f"Found {ncavities} cavities")

# 4. Calculate surface, volume, and area
surface, volumes, areas = pyKVFinder.spatial(cavities, step=0.6)
# volumes = {"KAA": 125.3, "KAB": 89.7, ...}

# 5. Obtain cavity residues
residues = pyKVFinder.constitutional(cavities, atomic, vertices, step=0.6)
# residues = {"KAA": [["83", "A", "LEU"], ["81", "A", "GLU"], ...], ...}

# 6. Export the result
pyKVFinder.export("cavities.pdb", cavities, surface, vertices, step=0.6)
```

## Key Parameters

| Parameter | Default | Description | Comparison with fpocket |
|------|--------|------|----------------|
| `step` | 0.6 Å | Grid resolution; smaller values give finer resolution | No corresponding grid parameter |
| `probe_in` | 1.4 Å | Inner probe defining the protein surface | Analogous to fpocket `-m` (minimum sphere radius) |
| `probe_out` | 4.0 Å | Outer probe defining the solvent-accessible boundary | Analogous to fpocket `-M` (maximum sphere radius) |
| `removal_distance` | 2.4 Å | Removes shallow surface depressions | No direct equivalent |
| `volume_cutoff` | 5.0 Å³ | Minimum cavity volume | Analogous filtering role to fpocket `-i` (minimum sphere count), with different units |

## Output

- `ncavities`: integer count of detected cavities.
- `cavities`: 3D NumPy array of integer voxel labels; cavity integers map to
  named IDs such as "KAA" and "KAB" in exported records.
- `volumes`: dictionary mapping cavity IDs to volumes in Å³.
- `residues`: dictionary mapping cavity IDs to lists of `[resnum, chain, resname]`.

## Comparison with fpocket

| Feature | pyKVFinder | fpocket |
|------|-----------|---------|
| Installation | `pip install` | Conda or compile C source |
| Interface | Python API + NumPy arrays | C executable |
| Method | Grid scanning | Voronoi + alpha spheres |
| Illustrative runtime | Moderate (~5s / protein) | Fast (<1s) |
| Druggability score | No | Yes |
| Volume calculation | Grid counting | Monte Carlo estimate |
| Residue assignment | Yes | Yes |
| Visualization output | PDB format | PDB + PQR |

The runtime figures are illustrative, not a benchmark guarantee; actual cost
depends on protein size, grid resolution, and hardware. Grid-counted volumes
depend on the selected resolution.

## Pitfalls

- **Large step (>1.0)**: may miss small cavities or underestimate their volume.
- **Small probe_out (<3.0)**: may fail to detect shallow surface pockets.
- **Large proteins (>5000 atoms)**: grid construction can use substantial memory;
  `step=0.8` can reduce this cost when used through an interface that exposes it.
  The current pipeline adapter fixes `step=0.6`; do not invent a config override.
- **No druggability score**: pyKVFinder does not assess whether a pocket is suitable
  for drug binding; use an additional method for that assessment.

## References

- Guerra et al., *BMC Bioinformatics* 22:607 (2021) — pyKVFinder
- Guerra et al., *SoftwareX* 12:100606 (2020) — KVFinder-web
