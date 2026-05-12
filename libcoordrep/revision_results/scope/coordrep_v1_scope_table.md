# CoordRep v1 Scope Table

## Supported in CoordRep v1

| Feature | Details |
|---------|--------|
| Mononuclear transition-metal complexes | 3d/4d/5d metals, single metal center |
| Atom-resolved eta-1 donors | Each donor atom tracked by element and rank |
| Monodentate ligands | Single-donor ligands (e.g., Cl, NH3, H2O) |
| Bidentate eta-1 chelating ligands | e.g., en, bpy, acac, oxalate |
| Multidentate eta-1 chelating ligands | dent >= 3 (e.g., EDTA, tpy, porphyrin) |
| Coordination number 2-14 | Full CN range from linear to high-coordinate |
| cis/trans, fac/mer, donor-relation tokens | Stereochemical constraint grammar |
| CoordRep-ID L0-L3 | Hierarchical identity: StateKey, ShapeID, TopoID, ConnID |

## Explicitly Out of v1 Scope

| Exclusion | Reason |
|-----------|--------|
| Multinuclear complexes | Multiple metal centers, bridging ligands |
| MOFs / coordination polymers | Periodic / extended lattice structures |
| Polyoxometalates | Large cluster metal-oxide frameworks |
| eta-n organometallics | Metallocenes, allyl, arene complexes (eta > 1) |
| Positional disorder / partial occupancy | Crystallographic disorder |
| Ambiguous donor assignment | Weak interactions below distance threshold |

## Future Extension Route

| Extension | Description |
|-----------|------------|
| CoordRep-Multi | Multi-metal graph + bridging ligand relations |
| CoordRep-Periodic | Metal node / linker / periodic topology |
| CoordRep-Haptic | eta-n fragment token + centroid / slippage geometry |
| CoordRep-Occ | Occupancy-aware disorder representation |

---
*CoordRep v1 focuses on mononuclear eta-1 coordination snapshots, including monodentate and multidentate chelating ligands.*
