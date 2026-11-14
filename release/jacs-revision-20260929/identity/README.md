# Constructed L1 stereochemistry and chelate controls

`rc3_constructed_collision_controls.py` executes the published 1.1.2rc3
shape, relation, exact canonicalization, serialization, and identity-key
path from the adjacent frozen release. It constructs idealized
cis/trans-[Pt(NH3)2Cl2], fac/mer-[Co(NH3)3Cl3], and octahedral
[Co(glycinato)(NH3)4]2+ first spheres. The generated JSON and CSV are
included so that the exact keys and field differences can be inspected
without executing code.

The Pt pair and the Co pair each give **two L0 strings and one L1 key**;
L2 and L3 also coincide within each pair. The common Pt key is:

```text
Pt|+2|CN4|SP/Td.ideal.D3|T2_C4|SMILES:[Cl-:1];SMILES:[Cl-:1];SMILES:[NH3:1];SMILES:[NH3:1]
```

The common Co key is:

```text
Co|+3|CN6|Oh/TPr.ideal.D3|T3_C12|SMILES:[Cl-:1];SMILES:[Cl-:1];SMILES:[Cl-:1];SMILES:[NH3:1];SMILES:[NH3:1];SMILES:[NH3:1]
```

The complete L0 relation maps distinguish the pairs without adding manual
cis/trans or fac/mer labels. A distinct supported explicit fac/mer marker
can survive L1, but the current detector emits that marker only for its
defined tridentate-ligand pattern; it does not emit one for the classic
monodentate Co pair above. Thus an ordinary fac/mer distinction is **not**
generally guaranteed at L1. In the chelate control, one glycinato ligand
has two pointers, `L1:N:1` and `L1:O:1`. L0 retains their same-ligand cis
relation; L1 retains the ligand payload and aggregate relation count but
not either donor pointer.

These constructed first spheres are not full CSD-ingestion tests or
estimates of stereoisomer prevalence. Their nominal ligand SMILES and
attachment descriptors are explicitly supplied to test the serializer and
key extractor under controlled conditions.

From the repository root, with Python, NumPy, NetworkX, and RDKit available:

```powershell
python -B release/jacs-revision-20260929/identity/rc3_constructed_collision_controls.py
```

The script checks that its constructed pairs share their L1–L3 keys, differ
at L0, preserve the expected element-resolved trans relations, and remain
invariant when ligand traversal is reversed. It overwrites only its own
generated JSON and CSV in this directory.
