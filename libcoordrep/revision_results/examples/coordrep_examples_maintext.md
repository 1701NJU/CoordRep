# Full Annotated CoordRep Examples

These examples demonstrate the complete CoordRep representation
for real coordination complexes at varying complexity levels.

## A_monodentate_square_planar

**Description:** Monodentate-only square-planar Pt(II) CN=4, source: tmQM ABAMIA

**Metal:** Pt  |  **CN:** 4  |  **Best shape:** Td

**Full CoordRep string:**

```
[Metal:Pt|ox:+2|d:d8|CN:4]<ShapeBest:Td|Class:dist|Delta:2|V:3.14,7.62>{trans:L1:P:1--L3:Cl:1}{trans:L2:P:1--L4:Cl:1}|L1=P|L2=P|L3=Cl|L4=Cl|
```

**Pretty-printed (multiline):**

```
[Metal:Pt|ox:+2|d:d8|CN:4]
<ShapeBest:Td|Class:dist|Delta:2|V:3.14,7.62>
{trans:L1:P:1--L3:Cl:1}
{trans:L2:P:1--L4:Cl:1}
|L1=P|
|L2=P|
|L3=Cl|
|L4=Cl|
```

**Annotation table:**

| Field | Value |
|-------|-------|
| Metal block | `[Metal:Pt|ox:+2|d:d8|CN:4]` |
| Shape block | `<ShapeBest:Td|Class:dist|Delta:2|V:3.14,7.62>` |
| Stereo constraint 1 | `{trans:L1:P:1--L3:Cl:1}` |
| Stereo constraint 2 | `{trans:L2:P:1--L4:Cl:1}` |
| Ligand L1 | `P` |
| Ligand L2 | `P` |
| Ligand L3 | `Cl` |
| Ligand L4 | `Cl` |
| Donor sites | L1:P:1, L3:Cl:1, L2:P:1, L4:Cl:1 |
| Valid | True |

**Identity keys:**

- **L0 (StateKey):** `[Metal:Pt|ox:+2|d:d8|CN:4]<ShapeBest:Td|Class:dist|Delta:2|V:3.14,7.62>{trans:L1…`
- **L1 (ShapeID):** `Pt|+2|Td/SP.dist.D2|T2|Cl;Cl;P;P…`
- **L2 (TopoID):** `Pt|CN4|Td|Cl;Cl;P;P…`
- **L3 (ConnID):** `Pt|CN4|Cl;Cl;P;P…`

---

## B_bidentate_chelating

**Description:** Bidentate chelating Pd CN=4, dents=[2, 2], source: tmQM ABAFOZ

**Metal:** Pd  |  **CN:** 4  |  **Best shape:** Td

**Full CoordRep string:**

```
[Metal:Pd|ox:+2|d:d8|CN:4]<ShapeBest:Td|Class:good|Delta:3|V:2.82,8.22>{trans:L1:N:1--L2:O:1}{trans:L1:N:2--L2:O:2}|L1=C18H24N2|L2=C22H16O2|
```

**Pretty-printed (multiline):**

```
[Metal:Pd|ox:+2|d:d8|CN:4]
<ShapeBest:Td|Class:good|Delta:3|V:2.82,8.22>
{trans:L1:N:1--L2:O:1}
{trans:L1:N:2--L2:O:2}
|L1=C18H24N2|
|L2=C22H16O2|
```

**Annotation table:**

| Field | Value |
|-------|-------|
| Metal block | `[Metal:Pd|ox:+2|d:d8|CN:4]` |
| Shape block | `<ShapeBest:Td|Class:good|Delta:3|V:2.82,8.22>` |
| Stereo constraint 1 | `{trans:L1:N:1--L2:O:1}` |
| Stereo constraint 2 | `{trans:L1:N:2--L2:O:2}` |
| Ligand L1 | `C18H24N2` |
| Ligand L2 | `C22H16O2` |
| Donor sites | L1:N:1, L2:O:1, L1:N:2, L2:O:2 |
| Valid | True |

**Identity keys:**

- **L0 (StateKey):** `[Metal:Pd|ox:+2|d:d8|CN:4]<ShapeBest:Td|Class:good|Delta:3|V:2.82,8.22>{trans:L1…`
- **L1 (ShapeID):** `Pd|+2|Td/SP.good.D3|T2|C18H24N2;C22H16O2…`
- **L2 (TopoID):** `Pd|CN4|Td|C18H24N2;C22H16O2…`
- **L3 (ConnID):** `Pd|CN4|C18H24N2;C22H16O2…`

---

## C_octahedral_multi_stereo

**Description:** Octahedral Ir CN=6, dents=[2, 2, 2], multiple trans constraints, source: tmQM ABAYUA

**Metal:** Ir  |  **CN:** 6  |  **Best shape:** Oh

**Full CoordRep string:**

```
[Metal:Ir|CN:6]<ShapeBest:Oh|Class:dist|Delta:2|V:4.37,7.88>{trans:L1:C:1--L3:N:1}{trans:L1:N:1--L2:N:1}{trans:L2:C:1--L3:O:1}|L1=C11F2H6N|L2=C11F2H6N|L3=C15H14NO|
```

**Pretty-printed (multiline):**

```
[Metal:Ir|CN:6]
<ShapeBest:Oh|Class:dist|Delta:2|V:4.37,7.88>
{trans:L1:C:1--L3:N:1}
{trans:L1:N:1--L2:N:1}
{trans:L2:C:1--L3:O:1}
|L1=C11F2H6N|
|L2=C11F2H6N|
|L3=C15H14NO|
```

**Annotation table:**

| Field | Value |
|-------|-------|
| Metal block | `[Metal:Ir|CN:6]` |
| Shape block | `<ShapeBest:Oh|Class:dist|Delta:2|V:4.37,7.88>` |
| Stereo constraint 1 | `{trans:L1:C:1--L3:N:1}` |
| Stereo constraint 2 | `{trans:L1:N:1--L2:N:1}` |
| Stereo constraint 3 | `{trans:L2:C:1--L3:O:1}` |
| Ligand L1 | `C11F2H6N` |
| Ligand L2 | `C11F2H6N` |
| Ligand L3 | `C15H14NO` |
| Donor sites | L1:C:1, L3:N:1, L1:N:1, L2:N:1, L2:C:1, L3:O:1 |
| Valid | True |

**Identity keys:**

- **L0 (StateKey):** `[Metal:Ir|CN:6]<ShapeBest:Oh|Class:dist|Delta:2|V:4.37,7.88>{trans:L1:C:1--L3:N:…`
- **L1 (ShapeID):** `Ir|?|Oh/TPr.dist.D2|T3|C11F2H6N;C11F2H6N;C15H14NO…`
- **L2 (TopoID):** `Ir|CN6|Oh|C11F2H6N;C11F2H6N;C15H14NO…`
- **L3 (ConnID):** `Ir|CN6|C11F2H6N;C11F2H6N;C15H14NO…`

---

