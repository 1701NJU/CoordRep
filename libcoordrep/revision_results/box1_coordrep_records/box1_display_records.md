# Box 1. Representative complete CoordRep records

## A. Monodentate square-planar record: cis-[Pt(MeNH2)2I2]

```
CoordRep-v1
| M       = [Metal:Pt|ox:+2|d:d8|CN:4]
| Shape   = <ShapeBest:SP|Class:ideal|Delta:2|V:3.25,0.47>
| Stereo  = {trans:L1:N:1--L3:I:1}{trans:L2:N:1--L4:I:1}
| Ligands = |L1=[H]N([H])C([H])([H])[H]|L2=[H]N([H])C([H])([H])[H]|L3=I|L4=I|
| ID      = [L0_hash:ed29845b|L1:Pt|+2|SP/Td.ideal.D2|T2|I;I;[H]N([H])C([H])([H])[H];[H]N([H]|L2:Pt|CN4|SP|I;I;[H]N([H])C([H])([H])[H];[H]N([H])C([H])([H])[H|L3:Pt|CN4|I;I;[H]N([H])C([H])([H])[H];[H]N([H])C([H])([H])[H]]
```

Source: CSD IYIPOV  
Serializer: coordrep_v1

## B. Bidentate chelate record: [Co(en)3]3+

```
CoordRep-v1
| M       = [Metal:Co|ox:+3|d:d6|CN:6]
| Shape   = <ShapeBest:Oh|Class:dist|Delta:2|V:3.64,7.04>
| Stereo  = {trans:L1:N:1--L2:N:2}{trans:L1:N:2--L3:N:1}{trans:L2:N:1--L3:N:2}
| Ligands = |L1=C2D8N2|L2=C2D8N2|L3=C2D8N2|
| ID      = [L0_hash:71024466|L1:Co|+3|Oh/TPr.dist.D2|T3|C2D8N2;C2D8N2;C2D8N2|L2:Co|CN6|Oh|C2D8N2;C2D8N2;C2D8N2|L3:Co|CN6|C2D8N2;C2D8N2;C2D8N2]
```

Source: CSD JOXSOI  
Serializer: coordrep_v1

## C_fac. Octahedral fac stereochemical record: fac-[Co(dien)(CN)3]

```
CoordRep-v1
| M       = [Metal:Co|ox:+3|d:d6|CN:6]
| Shape   = <ShapeBest:Oh|Class:dist|Delta:1|V:5.76,6.98>
| Stereo  = {trans:L1:N:1--L2:C:1}{trans:L1:N:2--L3:C:1}{trans:L1:N:3--L4:C:1}{fm:fac}
| Ligands = |L1=[H]N([H])C([H])([H])C([H])([H])N([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]|L2=C#N|L3=C#N|L4=C#N|
| ID      = [L0_hash:83759f7d|L1:Co|+3|Oh/TPr.dist.D1|T3_fac|C#N;C#N;C#N;[H]N([H])C([H])([H])|L2:Co|CN6|Oh|C#N;C#N;C#N;[H]N([H])C([H])([H])C([H])([H])N([H])C|L3:Co|CN6|C#N;C#N;C#N;[H]N([H])C([H])([H])C([H])([H])N([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]]
```

Source: CSD IDARAG  
Serializer: coordrep_v1

---

## Supporting Information records

### C_mer. Octahedral mer stereochemical record (SI): mer-[Co(dien)(CN)3]

```
CoordRep-v1
| M       = [Metal:Co|ox:+3|d:d6|CN:6]
| Shape   = <ShapeBest:TPr|Class:dist|Delta:0|V:5.67,5.34>
| Stereo  = {trans:L1:N:1--L1:N:2}{trans:L1:N:3--L2:C:1}{trans:L3:C:1--L4:C:1}{fm:mer}
| Ligands = |L1=[H]N([H])C([H])([H])C([H])([H])N([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]|L2=C#N|L3=C#N|L4=C#N|
| ID      = [L0_hash:6a1a2c27|L1:Co|+3|TPr/Oh_boundary.dist.D0|T3_mer|C#N;C#N;C#N;[H]N([H])C(|L2:Co|CN6|TPr|C#N;C#N;C#N;[H]N([H])C([H])([H])C([H])([H])N([H])|L3:Co|CN6|C#N;C#N;C#N;[H]N([H])C([H])([H])C([H])([H])N([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]]
```

Source: CSD IDAREK  
Serializer: coordrep_v1

---

The records are direct outputs of the CoordRep serializer and are line-broken only for readability; full machine-readable records and hashes are provided in the Supporting Information.
