# CoordRep-Haptic-v0: Prototype π/Haptic Records
**Scope:** Feasibility demonstration only. CoordRep v1 remains restricted to mononuclear, atom-resolved η1 coordination snapshots.
**Cases:** 8 curated structures (η2–η6 coverage)

## haptic_001: Ferrocene, bis(η5-cyclopentadienyl)iron(II)
- **Refcode:** FEROCE01
- **Metal:** Fe
- **Haptic class:** metallocene (η5)
- **Site:** pi_fragment, 5 atoms
- **Notes:** Canonical metallocene — both rings η5

```
[Metal: Fe|ox:+2|d:d6]
[Sites: S1={type:pi_fragment|ligand:L1|atoms:C1,C2,C3,C4,C5|eta:5|centroid:Cp1}; S2={type:pi_fragment|ligand:L2|atoms:C6,C7,C8,C9,C10|eta:5|centroid:Cp2}]
[Hapticity: {S1->M1|eta:5|mode:Cp}; {S2->M1|eta:5|mode:Cp}]
[Ligands: L1=C1=CC=CC1|L2=C1=CC=CC1]
[ID: HapticState=Fe|d6|eta5-Cp,eta5-Cp]
```

## haptic_002: Potassium trichloro(η2-ethylene)platinate(II) (Zeise's salt)
- **Refcode:** ZEISSE01
- **Metal:** Pt
- **Haptic class:** eta2_alkene (η2)
- **Site:** pi_fragment, 2 atoms
- **Notes:** Zeise's salt — prototypical η2-alkene coordination

```
[Metal: Pt|ox:+2|d:d8|CN:4]
[Sites: S1={type:pi_fragment|ligand:L1|atoms:C1,C2|eta:2|centroid:C=C_mid}]
[Hapticity: {S1->M1|eta:2|mode:alkene}]
[AtomDonors: {L2:Cl:1->M1}; {L3:Cl:1->M1}; {L4:Cl:1->M1}]
[Ligands: L1=C=C|L2=[Cl-]|L3=[Cl-]|L4=[Cl-]]
[ID: HapticState=Pt|d8|eta2-C2H4,Cl,Cl,Cl]
```

## haptic_003: Dichloro(η3-allyl)palladium(II) dimer
- **Refcode:** PDALLY
- **Metal:** Pd
- **Haptic class:** eta3_allyl (η3)
- **Site:** pi_fragment, 3 atoms
- **Notes:** η3-allyl complex — 3-carbon π fragment

```
[Metal: Pd|ox:+2|d:d8]
[Sites: S1={type:pi_fragment|ligand:L1|atoms:C1,C2,C3|eta:3|centroid:allyl_mid}]
[Hapticity: {S1->M1|eta:3|mode:allyl}]
[AtomDonors: {L2:Cl:1->M1}; {L3:Cl:1->M1}]
[Ligands: L1=C=CC|L2=[Cl-]|L3=[Cl-]]
[ID: HapticState=Pd|d8|eta3-allyl,Cl,Cl]
```

## haptic_004: Tricarbonyl(η6-benzene)chromium(0)
- **Refcode:** BZCRCB01
- **Metal:** Cr
- **Haptic class:** eta6_arene (η6)
- **Site:** pi_fragment, 6 atoms
- **Notes:** Piano-stool η6-arene complex

```
[Metal: Cr|ox:0|d:d6]
[Sites: S1={type:pi_fragment|ligand:L1|atoms:C1,C2,C3,C4,C5,C6|eta:6|centroid:Bz_mid}]
[Hapticity: {S1->M1|eta:6|mode:arene}]
[AtomDonors: {L2:C:1->M1}; {L3:C:1->M1}; {L4:C:1->M1}]
[Ligands: L1=c1ccccc1|L2=[C-]#[O+]|L3=[C-]#[O+]|L4=[C-]#[O+]]
[ID: HapticState=Cr|d6|eta6-benzene,CO,CO,CO]
```

## haptic_005: (η5-Pentamethylcyclopentadienyl)(η2:η2-cycloocta-1,5-diene)chlororuthenium(II)
- **Refcode:** RUCPCL
- **Metal:** Ru
- **Haptic class:** mixed_haptic (η5)
- **Site:** pi_fragment, 5 atoms
- **Notes:** Mixed hapticity: η5-Cp* + two η2 from cod + σ-Cl

```
[Metal: Ru|ox:+2|d:d6]
[Sites: S1={type:pi_fragment|ligand:L1|atoms:C1,C2,C3,C4,C5|eta:5|centroid:Cp*1}; S2={type:pi_fragment|ligand:L2|atoms:C11,C12|eta:2|centroid:cod_a}; S3={type:pi_fragment|ligand:L2|atoms:C15,C16|eta:2|centroid:cod_b}]
[Hapticity: {S1->M1|eta:5|mode:Cp}; {S2->M1|eta:2|mode:alkene}; {S3->M1|eta:2|mode:alkene}]
[AtomDonors: {L3:Cl:1->M1}]
[Ligands: L1=CC1=C(C)C(C)=C1C|L2=C1CC=CCCC=C1|L3=[Cl-]]
[ID: HapticState=Ru|d6|eta5-Cp*,eta2-cod,eta2-cod,Cl]
```

## haptic_006: Bis(triphenylphosphine)(η2-ethylene)nickel(0)
- **Refcode:** NIEPPH
- **Metal:** Ni
- **Haptic class:** eta2_alkene (η2)
- **Site:** pi_fragment, 2 atoms
- **Notes:** Trigonal planar Ni(0) with one η2-alkene

```
[Metal: Ni|ox:0|d:d10]
[Sites: S1={type:pi_fragment|ligand:L1|atoms:C1,C2|eta:2|centroid:C=C_mid}]
[Hapticity: {S1->M1|eta:2|mode:alkene}]
[AtomDonors: {L2:P:1->M1}; {L3:P:1->M1}]
[Ligands: L1=C=C|L2=c1ccc(cc1)P(c1ccccc1)c1ccccc1|L3=c1ccc(cc1)P(c1ccccc1)c1ccccc1]
[ID: HapticState=Ni|d10|eta2-C2H4,PPh3,PPh3]
```

## haptic_007: Tricarbonyl(η5-cyclopentadienyl)manganese(I)
- **Refcode:** CPMNCO
- **Metal:** Mn
- **Haptic class:** half_sandwich (η5)
- **Site:** pi_fragment, 5 atoms
- **Notes:** Half-sandwich CpMn(CO)3 — classic piano stool

```
[Metal: Mn|ox:+1|d:d6]
[Sites: S1={type:pi_fragment|ligand:L1|atoms:C1,C2,C3,C4,C5|eta:5|centroid:Cp1}]
[Hapticity: {S1->M1|eta:5|mode:Cp}]
[AtomDonors: {L2:C:1->M1}; {L3:C:1->M1}; {L4:C:1->M1}]
[Ligands: L1=C1=CC=CC1|L2=[C-]#[O+]|L3=[C-]#[O+]|L4=[C-]#[O+]]
[ID: HapticState=Mn|d6|eta5-Cp,CO,CO,CO]
```

## haptic_008: Dichloro-bis(η5-cyclopentadienyl)titanium(IV)
- **Refcode:** TICPDC01
- **Metal:** Ti
- **Haptic class:** bent_metallocene (η5)
- **Site:** pi_fragment, 5 atoms
- **Notes:** Bent metallocene Cp2TiCl2 — Ziegler-Natta catalyst precursor

```
[Metal: Ti|ox:+4|d:d0]
[Sites: S1={type:pi_fragment|ligand:L1|atoms:C1,C2,C3,C4,C5|eta:5|centroid:Cp1}; S2={type:pi_fragment|ligand:L2|atoms:C6,C7,C8,C9,C10|eta:5|centroid:Cp2}]
[Hapticity: {S1->M1|eta:5|mode:Cp}; {S2->M1|eta:5|mode:Cp}]
[AtomDonors: {L3:Cl:1->M1}; {L4:Cl:1->M1}]
[Ligands: L1=C1=CC=CC1|L2=C1=CC=CC1|L3=[Cl-]|L4=[Cl-]]
[ID: HapticState=Ti|d0|eta5-Cp,eta5-Cp,Cl,Cl]
```

