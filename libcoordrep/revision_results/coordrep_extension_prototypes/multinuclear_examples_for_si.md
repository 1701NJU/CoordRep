# CoordRep-Multi-v0: Prototype Multinuclear Records
**Scope:** Feasibility demonstration only. CoordRep v1 remains restricted to mononuclear, atom-resolved η1 coordination snapshots.
**Cases:** 15 curated structures (metal_count = 2–4, bridge types: μ2/μ3)

## multi_001: Di-μ2-chlorido-bis[(2,2'-bipyridine)copper(II)]
- **Refcode:** BIHWOL
- **Metals:** Cu,Cu (×2)
- **Bridges:** mu2-Cl (×2)
- **Local geometries:** CN = 5,5, shapes = SPY,SPY
- **Notes:** Classic di-μ2-Cl bridged dinuclear Cu(II)

```
[Metals: M1=Cu; M2=Cu]
[MetalGraph: {M1--M2|relation:bridged|dMM:3.42}]
[LocalSphere: M1=<ShapeBest:SPY|CN:5>{L1:N:1,L1:N:2,L3:Cl:1,L4:Cl:1,L5:Cl:1}; M2=<ShapeBest:SPY|CN:5>{L2:N:1,L2:N:2,L3:Cl:1,L4:Cl:1,L6:Cl:1}]
[Bridges: {L3:Cl->M1,M2|mu:2}; {L4:Cl->M1,M2|mu:2}]
[Ligands: L1=c1ccnc(c1)-c1ccccn1|L2=c1ccnc(c1)-c1ccccn1|L3=[Cl-]|L4=[Cl-]|L5=[Cl-]|L6=[Cl-]]
[ID: LocalID[M1]=Cu|SPY|5|bipy,Cl,Cl,Cl; LocalID[M2]=Cu|SPY|5|bipy,Cl,Cl,Cl; GlobalID=MULTI_001]
```

## multi_002: μ-Oxo-bis[trichloroiron(III)]
- **Refcode:** YOMHEW
- **Metals:** Fe,Fe (×2)
- **Bridges:** mu2-O (×1)
- **Local geometries:** CN = 4,4, shapes = Td,Td
- **Notes:** Fe2O core with Td local geometry on each Fe

```
[Metals: M1=Fe; M2=Fe]
[MetalGraph: {M1--M2|relation:bridged|dMM:3.56}]
[LocalSphere: M1=<ShapeBest:Td|CN:4>{L3:O:1,L4:Cl:1,L5:Cl:1,L6:Cl:1}; M2=<ShapeBest:Td|CN:4>{L3:O:1,L7:Cl:1,L8:Cl:1,L9:Cl:1}]
[Bridges: {L3:O->M1,M2|mu:2}]
[Ligands: L3=[O-2]|L4=[Cl-]|L5=[Cl-]|L6=[Cl-]|L7=[Cl-]|L8=[Cl-]|L9=[Cl-]]
[ID: LocalID[M1]=Fe|Td|4|O,Cl,Cl,Cl; LocalID[M2]=Fe|Td|4|O,Cl,Cl,Cl; GlobalID=MULTI_002]
```

## multi_003: Tetrakis(μ-acetato)dirhodium(II)
- **Refcode:** RHACAC03
- **Metals:** Rh,Rh (×2)
- **Bridges:** mu2-OAc (×4)
- **Local geometries:** CN = 5,5, shapes = SPY,SPY
- **Notes:** Rh2(OAc)4 paddlewheel with M-M single bond

```
[Metals: M1=Rh; M2=Rh]
[MetalGraph: {M1--M2|relation:bonded|dMM:2.39}]
[LocalSphere: M1=<ShapeBest:SPY|CN:5>{L1:O:1,L2:O:1,L3:O:1,L4:O:1,M2}; M2=<ShapeBest:SPY|CN:5>{L1:O:2,L2:O:2,L3:O:2,L4:O:2,M1}]
[Bridges: {L1:O->M1,M2|mu:2}; {L2:O->M1,M2|mu:2}; {L3:O->M1,M2|mu:2}; {L4:O->M1,M2|mu:2}]
[Ligands: L1=CC(=O)[O-]|L2=CC(=O)[O-]|L3=CC(=O)[O-]|L4=CC(=O)[O-]]
[ID: LocalID[M1]=Rh|SPY|5|O4+Rh; LocalID[M2]=Rh|SPY|5|O4+Rh; GlobalID=MULTI_003]
```

## multi_004: μ3-Oxo-tris(μ-acetato)trichromium(III)
- **Refcode:** TRCRAC
- **Metals:** Cr,Cr,Cr (×3)
- **Bridges:** mu3-O;mu2-OAc (×4)
- **Local geometries:** CN = 6,6,6, shapes = Oh,Oh,Oh
- **Notes:** Basic chromium(III) acetate — triangular μ3-O core

```
[Metals: M1=Cr; M2=Cr; M3=Cr]
[MetalGraph: {M1--M2|relation:bridged|dMM:3.28}; {M1--M3|relation:bridged|dMM:3.30}; {M2--M3|relation:bridged|dMM:3.29}]
[LocalSphere: M1=<ShapeBest:Oh|CN:6>{L7:O:1,L1:O:1,L2:O:1,L5:O:1,L6:O:1,L8:O:1}; M2=<ShapeBest:Oh|CN:6>{L7:O:1,L1:O:2,L3:O:1,L5:O:2,L6:O:2,L9:O:1}; M3=<ShapeBest:Oh|CN:6>{L7:O:1,L2:O:2,L4:O:1,L3:O:2,L4:O:2,L10:O:1}]
[Bridges: {L7:O->M1,M2,M3|mu:3}; {L1:O->M1,M2|mu:2}; {L2:O->M1,M3|mu:2}; {L3:O->M2,M3|mu:2}]
[Ligands: L1=CC(=O)[O-]|L2=CC(=O)[O-]|L3=CC(=O)[O-]|L5=CC(=O)[O-]|L6=CC(=O)[O-]|L7=[O-2]|L8=O|L9=O|L10=O]
[ID: LocalID[M1]=Cr|Oh|6|O6; LocalID[M2]=Cr|Oh|6|O6; LocalID[M3]=Cr|Oh|6|O6; GlobalID=MULTI_004]
```

## multi_005: Di-μ-hydroxido-bis[diaquacopper(II)]
- **Refcode:** CUAQOH
- **Metals:** Cu,Cu (×2)
- **Bridges:** mu2-OH (×2)
- **Local geometries:** CN = 5,5, shapes = SPY,SPY
- **Notes:** Cu2(OH)2 diamond core with aqua ligands

```
[Metals: M1=Cu; M2=Cu]
[MetalGraph: {M1--M2|relation:bridged|dMM:2.94}]
[LocalSphere: M1=<ShapeBest:SPY|CN:5>{L1:O:1,L2:O:1,L3:O:1,L4:O:1,L5:O:1}; M2=<ShapeBest:SPY|CN:5>{L1:O:1,L2:O:1,L6:O:1,L7:O:1,L8:O:1}]
[Bridges: {L1:O->M1,M2|mu:2}; {L2:O->M1,M2|mu:2}]
[Ligands: L1=[OH-]|L2=[OH-]|L3=O|L4=O|L5=O|L6=O|L7=O|L8=O]
[ID: LocalID[M1]=Cu|SPY|5|O5; LocalID[M2]=Cu|SPY|5|O5; GlobalID=MULTI_005]
```

## multi_006: Di-μ-chlorido-bis[dichlorocobalt(II)]
- **Refcode:** COBRCL
- **Metals:** Co,Co (×2)
- **Bridges:** mu2-Cl (×2)
- **Local geometries:** CN = 4,4, shapes = Td,Td
- **Notes:** Edge-sharing CoCl4 tetrahedra

```
[Metals: M1=Co; M2=Co]
[MetalGraph: {M1--M2|relation:bridged|dMM:3.68}]
[LocalSphere: M1=<ShapeBest:Td|CN:4>{L3:Cl:1,L4:Cl:1,L5:Cl:1,L6:Cl:1}; M2=<ShapeBest:Td|CN:4>{L3:Cl:1,L4:Cl:1,L7:Cl:1,L8:Cl:1}]
[Bridges: {L3:Cl->M1,M2|mu:2}; {L4:Cl->M1,M2|mu:2}]
[Ligands: L3=[Cl-]|L4=[Cl-]|L5=[Cl-]|L6=[Cl-]|L7=[Cl-]|L8=[Cl-]]
[ID: LocalID[M1]=Co|Td|4|Cl4; LocalID[M2]=Co|Td|4|Cl4; GlobalID=MULTI_006]
```

## multi_007: Tetrakis(μ-benzoato)dimanganese(II)
- **Refcode:** MNBZAC
- **Metals:** Mn,Mn (×2)
- **Bridges:** mu2-OBz (×4)
- **Local geometries:** CN = 5,5, shapes = SPY,SPY
- **Notes:** Paddlewheel Mn2(OBz)4 without M-M bond

```
[Metals: M1=Mn; M2=Mn]
[MetalGraph: {M1--M2|relation:bridged|dMM:2.82}]
[LocalSphere: M1=<ShapeBest:SPY|CN:5>{L1:O:1,L2:O:1,L3:O:1,L4:O:1,L5:O:1}; M2=<ShapeBest:SPY|CN:5>{L1:O:2,L2:O:2,L3:O:2,L4:O:2,L6:O:1}]
[Bridges: {L1:O->M1,M2|mu:2}; {L2:O->M1,M2|mu:2}; {L3:O->M1,M2|mu:2}; {L4:O->M1,M2|mu:2}]
[Ligands: L1=OC(=O)c1ccccc1|L2=OC(=O)c1ccccc1|L3=OC(=O)c1ccccc1|L4=OC(=O)c1ccccc1|L5=O|L6=O]
[ID: LocalID[M1]=Mn|SPY|5|O5; LocalID[M2]=Mn|SPY|5|O5; GlobalID=MULTI_007]
```

## multi_008: Di-μ-chlorido-bis[dichloroplatinate(II)]
- **Refcode:** PTCLAM
- **Metals:** Pt,Pt (×2)
- **Bridges:** mu2-Cl (×2)
- **Local geometries:** CN = 4,4, shapes = SP,SP
- **Notes:** Pt2Cl6 edge-sharing square-planar dimer

```
[Metals: M1=Pt; M2=Pt]
[MetalGraph: {M1--M2|relation:bridged|dMM:3.44}]
[LocalSphere: M1=<ShapeBest:SP|CN:4>{L3:Cl:1,L4:Cl:1,L5:Cl:1,L6:Cl:1}; M2=<ShapeBest:SP|CN:4>{L3:Cl:1,L4:Cl:1,L7:Cl:1,L8:Cl:1}]
[Bridges: {L3:Cl->M1,M2|mu:2}; {L4:Cl->M1,M2|mu:2}]
[Ligands: L3=[Cl-]|L4=[Cl-]|L5=[Cl-]|L6=[Cl-]|L7=[Cl-]|L8=[Cl-]]
[ID: LocalID[M1]=Pt|SP|4|Cl4; LocalID[M2]=Pt|SP|4|Cl4; GlobalID=MULTI_008]
```

## multi_009: Tetrakis(μ3-hydroxido)tetrakis[(pyridine)copper(II)]
- **Refcode:** CUDPOM
- **Metals:** Cu,Cu,Cu,Cu (×4)
- **Bridges:** mu3-OH (×4)
- **Local geometries:** CN = 5,5,5,5, shapes = SPY,SPY,SPY,SPY
- **Notes:** Cu4(OH)4 cubane core — prototypical μ3-OH bridging

```
[Metals: M1=Cu; M2=Cu; M3=Cu; M4=Cu]
[MetalGraph: {M1--M2|relation:bridged|dMM:3.02}; {M1--M3|relation:bridged|dMM:3.05}; {M1--M4|relation:bridged|dMM:3.01}; {M2--M3|relation:bridged|dMM:3.04}; {M2--M4|relation:bridged|dMM:3.03}; {M3--M4|relation:bridged|dMM:3.06}]
[LocalSphere: M1=<ShapeBest:SPY|CN:5>{L5:O:1,L6:O:1,L7:O:1,L9:N:1,L13:Cl:1}; M2=<ShapeBest:SPY|CN:5>{L5:O:1,L6:O:1,L8:O:1,L10:N:1,L14:Cl:1}; M3=<ShapeBest:SPY|CN:5>{L5:O:1,L7:O:1,L8:O:1,L11:N:1,L15:Cl:1}; M4=<ShapeBest:SPY|CN:5>{L6:O:1,L7:O:1,L8:O:1,L12:N:1,L16:Cl:1}]
[Bridges: {L5:O->M1,M2,M3|mu:3}; {L6:O->M1,M2,M4|mu:3}; {L7:O->M1,M3,M4|mu:3}; {L8:O->M2,M3,M4|mu:3}]
[Ligands: L5=[OH-]|L6=[OH-]|L7=[OH-]|L8=[OH-]|L9=c1ccncc1|L10=c1ccncc1|L11=c1ccncc1|L12=c1ccncc1|L13=[Cl-]|L14=[Cl-]|L15=[Cl-]|L16=[Cl-]]
[ID: LocalID[M1]=Cu|SPY|5|O3NCl; LocalID[M2]=Cu|SPY|5|O3NCl; LocalID[M3]=Cu|SPY|5|O3NCl; LocalID[M4]=Cu|SPY|5|O3NCl; GlobalID=MULTI_009]
```

## multi_010: Tetrakis(μ-acetato)dizinc(II) dihydrate
- **Refcode:** ZNACET
- **Metals:** Zn,Zn (×2)
- **Bridges:** mu2-OAc (×4)
- **Local geometries:** CN = 5,5, shapes = SPY,SPY
- **Notes:** Zn2(OAc)4·2H2O paddlewheel analog

```
[Metals: M1=Zn; M2=Zn]
[MetalGraph: {M1--M2|relation:bridged|dMM:2.95}]
[LocalSphere: M1=<ShapeBest:SPY|CN:5>{L1:O:1,L2:O:1,L3:O:1,L4:O:1,L5:O:1}; M2=<ShapeBest:SPY|CN:5>{L1:O:2,L2:O:2,L3:O:2,L4:O:2,L6:O:1}]
[Bridges: {L1:O->M1,M2|mu:2}; {L2:O->M1,M2|mu:2}; {L3:O->M1,M2|mu:2}; {L4:O->M1,M2|mu:2}]
[Ligands: L1=CC(=O)[O-]|L2=CC(=O)[O-]|L3=CC(=O)[O-]|L4=CC(=O)[O-]|L5=O|L6=O]
[ID: LocalID[M1]=Zn|SPY|5|O5; LocalID[M2]=Zn|SPY|5|O5; GlobalID=MULTI_010]
```

## multi_011: Tetrakis(μ-acetato)dimolybdenum(II)
- **Refcode:** MOACET02
- **Metals:** Mo,Mo (×2)
- **Bridges:** mu2-OAc (×4)
- **Local geometries:** CN = 5,5, shapes = SPY,SPY
- **Notes:** Mo2(OAc)4 with Mo≡Mo quadruple bond

```
[Metals: M1=Mo; M2=Mo]
[MetalGraph: {M1--M2|relation:bonded|dMM:2.09}]
[LocalSphere: M1=<ShapeBest:SPY|CN:5>{L1:O:1,L2:O:1,L3:O:1,L4:O:1,M2}; M2=<ShapeBest:SPY|CN:5>{L1:O:2,L2:O:2,L3:O:2,L4:O:2,M1}]
[Bridges: {L1:O->M1,M2|mu:2}; {L2:O->M1,M2|mu:2}; {L3:O->M1,M2|mu:2}; {L4:O->M1,M2|mu:2}]
[Ligands: L1=CC(=O)[O-]|L2=CC(=O)[O-]|L3=CC(=O)[O-]|L4=CC(=O)[O-]]
[ID: LocalID[M1]=Mo|SPY|5|O4+Mo; LocalID[M2]=Mo|SPY|5|O4+Mo; GlobalID=MULTI_011]
```

## multi_012: Di-μ-phenoxido-bis[bis(pyridine)nickel(II)]
- **Refcode:** NIPHOX
- **Metals:** Ni,Ni (×2)
- **Bridges:** mu2-OPh (×2)
- **Local geometries:** CN = 6,6, shapes = Oh,Oh
- **Notes:** Bridging phenoxide dinuclear Ni(II) with Oh geometry

```
[Metals: M1=Ni; M2=Ni]
[MetalGraph: {M1--M2|relation:bridged|dMM:3.12}]
[LocalSphere: M1=<ShapeBest:Oh|CN:6>{L3:O:1,L4:O:1,L5:N:1,L6:N:1,L9:Cl:1,L10:Cl:1}; M2=<ShapeBest:Oh|CN:6>{L3:O:1,L4:O:1,L7:N:1,L8:N:1,L11:Cl:1,L12:Cl:1}]
[Bridges: {L3:O->M1,M2|mu:2}; {L4:O->M1,M2|mu:2}]
[Ligands: L3=Oc1ccccc1|L4=Oc1ccccc1|L5=c1ccncc1|L6=c1ccncc1|L7=c1ccncc1|L8=c1ccncc1|L9=[Cl-]|L10=[Cl-]|L11=[Cl-]|L12=[Cl-]]
[ID: LocalID[M1]=Ni|Oh|6|O2N2Cl2; LocalID[M2]=Ni|Oh|6|O2N2Cl2; GlobalID=MULTI_012]
```

## multi_013: μ3-Oxo-hexakis(μ-carboxylato)triiron(III)
- **Refcode:** FEOXTM
- **Metals:** Fe,Fe,Fe (×3)
- **Bridges:** mu3-O;mu2-OAc (×7)
- **Local geometries:** CN = 6,6,6, shapes = Oh,Oh,Oh
- **Notes:** Fe3O oxo-acetate triangle

```
[Metals: M1=Fe; M2=Fe; M3=Fe]
[MetalGraph: {M1--M2|relation:bridged|dMM:3.30}; {M1--M3|relation:bridged|dMM:3.32}; {M2--M3|relation:bridged|dMM:3.31}]
[LocalSphere: M1=<ShapeBest:Oh|CN:6>{L7:O:1,L1:O:1,L2:O:1,L3:O:1,L4:O:1,L8:O:1}; M2=<ShapeBest:Oh|CN:6>{L7:O:1,L1:O:2,L2:O:2,L5:O:1,L6:O:1,L9:O:1}; M3=<ShapeBest:Oh|CN:6>{L7:O:1,L3:O:2,L4:O:2,L5:O:2,L6:O:2,L10:O:1}]
[Bridges: {L7:O->M1,M2,M3|mu:3}; {L1:O->M1,M2|mu:2}; {L2:O->M1,M2|mu:2}; {L3:O->M1,M3|mu:2}; {L4:O->M1,M3|mu:2}; {L5:O->M2,M3|mu:2}; {L6:O->M2,M3|mu:2}]
[Ligands: L1=CC(=O)[O-]|L2=CC(=O)[O-]|L3=CC(=O)[O-]|L4=CC(=O)[O-]|L5=CC(=O)[O-]|L6=CC(=O)[O-]|L7=[O-2]|L8=O|L9=O|L10=O]
[ID: LocalID[M1]=Fe|Oh|6|O6; LocalID[M2]=Fe|Oh|6|O6; LocalID[M3]=Fe|Oh|6|O6; GlobalID=MULTI_013]
```

## multi_014: Bis(μ-dppm)dicopper(I) dichloride
- **Refcode:** CUDPPM
- **Metals:** Cu,Cu (×2)
- **Bridges:** mu2-dppm (×2)
- **Local geometries:** CN = 3,3, shapes = TP,TP
- **Notes:** dppm-bridged dicopper(I), low CN

```
[Metals: M1=Cu; M2=Cu]
[MetalGraph: {M1--M2|relation:bridged|dMM:3.18}]
[LocalSphere: M1=<ShapeBest:TP|CN:3>{L1:P:1,L2:P:1,L3:Cl:1}; M2=<ShapeBest:TP|CN:3>{L1:P:2,L2:P:2,L4:Cl:1}]
[Bridges: {L1:P->M1,M2|mu:2}; {L2:P->M1,M2|mu:2}]
[Ligands: L1=c1ccc(cc1)P(c1ccccc1)CP(c1ccccc1)c1ccccc1|L2=c1ccc(cc1)P(c1ccccc1)CP(c1ccccc1)c1ccccc1|L3=[Cl-]|L4=[Cl-]]
[ID: LocalID[M1]=Cu|TP|3|P2Cl; LocalID[M2]=Cu|TP|3|P2Cl; GlobalID=MULTI_014]
```

## multi_015: Tetrakis(μ3-methoxido)tetramanganese(II) cluster
- **Refcode:** MNOXCB
- **Metals:** Mn,Mn,Mn,Mn (×4)
- **Bridges:** mu3-OMe (×4)
- **Local geometries:** CN = 6,6,6,6, shapes = Oh,Oh,Oh,Oh
- **Notes:** Mn4(OMe)4 cubane core

```
[Metals: M1=Mn; M2=Mn; M3=Mn; M4=Mn]
[MetalGraph: {M1--M2|relation:bridged|dMM:3.18}; {M1--M3|relation:bridged|dMM:3.20}; {M1--M4|relation:bridged|dMM:3.19}; {M2--M3|relation:bridged|dMM:3.21}; {M2--M4|relation:bridged|dMM:3.17}; {M3--M4|relation:bridged|dMM:3.22}]
[LocalSphere: M1=<ShapeBest:Oh|CN:6>{L5:O:1,L6:O:1,L7:O:1,L9:Cl:1,L10:Cl:1,L11:Cl:1}; M2=<ShapeBest:Oh|CN:6>{L5:O:1,L6:O:1,L8:O:1,L12:Cl:1,L13:Cl:1,L14:Cl:1}; M3=<ShapeBest:Oh|CN:6>{L5:O:1,L7:O:1,L8:O:1,L15:Cl:1,L16:Cl:1,L17:Cl:1}; M4=<ShapeBest:Oh|CN:6>{L6:O:1,L7:O:1,L8:O:1,L18:Cl:1,L19:Cl:1,L20:Cl:1}]
[Bridges: {L5:O->M1,M2,M3|mu:3}; {L6:O->M1,M2,M4|mu:3}; {L7:O->M1,M3,M4|mu:3}; {L8:O->M2,M3,M4|mu:3}]
[Ligands: L5=CO|L6=CO|L7=CO|L8=CO|L9=[Cl-]|L10=[Cl-]|L11=[Cl-]|L12=[Cl-]|L13=[Cl-]|L14=[Cl-]|L15=[Cl-]|L16=[Cl-]|L17=[Cl-]|L18=[Cl-]|L19=[Cl-]|L20=[Cl-]]
[ID: LocalID[M1]=Mn|Oh|6|O3Cl3; LocalID[M2]=Mn|Oh|6|O3Cl3; LocalID[M3]=Mn|Oh|6|O3Cl3; LocalID[M4]=Mn|Oh|6|O3Cl3; GlobalID=MULTI_015]
```

