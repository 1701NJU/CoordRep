# CSD Application Casebook

Auto-generated showcase cases for JACS revision.
Each case type has 3–5 candidates; 1 selected for main text, rest for SI.

## Case Type A: Refcode-Family Identity Ladder

### Main Text Selection

```json
{
  "case_type": "A",
  "family": "ACUWOK",
  "n_entries": 12,
  "refcodes": [
    "ACUWOK",
    "ACUWOK01",
    "ACUWOK02",
    "ACUWOK03",
    "ACUWOK04",
    "ACUWOK05"
  ],
  "metal": "Pt",
  "cn": 4,
  "n_unique_L0": 12,
  "n_unique_L1": 2,
  "n_unique_L2": 2,
  "n_unique_L3": 1,
  "shapes": [
    "SP",
    "Td"
  ],
  "has_boundary": false,
  "L0_examples": [
    "[Metal:Pt|ox:+2|d:d8|CN:4]<ShapeBest:Td|Class:good|Delta:3|V:2.92,8.34>{trans:L1…",
    "[Metal:Pt|ox:+2|d:d8|CN:4]<ShapeBest:Td|Class:good|Delta:3|V:2.91,8.39>{trans:L1…",
    "[Metal:Pt|ox:+2|d:d8|CN:4]<ShapeBest:Td|Class:good|Delta:3|V:2.91,8.37>{trans:L1…"
  ],
  "L3_common": "Pt|CN4|Cl;Cl;[H]C1([H])SC([H])([H])C([H])([H])SC([H])([H])C([H])([H])SC1([H])[H]",
  "why_L0_differs": "Geometric detail (stereo/shape) varies across 12 conformers",
  "why_L3_links": "Metal + CN + canonical ligand set is invariant → same chemical identity"
}
```

### SI Selections (4 cases)

**SI Case A-1:**
```json
{
  "case_type": "A",
  "family": "CDMSOR",
  "n_entries": 10,
  "refcodes": [
    "CDMSOR01",
    "CDMSOR02",
    "CDMSOR03",
    "CDMSOR04",
    "CDMSOR05",
    "CDMSOR06"
  ],
  "metal": "Ru",
  "cn": 6,
  "n_unique_L0": 10,
  "n_unique_L1": 4,
  "n_unique_L2": 2,
  "n_unique_L3": 1,
  "shapes": [
    "Oh",
    "TPr"
  ],
  "has_boundary": true,
  "L0_examples": [
    "[Metal:Ru|ox:+2|d:d6|CN:6]<ShapeBest:Oh|Class:dist|Delta:1|V:3.90,5.70>{trans:L1…",
    "[Metal:Ru|ox:+2|d:d6|CN:6]<ShapeBest:Oh|Class:dist|Delta:1|V:3.83,5.46>{trans:L1…",
    "[Metal:Ru|ox:+2|d:d6|CN:6]<ShapeBest:Oh|Class:dist|Delta:2|V:5.49,8.66>{trans:L1…"
  ],
  "L3_common": "Ru|CN6|Cl;Cl;[H]C([H])([H])S(=O)C([H])([H])[H];[H]C([H])([H])S(=O)C([H])([H])[H]…",
  "why_L0_differs": "Geometric detail (stereo/shape) varies across 10 conformers",
  "why_L3_links": "Metal + CN + canonical ligand set is invariant → same chemical identity"
}
```

**SI Case A-2:**
```json
{
  "case_type": "A",
  "family": "DIBFAY",
  "n_entries": 7,
  "refcodes": [
    "DIBFAY",
    "DIBFAY01",
    "DIBFAY02",
    "DIBFAY05",
    "DIBFAY08",
    "DIBFAY10"
  ],
  "metal": "Re",
  "cn": 6,
  "n_unique_L0": 7,
  "n_unique_L1": 4,
  "n_unique_L2": 2,
  "n_unique_L3": 1,
  "shapes": [
    "Oh",
    "TPr"
  ],
  "has_boundary": true,
  "L0_examples": [
    "[Metal:Re|CN:6]<ShapeBest:Oh|Class:dist|Delta:2|V:3.83,6.30>{trans:L1:N:1--L2:C:…",
    "[Metal:Re|CN:6]<ShapeBest:TPr|Class:dist|Delta:1|V:5.22,4.51>{trans:L1:N:1--L2:C…",
    "[Metal:Re|CN:6]<ShapeBest:Oh|Class:dist|Delta:3|V:3.71,9.80>{trans:L1:N:1--L2:C:…"
  ],
  "L3_common": "Re|CN6|C#O;C#O;C#O;O=N=O;[H]C1NC(C2NC([H])C([H])C(OC([H])([H])[H])C2[H])C([H])C(…",
  "why_L0_differs": "Geometric detail (stereo/shape) varies across 7 conformers",
  "why_L3_links": "Metal + CN + canonical ligand set is invariant → same chemical identity"
}
```

**SI Case A-3:**
```json
{
  "case_type": "A",
  "family": "COVLUX",
  "n_entries": 6,
  "refcodes": [
    "COVLUX",
    "COVLUX01",
    "COVLUX02",
    "COVLUX03",
    "COVLUX05",
    "COVLUX06"
  ],
  "metal": "Ni",
  "cn": 4,
  "n_unique_L0": 6,
  "n_unique_L1": 2,
  "n_unique_L2": 2,
  "n_unique_L3": 1,
  "shapes": [
    "SP",
    "Td"
  ],
  "has_boundary": false,
  "L0_examples": [
    "[Metal:Ni|ox:+2|d:d8|CN:4]<ShapeBest:Td|Class:good|Delta:3|V:2.92,8.04>{trans:L1…",
    "[Metal:Ni|ox:+2|d:d8|CN:4]<ShapeBest:SP|Class:ideal|Delta:2|V:2.91,0.03>{trans:L…",
    "[Metal:Ni|ox:+2|d:d8|CN:4]<ShapeBest:SP|Class:ideal|Delta:2|V:2.92,0.04>{trans:L…"
  ],
  "L3_common": "Ni|CN4|S=c1sc(S)c(S)s1;S=c1sc(S)c(S)s1",
  "why_L0_differs": "Geometric detail (stereo/shape) varies across 6 conformers",
  "why_L3_links": "Metal + CN + canonical ligand set is invariant → same chemical identity"
}
```

**SI Case A-4:**
```json
{
  "case_type": "A",
  "family": "AWEPAW",
  "n_entries": 5,
  "refcodes": [
    "AWEPAW",
    "AWEPAW01",
    "AWEPAW02",
    "AWEPAW03",
    "AWEPAW04"
  ],
  "metal": "Cu",
  "cn": 4,
  "n_unique_L0": 5,
  "n_unique_L1": 2,
  "n_unique_L2": 2,
  "n_unique_L3": 1,
  "shapes": [
    "SP",
    "Td"
  ],
  "has_boundary": false,
  "L0_examples": [
    "[Metal:Cu|ox:+2|d:d9|CN:4]<ShapeBest:Td|Class:good|Delta:3|V:2.69,8.79>{trans:L1…",
    "[Metal:Cu|ox:+2|d:d9|CN:4]<ShapeBest:SP|Class:ideal|Delta:2|V:2.67,0.04>{trans:L…",
    "[Metal:Cu|ox:+2|d:d9|CN:4]<ShapeBest:SP|Class:ideal|Delta:2|V:2.67,0.03>{trans:L…"
  ],
  "L3_common": "Cu|CN4|[H]C(O)C([H])C(O)C(C([H])([H])[H])(C([H])([H])[H])C([H])([H])[H];[H]C(O)C…",
  "why_L0_differs": "Geometric detail (stereo/shape) varies across 5 conformers",
  "why_L3_links": "Metal + CN + canonical ligand set is invariant → same chemical identity"
}
```

---

## Case Type B: Boundary Geometry Case

### Main Text Selection

```json
{
  "case_type": "B",
  "refcode": "AFOSIA",
  "family": "AFOSIA",
  "metal": "Cr",
  "cn": 6,
  "shape_top1": "TPr",
  "shape_top2": "Oh",
  "shape_class": "dist",
  "V_top1": 4.49,
  "V_top2": 4.48,
  "V_competition": 0.01,
  "explanation": "CN=6 Cr complex: CShM(TPr)=4.49 vs CShM(Oh)=4.48 (Δ=0.01). Class=dist. Close competition means discrete shape label is unreliable → continuous boundary identity needed."
}
```

### SI Selections (4 cases)

**SI Case B-1:**
```json
{
  "case_type": "B",
  "refcode": "ALACOH",
  "family": "ALACOH",
  "metal": "Ru",
  "cn": 6,
  "shape_top1": "Oh",
  "shape_top2": "TPr",
  "shape_class": "dist",
  "V_top1": 4.33,
  "V_top2": 4.34,
  "V_competition": 0.01,
  "explanation": "CN=6 Ru complex: CShM(Oh)=4.33 vs CShM(TPr)=4.34 (Δ=0.01). Class=dist. Close competition means discrete shape label is unreliable → continuous boundary identity needed."
}
```

**SI Case B-2:**
```json
{
  "case_type": "B",
  "refcode": "APZBFE",
  "family": "APZBFE",
  "metal": "Fe",
  "cn": 6,
  "shape_top1": "Oh",
  "shape_top2": "TPr",
  "shape_class": "dist",
  "V_top1": 5.48,
  "V_top2": 5.49,
  "V_competition": 0.01,
  "explanation": "CN=6 Fe complex: CShM(Oh)=5.48 vs CShM(TPr)=5.49 (Δ=0.01). Class=dist. Close competition means discrete shape label is unreliable → continuous boundary identity needed."
}
```

**SI Case B-3:**
```json
{
  "case_type": "B",
  "refcode": "ASEJOA",
  "family": "ASEJOA",
  "metal": "Cu",
  "cn": 6,
  "shape_top1": "TPr",
  "shape_top2": "Oh",
  "shape_class": "dist",
  "V_top1": 5.48,
  "V_top2": 5.47,
  "V_competition": 0.01,
  "explanation": "CN=6 Cu complex: CShM(TPr)=5.48 vs CShM(Oh)=5.47 (Δ=0.01). Class=dist. Close competition means discrete shape label is unreliable → continuous boundary identity needed."
}
```

**SI Case B-4:**
```json
{
  "case_type": "B",
  "refcode": "ASPTCR",
  "family": "ASPTCR",
  "metal": "Re",
  "cn": 5,
  "shape_top1": "SPY",
  "shape_top2": "TBP",
  "shape_class": "dist",
  "V_top1": 5.89,
  "V_top2": 5.88,
  "V_competition": 0.01,
  "explanation": "CN=5 Re complex: CShM(SPY)=5.89 vs CShM(TBP)=5.88 (Δ=0.01). Class=dist. Close competition means discrete shape label is unreliable → continuous boundary identity needed."
}
```

---

## Case Type C: Tool B Repair Case (MLM vs Rule)

### Main Text Selection

```json
{
  "case_type": "C",
  "refcode": "AGOTIA",
  "family": "AGOTIA",
  "metal": "Re",
  "cn": 7,
  "corruption_type": "missing_bracket",
  "clean_excerpt": "[Metal:Re|CN:7]{trans:L1:O:1--L3:O:1}{trans:L2:O:2--L4:Cl:1}|L1=[H]c1oc(C([H])([H])[H])c(O)c(=O)c1[H]|L2=[H]c1oc(C([H])(…",
  "corrupted_excerpt": "[Metal:Re|CN:7]{trans:L1:O:1--L3:O:1}{trans:L2:O:2--L4:Cl:1}|L1=[H]c1oc(C([H])([H)[H])c(O)c(=O)c1[H]|L2=[H]c1oc(C([H])([…",
  "rule_repaired_excerpt": "[Metal:Re|CN:7]{trans:L1:O:1--L3:O:1}{trans:L2:O:2--L4:Cl:1}|L1=[H]c1oc(C([H])([H][H])c(O)c(=O)c1[H]|L2=[H]c1oc(C([H])([…",
  "rule_valid": true,
  "rule_exact": false,
  "rule_field_score": "5/6",
  "rule_field_detail": {
    "metal": true,
    "cn": true,
    "ox": true,
    "shape": true,
    "ligands": false,
    "stereo": true
  },
  "corruption_position": 81,
  "has_stereo": true,
  "string_length": 154,
  "explanation": "missing_bracket corruption at position 81 makes string invalid. CoordRep grammar-aware rule repair restores validity and preserves 5/6 semantic fields. This deterministic repair is enabled by CoordRep's explicit bracket/delimiter/field grammar."
}
```

### SI Selections (4 cases)

**SI Case C-1:**
```json
{
  "case_type": "C",
  "refcode": "AWACEJ",
  "family": "AWACEJ",
  "metal": "Pd",
  "cn": 4,
  "corruption_type": "missing_bracket",
  "clean_excerpt": "[Metal:Pd|ox:+2|d:d8|CN:4]<ShapeBest:SP|Class:ideal|Delta:2|V:3.28,0.49>{trans:L1:N:1--L2:Cl:1}{trans:L1:N:2--L1:S:1}{fm…",
  "corrupted_excerpt": "[Metal:Pd|ox:+2|d:d8|CN:4]<ShapeBest:SP|Class:ideal|Delta:2|V:3.28,0.49>{trans:L1:N:1--L2:Cl:1}{trans:L1:N:2--L1:S:1}{fm…",
  "rule_repaired_excerpt": "[Metal:Pd|ox:+2|d:d8|CN:4]<ShapeBest:SP|Class:ideal|Delta:2|V:3.28,0.49>{trans:L1:N:1--L2:Cl:1}{trans:L1:N:2--L1:S:1}{fm…",
  "rule_valid": true,
  "rule_exact": false,
  "rule_field_score": "5/6",
  "rule_field_detail": {
    "metal": true,
    "cn": true,
    "ox": true,
    "shape": true,
    "ligands": false,
    "stereo": true
  },
  "corruption_position": 163,
  "has_stereo": true,
  "string_length": 249,
  "explanation": "missing_bracket corruption at position 163 makes string invalid. CoordRep grammar-aware rule repair restores validity and preserves 5/6 semantic fields. This deterministic repair is enabled by CoordRep's explicit bracket/delimiter/field grammar."
}
```

**SI Case C-2:**
```json
{
  "case_type": "C",
  "refcode": "DAJJAC",
  "family": "DAJJAC",
  "metal": "Rh",
  "cn": 6,
  "corruption_type": "missing_bracket",
  "clean_excerpt": "[Metal:Rh|ox:+3|d:d6|CN:6]<ShapeBest:Oh|Class:ideal|Delta:3|V:0.96,6.40>{trans:L1:C:1--L2:Cl:1}{trans:L1:N:1--L3:Cl:1}|L…",
  "corrupted_excerpt": "[Metal:Rh|ox:+3|d:d6|CN:6]<ShapeBest:Oh|Class:ideal|Delta:3|V:0.96,6.40>{trans:L1:C:1--L2:Cl:1}{trans:L1:N:1--L3:Cl:1}|L…",
  "rule_repaired_excerpt": "[Metal:Rh|ox:+3|d:d6|CN:6]<ShapeBest:Oh|Class:ideal|Delta:3|V:0.96,6.40>{trans:L1:C:1--L2:Cl:1}{trans:L1:N:1--L3:Cl:1}|L…",
  "rule_valid": true,
  "rule_exact": false,
  "rule_field_score": "5/6",
  "rule_field_detail": {
    "metal": true,
    "cn": true,
    "ox": true,
    "shape": true,
    "ligands": false,
    "stereo": true
  },
  "corruption_position": 154,
  "has_stereo": true,
  "string_length": 257,
  "explanation": "missing_bracket corruption at position 154 makes string invalid. CoordRep grammar-aware rule repair restores validity and preserves 5/6 semantic fields. This deterministic repair is enabled by CoordRep's explicit bracket/delimiter/field grammar."
}
```

**SI Case C-3:**
```json
{
  "case_type": "C",
  "refcode": "CIZJAZ",
  "family": "CIZJAZ",
  "metal": "Zn",
  "cn": 6,
  "corruption_type": "missing_bracket",
  "clean_excerpt": "[Metal:Zn|ox:+2|d:d10|CN:6]<ShapeBest:Oh|Class:dist|Delta:2|V:3.70,8.51>{trans:L1:N:1--L2:N:1}|L1=[H]OC1C([H])C([H])C([H…",
  "corrupted_excerpt": "[Metal:Zn|ox:+2|d:d10|CN:6]<ShapeBest:Oh|Class:dist|Delta:2|V:3.70,8.51>{trans:L1:N:1--L2:N:1}|L1=[H]OC1C([H])C([H])C([H…",
  "rule_repaired_excerpt": "[Metal:Zn|ox:+2|d:d10|CN:6]<ShapeBest:Oh|Class:dist|Delta:2|V:3.70,8.51>{trans:L1:N:1--L2:N:1}|L1=[H]OC1C([H])C([H])C([H…",
  "rule_valid": true,
  "rule_exact": false,
  "rule_field_score": "5/6",
  "rule_field_detail": {
    "metal": true,
    "cn": true,
    "ox": true,
    "shape": true,
    "ligands": false,
    "stereo": true
  },
  "corruption_position": 147,
  "has_stereo": true,
  "string_length": 259,
  "explanation": "missing_bracket corruption at position 147 makes string invalid. CoordRep grammar-aware rule repair restores validity and preserves 5/6 semantic fields. This deterministic repair is enabled by CoordRep's explicit bracket/delimiter/field grammar."
}
```

**SI Case C-4:**
```json
{
  "case_type": "C",
  "refcode": "DARYAZ",
  "family": "DARYAZ",
  "metal": "Pd",
  "cn": 5,
  "corruption_type": "missing_bracket",
  "clean_excerpt": "[Metal:Pd|ox:+2|d:d8|CN:5]<ShapeBest:TBP|Class:irreg|Delta:2|V:8.64,13.19>{trans:L1:N:1--L1:N:4}{trans:L1:N:2--L1:N:3}|L…",
  "corrupted_excerpt": "[Metal:Pd|ox:+2|d:d8|CN:5]<ShapeBest:TBP|Class:irreg|Delta:2|V:8.64,13.19>{trans:L1:N:1--L1:N:4}{trans:L1:N:2--L1:N:3}|L…",
  "rule_repaired_excerpt": "[Metal:Pd|ox:+2|d:d8|CN:5]<ShapeBest:TBP|Class:irreg|Delta:2|V:8.64,13.19>{trans:L1:N:1--L1:N:4}{trans:L1:N:2--L1:N:3}|L…",
  "rule_valid": true,
  "rule_exact": false,
  "rule_field_score": "5/6",
  "rule_field_detail": {
    "metal": true,
    "cn": true,
    "ox": true,
    "shape": true,
    "ligands": false,
    "stereo": true
  },
  "corruption_position": 146,
  "has_stereo": true,
  "string_length": 285,
  "explanation": "missing_bracket corruption at position 146 makes string invalid. CoordRep grammar-aware rule repair restores validity and preserves 5/6 semantic fields. This deterministic repair is enabled by CoordRep's explicit bracket/delimiter/field grammar."
}
```

---

## Case Type D: Stereo Semantic Consistency

### Main Text Selection

```json
{
  "case_type": "D",
  "refcode": "CIJWUO",
  "family": "CIJWUO",
  "metal": "Cu",
  "cn": 5,
  "real_stereo_token": "{trans:L1:N:1--L2:N:1}",
  "decoy_stereo_token": "{cis:L1:N:1--L2:N:1}",
  "n_stereo_relations": 1,
  "explanation": "Stereo flip (trans→cis): semantic field-level change that affects ligand arrangement without changing raw geometry descriptors (bond lengths, angles). CoordRep-Ranker can detect via learned stereo compatibility.",
  "ranker_score_real": -2.479,
  "ranker_score_decoy": -3.1637,
  "ranker_correct": true
}
```

### SI Selections (4 cases)

**SI Case D-1:**
```json
{
  "case_type": "D",
  "refcode": "BATTIZ",
  "family": "BATTIZ",
  "metal": "Pd",
  "cn": 4,
  "real_stereo_token": "{trans:L1:N:1--L2:Cl:1}",
  "decoy_stereo_token": "{cis:L1:N:1--L2:Cl:1}",
  "n_stereo_relations": 2,
  "explanation": "Stereo flip (trans→cis): semantic field-level change that affects ligand arrangement without changing raw geometry descriptors (bond lengths, angles). CoordRep-Ranker can detect via learned stereo compatibility.",
  "ranker_score_real": -4.653,
  "ranker_score_decoy": -4.979,
  "ranker_correct": true
}
```

**SI Case D-2:**
```json
{
  "case_type": "D",
  "refcode": "CIRRUS",
  "family": "CIRRUS",
  "metal": "Pt",
  "cn": 4,
  "real_stereo_token": "{trans:L1:C:1--L2:C:1}",
  "decoy_stereo_token": "{cis:L1:C:1--L2:C:1}",
  "n_stereo_relations": 2,
  "explanation": "Stereo flip (trans→cis): semantic field-level change that affects ligand arrangement without changing raw geometry descriptors (bond lengths, angles). CoordRep-Ranker can detect via learned stereo compatibility.",
  "ranker_score_real": 3.9317,
  "ranker_score_decoy": 3.5915,
  "ranker_correct": true
}
```

**SI Case D-3:**
```json
{
  "case_type": "D",
  "refcode": "CEHZIZ",
  "family": "CEHZIZ",
  "metal": "Mn",
  "cn": 6,
  "real_stereo_token": "{trans:L1:O:1--L2:O:1}",
  "decoy_stereo_token": "{cis:L1:O:1--L2:O:1}",
  "n_stereo_relations": 3,
  "explanation": "Stereo flip (trans→cis): semantic field-level change that affects ligand arrangement without changing raw geometry descriptors (bond lengths, angles). CoordRep-Ranker can detect via learned stereo compatibility.",
  "ranker_score_real": -0.7165,
  "ranker_score_decoy": -0.8132,
  "ranker_correct": true
}
```

**SI Case D-4:**
```json
{
  "case_type": "D",
  "refcode": "BECDUK",
  "family": "BECDUK",
  "metal": "Cu",
  "cn": 2,
  "real_stereo_token": "{trans:L1:C:1--L2:N:1}",
  "decoy_stereo_token": "{cis:L1:C:1--L2:N:1}",
  "n_stereo_relations": 1,
  "explanation": "Stereo flip (trans→cis): semantic field-level change that affects ligand arrangement without changing raw geometry descriptors (bond lengths, angles). CoordRep-Ranker can detect via learned stereo compatibility.",
  "ranker_score_real": -10.5309,
  "ranker_score_decoy": -11.1699,
  "ranker_correct": true
}
```

---

## Summary

| Case Type | Main Text | SI Count | Key Finding |
|-----------|-----------|----------|-------------|
| A | ACUWOK | 4 | L0 varies across conformers but L3 identity is stable → CoordRep-ID links families |
| B | AFOSIA | 4 | Boundary geometry requires continuous shape descriptor, not discrete label |
| C | AGOTIA | 4 | CoordRep-MLM repairs token-level corruption that rules cannot restore exactly |
| D | CIJWUO | 4 | CoordRep-Ranker detects stereo semantic inconsistency via learned field compatibility |
