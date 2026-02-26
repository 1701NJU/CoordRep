# Supplementary Information

## CoordRep: A Canonical Representation for Coordination Complexes Enabling Chemical Language Modeling

---

## Table of Contents

- [S1. Data Sources, Filtering Rules, and Statistics](#s1-data-sources-filtering-rules-and-statistics)
- [S2. Transparent Conversion Pipeline](#s2-transparent-conversion-pipeline)
- [S3. Continuous Shape Measure (CShM) Definitions](#s3-continuous-shape-measure-cshm-definitions)
- [S4. Canonicalization and Determinism](#s4-canonicalization-and-determinism)
- [S5. Tokenizer Specification](#s5-tokenizer-specification)
- [S6. Masking Protocol and Task Definitions](#s6-masking-protocol-and-task-definitions)
- [S7. Synthetic Corruption Suite (Tool B)](#s7-synthetic-corruption-suite-tool-b)
- [S8. Baselines](#s8-baselines)
- [S9. Training Details](#s9-training-details)
- [S10. Reproducibility Package](#s10-reproducibility-package)

---

## S1. Data Sources, Filtering Rules, and Statistics

### S1.1 Data Sources

| Dataset | Version | Access Date | License | Records |
|---------|---------|-------------|---------|---------|
| tmQM | v1.0 | 2024-01-15 | CC BY 4.0 | 86,665 |
| COD (Crystallography Open Database) | rev. 2024.01 | 2024-01-20 | Public Domain | 523,142 |

**tmQM Reference**: Balcells, D. & Skjelstad, B. B. tmQM Dataset—Quantum Geometries and Properties of 86k Transition Metal Complexes. *J. Chem. Inf. Model.* 60, 6135–6146 (2020).

**COD Reference**: Gražulis, S. et al. Crystallography Open Database—an open-access collection of crystal structures. *J. Appl. Cryst.* 42, 726–729 (2009).

### S1.2 Filtering Rules

We apply the following sequential filtering rules to ensure data quality:

| Rule ID | Rule Description | Rationale |
|---------|------------------|-----------|
| F1 | Single metal center only | Exclude multinuclear complexes (data complexity) |
| F2 | Hapticity ≤ 1 for all ligands | Exclude π-coordinated ligands (ambiguous representation) |
| F3 | Coordination number 2 ≤ CN ≤ 14 | Focus on common coordination environments |
| F4 | No fragmented structures | Ensure complete coordination sphere |
| F5 | No positional disorder | Avoid ambiguous atomic positions |
| F6 | R-factor ≤ 0.10 (COD only) | Quality threshold for experimental structures |
| F7 | No missing donor atoms | Ensure complete ligand specification |
| F8 | Valid SMILES generation | Ensure RDKit compatibility |

### S1.3 Data Flow Statistics

```
tmQM Raw:                    86,665 structures
  ├─ After F1 (single metal): 82,341 (95.0%)
  ├─ After F2 (hapticity):    71,856 (82.9%)
  ├─ After F3 (CN range):     70,124 (80.9%)
  ├─ After F4-F7:             63,847 (73.7%)
  └─ Final (valid SMILES):    61,234 (70.7%)

COD Raw:                     523,142 structures
  ├─ Transition metals only: 127,845 (24.4%)
  ├─ After F1-F6:             42,156 (8.1%)
  └─ Final (valid CoordRep):  38,721 (7.4%)

Combined Dataset:            99,955 structures
```

### S1.4 Stratified Split

We perform stratified sampling to ensure balanced representation across metal elements and coordination numbers:

```python
# Stratification procedure
from sklearn.model_selection import train_test_split

# Create joint stratification key
df['strat_key'] = df['metal_element'] + '_CN' + df['cn'].astype(str)

# Split with stratification
train, temp = train_test_split(df, test_size=0.2, stratify=df['strat_key'], random_state=42)
val, test = train_test_split(temp, test_size=0.5, stratify=temp['strat_key'], random_state=42)

# Final sizes
# Train: 79,964 (80%)
# Val:    9,996 (10%)  
# Test:   9,995 (10%)
```

**Stratification Statistics:**

| Metal Group | Train | Val | Test | Total |
|-------------|-------|-----|------|-------|
| 3d metals | 48,234 | 6,029 | 6,029 | 60,292 |
| 4d metals | 19,856 | 2,482 | 2,482 | 24,820 |
| 5d metals | 11,874 | 1,485 | 1,484 | 14,843 |

| CN | Train | Val | Test |
|----|-------|-----|------|
| 4 | 8,547 | 1,068 | 1,068 |
| 5 | 12,794 | 1,599 | 1,599 |
| 6 | 31,986 | 3,998 | 3,998 |
| 7 | 10,395 | 1,299 | 1,300 |
| 8 | 9,596 | 1,200 | 1,199 |
| Other | 6,646 | 832 | 831 |

---

## S2. Transparent Conversion Pipeline

### S2.1 Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        Raw Structure → CoordRep Pipeline                      │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                               │
│  ┌─────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────────┐   │
│  │ CIF/XYZ │───▶│   Parser    │───▶│ Metal ID &  │───▶│ First Coord.    │   │
│  │   /SDF  │    │ (gemmi/RDKit)│    │ Validation  │    │ Sphere Detect   │   │
│  └─────────┘    └─────────────┘    └─────────────┘    └────────┬────────┘   │
│                                                                  │            │
│                                                                  ▼            │
│  ┌─────────────────┐    ┌─────────────────┐    ┌─────────────────────────┐  │
│  │ CoordRep String │◀───│  Canonicalize   │◀───│  Ligand Extraction &    │  │
│  │    Output       │    │  & Tokenize     │    │  CShM Calculation       │  │
│  └─────────────────┘    └─────────────────┘    └─────────────────────────┘  │
│                                                                               │
└─────────────────────────────────────────────────────────────────────────────┘
```

### S2.2 Input Formats and Parsers

| Format | Parser | Version | Notes |
|--------|--------|---------|-------|
| CIF | gemmi | 0.6.4 | Primary format for COD |
| XYZ | custom | - | Coordinate-only format |
| SDF/MOL | RDKit | 2023.09.1 | For tmQM structures |

### S2.3 Metal Identification

We identify transition metals using the following criteria:

```python
TRANSITION_METALS = {
    # 3d metals
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    # 4d metals
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    # 5d metals  
    'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    # Lanthanides (selected)
    'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu'
}

def identify_metal_center(atoms):
    """Identify metal center(s) in structure."""
    metals = [a for a in atoms if a.element in TRANSITION_METALS]
    if len(metals) == 0:
        raise ValueError("No transition metal found")
    if len(metals) > 1:
        raise ValueError("Multinuclear complex excluded")
    return metals[0]
```

### S2.4 First Coordination Sphere Detection

The coordination sphere is determined using a distance-based criterion:

```
d(M−X) ≤ α × (r_M + r_X) + δ
```

Where:
- `d(M−X)`: Distance between metal M and potential donor atom X
- `r_M`, `r_X`: Covalent radii (from CSD/Cordero dataset)
- `α = 1.15`: Scaling factor
- `δ = 0.4 Å`: Tolerance parameter

**Sensitivity Analysis:**

| α | δ (Å) | Mean CN | False Positives | False Negatives |
|---|-------|---------|-----------------|-----------------|
| 1.10 | 0.3 | 5.42 | 2.1% | 8.3% |
| **1.15** | **0.4** | **5.87** | **3.4%** | **2.1%** |
| 1.20 | 0.5 | 6.34 | 7.8% | 0.9% |

We select α=1.15, δ=0.4 Å as the optimal balance.

### S2.5 Failure and Ambiguity Handling

| Failure Type | Detection | Resolution |
|--------------|-----------|------------|
| Multiple metal candidates | Count check | Exclude structure |
| Fragmented structure | Graph connectivity | Exclude structure |
| Missing hydrogens | Valence check | Add hydrogens (RDKit) |
| Positional disorder | Occupancy < 1.0 | Exclude structure |
| Ambiguous CN | Multiple valid spheres | Use stricter threshold |

### S2.6 Pseudocode

```python
def structure_to_coordrep(structure_file):
    # Step 1: Parse structure
    atoms, coords = parse_structure(structure_file)
    
    # Step 2: Identify metal center
    metal = identify_metal_center(atoms)
    
    # Step 3: Detect first coordination sphere
    donors = []
    for atom in atoms:
        if atom != metal:
            dist = distance(metal.coord, atom.coord)
            threshold = 1.15 * (covalent_radius[metal.element] + 
                               covalent_radius[atom.element]) + 0.4
            if dist <= threshold:
                donors.append(atom)
    
    # Step 4: Extract ligands (subgraph from each donor)
    ligands = extract_ligands(atoms, metal, donors)
    
    # Step 5: Calculate CShM
    cshm_results = calculate_cshm(metal.coord, [d.coord for d in donors])
    
    # Step 6: Canonicalize and generate CoordRep
    coordrep = canonicalize_and_encode(metal, donors, ligands, cshm_results)
    
    return coordrep
```

---

## S3. Continuous Shape Measure (CShM) Definitions

### S3.1 Mathematical Definition

The Continuous Shape Measure (CShM) quantifies the deviation of a coordination geometry from an ideal reference polyhedron:

```
S(P₀) = min[100 × Σᵢ|qᵢ - pᵢ|² / Σᵢ|qᵢ - q₀|²]
```

Where:
- `P₀`: Reference polyhedron (ideal geometry)
- `qᵢ`: Coordinates of actual donor atoms (after optimal alignment)
- `pᵢ`: Coordinates of vertices of reference polyhedron
- `q₀`: Centroid of actual donor atoms

The minimization is performed over all possible vertex permutations and rotations.

### S3.2 Alignment and Permutation

**Alignment procedure:**
1. Center both actual and reference coordinates at origin
2. Compute optimal rotation using Kabsch algorithm
3. Enumerate all vertex permutations (n! for CN=n, pruned using Hungarian algorithm)
4. Select permutation minimizing S(P₀)

```python
def compute_cshm(actual_coords, reference_coords):
    """
    Compute CShM for a given reference polyhedron.
    
    Args:
        actual_coords: (N, 3) array of donor atom coordinates
        reference_coords: (N, 3) array of ideal polyhedron vertices
    
    Returns:
        S: CShM value (0 = perfect match, larger = more distorted)
    """
    # Center coordinates
    actual_centered = actual_coords - actual_coords.mean(axis=0)
    ref_centered = reference_coords - reference_coords.mean(axis=0)
    
    # Normalize reference to unit size
    ref_scale = np.sqrt((ref_centered ** 2).sum() / len(ref_centered))
    ref_norm = ref_centered / ref_scale
    
    # Find optimal permutation and rotation
    min_S = float('inf')
    for perm in permutations(range(len(actual_coords))):
        perm_actual = actual_centered[list(perm)]
        
        # Kabsch alignment
        R = kabsch_rotation(perm_actual, ref_norm)
        aligned = perm_actual @ R
        
        # Compute S
        S = 100 * np.sum((aligned - ref_norm) ** 2) / np.sum(ref_norm ** 2)
        min_S = min(min_S, S)
    
    return min_S
```

### S3.3 Reference Polyhedra

| CN | Label | Geometry | Point Group | Vertices |
|----|-------|----------|-------------|----------|
| 2 | L-2 | Linear | D∞h | See Table S3.1 |
| 3 | TP-3 | Trigonal Planar | D3h | See Table S3.1 |
| 4 | **T-4** | **Tetrahedral** | **Td** | See Table S3.1 |
| 4 | **SP-4** | **Square Planar** | **D4h** | See Table S3.1 |
| 5 | **TBPY-5** | **Trigonal Bipyramidal** | **D3h** | See Table S3.1 |
| 5 | **SPY-5** | **Square Pyramidal** | **C4v** | See Table S3.1 |
| 6 | **OC-6** | **Octahedral** | **Oh** | See Table S3.1 |
| 6 | TPR-6 | Trigonal Prismatic | D3h | See Table S3.1 |
| 7 | PBPY-7 | Pentagonal Bipyramidal | D5h | See Table S3.1 |
| 7 | COC-7 | Capped Octahedral | C3v | See Table S3.1 |
| 8 | SAPR-8 | Square Antiprismatic | D4d | See Table S3.1 |
| 8 | CU-8 | Cubic | Oh | See Table S3.1 |
| 9 | TCA-9 | Tricapped Trigonal Prismatic | D3h | See Table S3.1 |

### S3.4 Reference Vertex Coordinates (Table S3.1)

**Tetrahedral (T-4, Td):**
```
V1: ( 1.0000,  0.0000, -0.7071)
V2: (-1.0000,  0.0000, -0.7071)
V3: ( 0.0000,  1.0000,  0.7071)
V4: ( 0.0000, -1.0000,  0.7071)
```

**Square Planar (SP-4, D4h):**
```
V1: ( 1.0000,  0.0000,  0.0000)
V2: ( 0.0000,  1.0000,  0.0000)
V3: (-1.0000,  0.0000,  0.0000)
V4: ( 0.0000, -1.0000,  0.0000)
```

**Octahedral (OC-6, Oh):**
```
V1: ( 1.0000,  0.0000,  0.0000)
V2: (-1.0000,  0.0000,  0.0000)
V3: ( 0.0000,  1.0000,  0.0000)
V4: ( 0.0000, -1.0000,  0.0000)
V5: ( 0.0000,  0.0000,  1.0000)
V6: ( 0.0000,  0.0000, -1.0000)
```

**Trigonal Bipyramidal (TBPY-5, D3h):**
```
V1: ( 0.0000,  0.0000,  1.0000)  # Axial
V2: ( 0.0000,  0.0000, -1.0000)  # Axial
V3: ( 1.0000,  0.0000,  0.0000)  # Equatorial
V4: (-0.5000,  0.8660,  0.0000)  # Equatorial
V5: (-0.5000, -0.8660,  0.0000)  # Equatorial
```

**Square Pyramidal (SPY-5, C4v):**
```
V1: ( 0.0000,  0.0000,  1.0000)  # Apical
V2: ( 1.0000,  0.0000,  0.0000)  # Basal
V3: ( 0.0000,  1.0000,  0.0000)  # Basal
V4: (-1.0000,  0.0000,  0.0000)  # Basal
V5: ( 0.0000, -1.0000,  0.0000)  # Basal
```

### S3.5 Normalization and Discretization

**Shape classification rule:**
```
Best Shape = argmin{S(P) for P in reference_polyhedra}
```

**Discretization thresholds:**
- S < 1.0: "Perfect" geometry
- 1.0 ≤ S < 3.0: "Near-ideal" geometry
- 3.0 ≤ S < 8.0: "Distorted" geometry
- S ≥ 8.0: "Severely distorted" or intermediate

**Token encoding:**
```
<ShapeBest:{LABEL}|Class:{CLASS}|Delta:{DELTA}|V:{S1},{S2}>
```
Where:
- `LABEL`: Best-matching polyhedron (e.g., "Oh", "SP")
- `CLASS`: "ideal" / "dist" / "severe"
- `DELTA`: Difference between top-2 CShM values
- `S1,S2`: Top-2 CShM values

---

## S4. Canonicalization and Determinism

### S4.1 Ligand Segmentation

Ligands are extracted as connected subgraphs rooted at donor atoms:

```python
def extract_ligands(mol_graph, metal, donors):
    """Extract ligands as subgraphs."""
    ligands = []
    visited = set()
    
    for donor in donors:
        if donor in visited:
            continue
        
        # BFS from donor, excluding metal
        ligand_atoms = bfs_subgraph(mol_graph, donor, exclude={metal})
        
        # Mark all atoms as visited
        visited.update(ligand_atoms)
        
        # Store ligand with donor info
        ligands.append({
            'atoms': ligand_atoms,
            'donor': donor,
            'denticity': sum(1 for d in donors if d in ligand_atoms)
        })
    
    return ligands
```

### S4.2 Ligand Ordering Priority Rules

Ligands are sorted using the following priority rules (applied in order):

| Priority | Rule | Description | Tie-break |
|----------|------|-------------|-----------|
| 1 | Donor element | Periodic table order: C < N < O < F < P < S < Cl < Br < I | → Rule 2 |
| 2 | Denticity | Monodentate (1) < Bidentate (2) < ... | → Rule 3 |
| 3 | Heavy atom count | Smaller first | → Rule 4 |
| 4 | Formal charge | More negative first (-2 < -1 < 0 < +1) | → Rule 5 |
| 5 | Molecular weight | Lighter first | → Rule 6 |
| 6 | Canonical SMILES | Lexicographic order (RDKit 2023.09.1) | → Rule 7 |
| 7 | Hash tie-break | SHA256 hash of SMILES | Deterministic |

```python
DONOR_PRIORITY = {'C': 0, 'N': 1, 'O': 2, 'F': 3, 'P': 4, 'S': 5, 'Cl': 6, 'Br': 7, 'I': 8}

def ligand_sort_key(ligand):
    """Generate sort key for ligand ordering."""
    return (
        DONOR_PRIORITY.get(ligand['donor_element'], 99),  # Rule 1
        ligand['denticity'],                               # Rule 2
        ligand['heavy_atom_count'],                        # Rule 3
        ligand['formal_charge'],                           # Rule 4
        ligand['molecular_weight'],                        # Rule 5
        ligand['canonical_smiles'],                        # Rule 6
        hashlib.sha256(ligand['canonical_smiles'].encode()).hexdigest()  # Rule 7
    )
```

### S4.3 Identical Ligand Disambiguation (Stereochemistry)

For identical ligands, we determine ordering based on spatial arrangement:

**Method: Donor-Metal-Donor Angle Classification**

```python
def classify_stereo_arrangement(metal_coord, donor_coords, ligand_assignments):
    """
    Classify stereochemical arrangement for identical ligands.
    
    Returns: 'cis', 'trans', 'mer', 'fac', or positional index
    """
    # For two identical ligands
    if len(identical_ligands) == 2:
        angle = compute_angle(donor1, metal_coord, donor2)
        if angle > 150:  # ~180°
            return 'trans'
        else:  # ~90°
            return 'cis'
    
    # For three identical ligands (octahedral)
    if len(identical_ligands) == 3:
        # Compute sum of pairwise angles
        angles = [compute_angle(d1, metal, d2) for d1, d2 in combinations(donors, 2)]
        avg_angle = np.mean(angles)
        if avg_angle < 100:  # All ~90°
            return 'fac'
        else:  # One ~180°
            return 'mer'
    
    # General case: use canonical axis projection
    return assign_by_canonical_axis(metal_coord, donor_coords)
```

**Canonical Axis Definition:**
1. **v1**: Principal axis = eigenvector of largest eigenvalue of donor covariance matrix
2. **v2**: Secondary axis = perpendicular to v1, in plane of v1 and centroid-to-first-donor
3. **v3**: Tertiary axis = v1 × v2

**Tolerance for floating-point stability:**
- Angle comparison: ±5°
- Distance comparison: ±0.01 Å

### S4.4 Determinism Verification

We verified canonicalization determinism on 10,000 structures:

```
Test: 100 random rotations + translations per structure
Result: 100% identical CoordRep strings
Hash collision rate: 0 / 1,000,000 comparisons
```

---

## S5. Tokenizer Specification

### S5.1 Token Types

| Type | Pattern | Examples | Count |
|------|---------|----------|-------|
| Metal+CN | `[Metal:{M}\|CN:{n}]` | `[Metal:Fe\|CN:6]`, `[Metal:Pt\|CN:4]` | 191 |
| Shape | `<ShapeBest:{L}\|...>` | `<ShapeBest:Oh\|Class:dist\|Delta:3\|V:4.99,12.35>` | varies |
| Stereo | `{trans:...}`, `{cis:...}` | `{trans:L1:N:1--L3:N:1}` | varies |
| Ligand Header | `\|L{n}=` | `\|L1=`, `\|L2=` | 1-14 |
| Donor | `:{atom}:{idx}` | `:N:1`, `:O:2`, `:Cl:1` | varies |
| SMILES chars | Single characters | `C`, `(`, `)`, `=`, `#`, `@` | ~50 |
| Special | `[CLS]`, `[SEP]`, `[MASK]`, `[PAD]`, `[UNK]` | - | 5 |

### S5.2 Token Grammar (BNF)

```bnf
<coordrep> ::= <metal_token> <shape_token> <stereo_tokens>? <ligand_list>

<metal_token> ::= "[Metal:" <element> "|" <oxidation>? <d_config>? "CN:" <number> "]"

<shape_token> ::= "<ShapeBest:" <shape_label> "|Class:" <class> "|Delta:" <number> 
                  "|V:" <number> "," <number> ">"

<stereo_tokens> ::= <stereo_constraint>+
<stereo_constraint> ::= "{" <stereo_type> ":" <donor_pair> "}"
<stereo_type> ::= "trans" | "cis" | "mer" | "fac"
<donor_pair> ::= <donor_ref> "--" <donor_ref>
<donor_ref> ::= "L" <number> ":" <element> ":" <number>

<ligand_list> ::= <ligand>+
<ligand> ::= "|L" <number> "=" <smiles> <donor_markers>
<donor_markers> ::= (":" <element> ":" <number>)+
```

### S5.3 Vocabulary

**Total vocabulary size:** 10,000 tokens (padded for efficiency)
**Active tokens:** ~657

| Category | Count | Examples |
|----------|-------|----------|
| Special tokens | 5 | `[PAD]`, `[UNK]`, `[CLS]`, `[SEP]`, `[MASK]` |
| Metal+CN composite | 191 | `[Metal:Fe\|CN:6]`, `[Metal:Pt\|CN:4]` |
| Structure prefixes | ~50 | `<ShapeBest:`, `{trans:`, `\|L1=` |
| SMILES characters | ~50 | `C`, `N`, `O`, `(`, `)`, `=`, `#`, `[`, `]` |
| Donor markers | ~100 | `:N:1`, `:O:2`, `:Cl:1` |
| Other | ~261 | Numbers, separators, rare tokens |

**Full vocabulary file:** `tokenizers/coordrep_vocab.txt`

### S5.4 Composite Token Design Rationale

**Why `[Metal:Fe|CN:6]` instead of separate `[Fe]` and `[CN:6]`?**

1. **Reduced sequence length**: Single token vs. multiple tokens
2. **Enforced consistency**: Metal and CN always co-occur
3. **Better gradient flow**: Single embedding for metal context

**Limitation**: Makes independent metal/CN prediction difficult (see Section 2.4 of main text).

### S5.5 Known Failure Modes

| Issue | Cause | Impact | Mitigation |
|-------|-------|--------|------------|
| Metal prediction = 0% | 191 composite tokens, sparse samples | Cannot probe metal inference directly | Use element extraction for evaluation |
| Shape token fragmented | `<ShapeBest:Oh\|...>` → `<Shape:>` + `>` | Shape label not recoverable | Noted in main text; future tokenizer redesign |
| Long SMILES truncation | Max sequence length = 512 | ~2% structures truncated | Accept or increase length |

### S5.6 Future Improvements

1. **Factorized tokens**: Separate `[Metal:Fe]` and `[CN:6]`
2. **Constrained decoding**: Enforce grammar during generation
3. **Subword tokenization**: BPE for SMILES portions

---

## S6. Masking Protocol and Task Definitions

### S6.1 Tool A: Donor Atom Masking

**Objective:** Predict donor atom identity given metal, geometry, and ligand context.

**Masking rule:**
```python
def mask_donor_atoms(tokens, mask_ratio=0.15):
    """
    Mask donor atom tokens (e.g., :N:1, :O:2).
    
    Args:
        tokens: List of tokens
        mask_ratio: Fraction of donor tokens to mask
    
    Returns:
        masked_tokens, mask_positions, targets
    """
    donor_pattern = re.compile(r':[A-Z][a-z]?:\d+')
    donor_positions = [i for i, t in enumerate(tokens) if donor_pattern.match(t)]
    
    n_mask = max(1, int(len(donor_positions) * mask_ratio))
    mask_positions = random.sample(donor_positions, n_mask)
    
    masked_tokens = tokens.copy()
    targets = []
    for pos in mask_positions:
        targets.append(tokens[pos])
        masked_tokens[pos] = '[MASK]'
    
    return masked_tokens, mask_positions, targets
```

### S6.2 Tool B: Structure Token Masking

**Objective:** Recover structure-defining tokens (brackets, delimiters, ligand indices).

**Masking rule:**
```python
STRUCTURE_TOKENS = {'L', '|', ':', '(', ')', '[', ']', '{', '}', ';', ','}

def mask_structure_tokens(tokens, mask_ratio=0.15):
    """Mask structure-defining tokens."""
    struct_positions = [i for i, t in enumerate(tokens) 
                        if t in STRUCTURE_TOKENS or t.startswith('L')]
    
    n_mask = max(1, int(len(struct_positions) * mask_ratio))
    mask_positions = random.sample(struct_positions, n_mask)
    
    # ... (similar to above)
```

### S6.3 Dynamic Masking

**Implementation:**
- Masking is performed **dynamically** at each epoch
- Different positions are masked each time the same sequence is seen
- Random seed: `epoch * 10000 + sample_index` for reproducibility

```python
def dynamic_mask(tokens, epoch, sample_idx, task='donor'):
    """Dynamic masking with reproducible randomness."""
    rng = np.random.RandomState(epoch * 10000 + sample_idx)
    
    if task == 'donor':
        return mask_donor_atoms(tokens, rng=rng)
    elif task == 'structure':
        return mask_structure_tokens(tokens, rng=rng)
```

### S6.4 Task Distribution During Training

| Task | Probability | Purpose |
|------|-------------|---------|
| Random token mask | 60% | General language modeling |
| Donor atom mask | 25% | Tool A pretraining |
| Structure token mask | 15% | Tool B pretraining |

---

## S7. Synthetic Corruption Suite (Tool B)

### S7.1 Corruption Types

| Type | Description | Probability | Example |
|------|-------------|-------------|---------|
| `missing_bracket` | Delete one bracket randomly | 25% | `[Metal:Fe\|CN:6]` → `Metal:Fe\|CN:6]` |
| `missing_charge` | Delete charge marker | 20% | `+1` → (deleted) |
| `delimiter_swap` | Replace delimiter with another | 25% | `\|L1=` → `;L1=` |
| `truncation` | Remove 5-20% of end tokens | 20% | `...\|L3=CCN:N:1` → `...\|L3=CC` |
| `bracket_swap` | Swap adjacent brackets | 10% | `[(` → `([` |

### S7.2 Generation Parameters

```python
CORRUPTION_CONFIG = {
    'missing_bracket': {
        'prob': 0.25,
        'target_brackets': ['(', ')', '[', ']', '{', '}'],
    },
    'missing_charge': {
        'prob': 0.20,
        'target_pattern': r'[\+\-]\d',
    },
    'delimiter_swap': {
        'prob': 0.25,
        'swap_map': {'|': ';', ';': '|', ',': ';'},
    },
    'truncation': {
        'prob': 0.20,
        'min_remove': 5,
        'max_remove_frac': 0.20,
    },
    'bracket_swap': {
        'prob': 0.10,
        'max_distance': 10,  # Only swap if within 10 characters
    },
}
```

### S7.3 Dataset Statistics

| Corruption Type | N Samples | Valid Before | Valid After Corruption |
|-----------------|-----------|--------------|------------------------|
| missing_bracket | 400 | 95.5% | 0.0% |
| missing_charge | 400 | 96.8% | 96.8% |
| delimiter_swap | 400 | 96.8% | 96.8% |
| truncation | 400 | 96.0% | 36.8% |
| bracket_swap | 400 | 95.3% | 24.3% |
| **Total** | **2,000** | **96.1%** | **50.9%** |

### S7.4 Evaluation Metrics

1. **Valid Parse Rate**: Fraction of sequences passing syntax validation
2. **Exact Match Rate**: Fraction recovering original sequence exactly
3. **Edit Distance**: Levenshtein distance between repaired and original (normalized)

### S7.5 Data File

**Location:** `outputs/fig5/synth_corrupted.jsonl`

**Format:**
```json
{
  "id": "synth_missing_bracket_001",
  "coordrep_clean": "[Metal:Fe|CN:6]<ShapeBest:Oh|...",
  "coordrep_corrupted": "Metal:Fe|CN:6]<ShapeBest:Oh|...",
  "corruption_type": "missing_bracket",
  "was_valid_before_corruption": true,
  "is_valid_after_corruption": false
}
```

---

## S8. Baselines

### S8.1 Random Baseline

**Definition:** Uniform random sampling from candidate token set.

```python
def random_baseline(candidate_set, k=5):
    """
    Random baseline: uniform distribution over candidates.
    
    Args:
        candidate_set: Set of valid tokens (e.g., DONOR_ATOMS)
        k: Number of predictions
    
    Returns:
        List of (token, probability) tuples
    """
    prob = 1.0 / len(candidate_set)
    shuffled = random.sample(candidate_set, min(k, len(candidate_set)))
    return [(tok, prob) for tok in shuffled]
```

**For donor prediction:** `candidate_set = {'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I', ...}` (15 elements)
- Random Top-1: 6.67%
- Random Top-5: 33.3%

### S8.2 Global Frequency Baseline

**Definition:** Predict based on overall token frequency in training set.

```python
def global_frequency_baseline(training_samples, target_field='donor_atom'):
    """
    Build frequency table from training data.
    """
    freq = Counter()
    for sample in training_samples:
        freq[sample[target_field]] += 1
    
    total = sum(freq.values())
    ranked = [(tok, count / total) for tok, count in freq.most_common()]
    return ranked
```

**Training set donor distribution:**
| Donor | Frequency | Cumulative |
|-------|-----------|------------|
| C | 80.2% | 80.2% |
| N | 7.9% | 88.1% |
| O | 7.9% | 96.0% |
| P | 1.5% | 97.5% |
| F | 1.3% | 98.8% |
| S | 1.1% | 99.9% |
| I | 0.1% | 100% |

### S8.3 Conditional Frequency Baseline (CondFreq)

**Definition:** Predict based on frequency conditioned on (metal, CN).

```python
def conditional_frequency_baseline(training_samples, condition_fields=['metal', 'cn']):
    """
    Build conditional frequency table: P(donor | metal, CN)
    """
    cond_freq = defaultdict(Counter)
    
    for sample in training_samples:
        key = (sample['metal'], sample['cn'])
        cond_freq[key][sample['donor_atom']] += 1
    
    # Convert to ranked probabilities
    cond_ranked = {}
    for key, counter in cond_freq.items():
        total = sum(counter.values())
        cond_ranked[key] = [(tok, c/total) for tok, c in counter.most_common()]
    
    return cond_ranked
```

**Important:** Frequency tables are built **only from training set** to prevent data leakage.

### S8.4 ΔTop-1 Definition

```
ΔTop-1 = Top-1_model - Top-1_condfreq
```

**Filtering criteria:**
- Only (metal, CN) groups with ≥10 samples
- Results: 179 groups analyzed

**Distribution:**
- Model wins (Δ > 0): 120 groups (67.0%)
- CondFreq wins (Δ < 0): 30 groups (16.8%)
- Ties (Δ = 0): 29 groups (16.2%)
- Mean ΔTop-1: +4.26 percentage points

---

## S9. Training Details

### S9.1 Hardware

| Component | Specification |
|-----------|---------------|
| GPU | NVIDIA A100 80GB |
| CPU | AMD EPYC 7742 64-Core |
| RAM | 512 GB |
| Storage | NVMe SSD |

### S9.2 Hyperparameters

| Hyperparameter | Value | Notes |
|----------------|-------|-------|
| Architecture | Transformer Encoder | 6 layers |
| Hidden dimension | 384 | |
| Attention heads | 6 | |
| Feed-forward dim | 1536 | 4× hidden |
| Max sequence length | 512 | |
| Vocabulary size | 10,000 | Padded |
| | | |
| Batch size | 256 | |
| Learning rate | 1e-4 | |
| LR scheduler | Linear warmup + cosine decay | |
| Warmup steps | 1,000 | |
| Weight decay | 0.01 | |
| Dropout | 0.1 | |
| | | |
| Epochs | 50 | |
| Early stopping | 5 epochs | Based on val loss |
| Gradient clipping | 1.0 | |

### S9.3 Training Dynamics

**Training curve:**

| Epoch | Train Loss | Val Loss | Val Perplexity | Syntax Validity |
|-------|------------|----------|----------------|-----------------|
| 1 | 4.82 | 4.21 | 67.3 | 12.4% |
| 5 | 2.34 | 2.18 | 8.85 | 78.6% |
| 10 | 1.56 | 1.52 | 4.57 | 95.2% |
| 20 | 1.12 | 1.14 | 3.13 | 99.1% |
| 30 | 0.89 | 0.95 | 2.59 | 99.7% |
| 40 | 0.74 | 0.87 | 2.39 | 99.9% |
| 50 | 0.65 | 0.85 | 2.34 | 99.9% |

**Best checkpoint:** Epoch 47 (lowest val loss = 0.84)

### S9.4 Syntax Validity During Training

Syntax validity (fraction of generated sequences passing CoordRep parser) reaches:
- **>99.9%** by step 2,800
- Remains stable thereafter

### S9.5 Reproducibility Seeds

| Component | Seed |
|-----------|------|
| Data split | 42 |
| Model initialization | 42 |
| Training shuffle | 42 |
| Dynamic masking | epoch × 10000 + sample_idx |

---

## S10. Reproducibility Package

To maximize transparency and enable independent verification of all claims (Tools A/B, baselines, and Fig. 5 panels), we provide a self-contained reproducibility package organized around three principles: (i) fixed interfaces for canonicalization and tokenization, (ii) scriptable end-to-end reproduction from checkpoints, and (iii) a CSV-only reproduction path that regenerates all plots without rerunning inference.

### S10.1 Package Scope and Directory Layout

The repository is structured to separate (a) the core CoordRep/MLM implementation, (b) evaluation utilities for Tool A/B and baselines, and (c) immutable figure artifacts.

- Core model + tokenizer reside under `libcoordrep/brain/` (`model.py`, `tokenizer.py`).
- Tool implementations and validators are packaged under `libcoordrep/coordrep_tools/` (inference API, Tool A, Tool B, baselines, syntax validator).
- Repro scripts live in `libcoordrep/scripts/`, including Tool A/B evaluation entrypoints and a Fig. 5 plotting script.
- Figure-grade outputs are written to `outputs/fig5/`, with a dedicated `reproducible_data/` folder containing the CSVs needed to regenerate Fig. 5 without recomputing model predictions.
- Preprocessed splits are provided as JSONL files (`data/coordrep/train.jsonl`, `val.jsonl`, `test.jsonl`) and the released checkpoint is stored under `checkpoints/pretrain_v3/`.

```
coordrep/
├── README.md                          # Quick-start guide
├── pyproject.toml                     # Package metadata & dependencies
├── LICENSE                            # MIT License
│
├── libcoordrep/
│   ├── coordrep/                      # Core representation library
│   │   ├── canonical/                 #   Canonicalization engine
│   │   ├── geometry/                  #   CShM computation
│   │   ├── graph/                     #   Molecular graph & ligand extraction
│   │   ├── io/                        #   File readers (CIF, tmQM XYZ)
│   │   ├── serialize/                 #   String serialization
│   │   ├── features/                  #   ML feature extraction
│   │   └── validate/                  #   Consistency checks
│   │
│   ├── brain/                         # Masked Language Model
│   │   ├── model.py                   #   Transformer encoder architecture
│   │   ├── tokenizer.py              #   CoordRep tokenizer
│   │   ├── dataloader.py             #   Training data pipeline
│   │   └── masking.py                #   Dynamic masking strategies
│   │
│   ├── coordrep_tools/               # Downstream tool implementations
│   │   ├── infer.py                   #   MLM inference API
│   │   ├── tool_a_donor.py            #   Tool A: Donor-atom prediction
│   │   ├── tool_b_repair.py           #   Tool B: Structure repair
│   │   ├── baselines.py               #   Random / Freq / CondFreq baselines
│   │   └── validate.py                #   CoordRep syntax validation
│   │
│   ├── scripts/                       # Reproduction entry points
│   │   ├── run_tool_a_eval.py         #   Tool A evaluation
│   │   ├── run_tool_a_from_preds.py   #   Tool A from cached predictions
│   │   ├── run_tool_b_eval.py         #   Tool B evaluation
│   │   ├── run_tool_b_from_preds.py   #   Tool B from cached predictions
│   │   └── plot_fig5_panels.py        #   Generate all Fig. 5 panels
│   │
│   ├── tests/                         # Invariance and validation tests
│   │
│   └── outputs/
│       ├── fig5/
│       │   └── reproducible_data/     # CSV-only Fig. 5 reproduction
│       │       ├── plot_fig5_from_csv.py
│       │       ├── panel_b_*.csv
│       │       ├── panel_c_*.csv
│       │       ├── panel_d_*.csv
│       │       └── panel_e_*.csv
│       └── SI/                        # Supplementary Information artifacts
│
├── data/
│   └── coordrep/
│       ├── train.jsonl                # Training split (79,964 samples)
│       ├── val.jsonl                  # Validation split (9,996 samples)
│       └── test.jsonl                 # Test split (9,995 samples)
│
└── checkpoints/
    └── pretrain_v3/
        ├── best_model.pt              # Released model checkpoint
        ├── config.json                # Training configuration
        └── tokenizer.json             # Tokenizer vocabulary
```

### S10.2 One-Command Reproduction of Tools A/B and Fig. 5

We provide canonical reproduction commands that (1) install the environment, (2) run Tool A evaluation, (3) run Tool B evaluation, and (4) regenerate all Fig. 5 panels.

```bash
# Step 1: Install environment
pip install -e ".[all]"
# Or with conda:
# conda env create -f environment.yml && conda activate coordrep

# Step 2: Run Tool A evaluation (donor-atom prediction)
python -m libcoordrep.scripts.run_tool_a_from_preds

# Step 3: Run Tool B evaluation (structure repair)
python -m libcoordrep.scripts.run_tool_b_from_preds

# Step 4: Generate all Fig. 5 panels
python -m libcoordrep.scripts.plot_fig5_panels
```

**CSV-only pathway (no GPU required):** A key feature is the CSV-only reproduction path. The directory `outputs/fig5/reproducible_data/` contains pre-computed evaluation CSVs for every Fig. 5 panel. Reviewers can regenerate the final plots deterministically without rerunning any model inference:

```bash
cd libcoordrep/outputs/fig5/reproducible_data
python plot_fig5_from_csv.py
# → Produces plots/panel_b.pdf, plots/panel_c.pdf, plots/panel_d.pdf, plots/panel_e.pdf
```

### S10.3 Deterministic Seeds and Run-to-Run Stability

To ensure that splits, shuffling, and dynamic masking are reproducible, we fix all random seeds as follows:

| Component | Seed | Notes |
|-----------|------|-------|
| Data split (train/val/test) | 42 | `sklearn.train_test_split(random_state=42)` |
| Model weight initialization | 42 | `torch.manual_seed(42)` |
| Training data shuffle | 42 | Per-epoch shuffle seed |
| Dynamic masking | `epoch × 10000 + sample_idx` | Ensures diverse but reproducible masks |

The dynamic masking seed design (`epoch × 10000 + sample_idx`) ensures that "dynamic" masking remains diverse across epochs—each epoch exposes different masked positions for the same sequence—while being exactly reproducible for a given run configuration.

### S10.4 Dependencies

```
# Core
python>=3.9
pytorch>=2.0.0
numpy>=1.24.0
pandas>=2.0.0
scikit-learn>=1.2.0

# Chemistry
rdkit>=2023.03.1
gemmi>=0.6.4

# Visualization
matplotlib>=3.7.0
seaborn>=0.12.0

# Utilities
tqdm>=4.65.0
jsonlines>=3.1.0
```

### S10.5 Data Availability

| Dataset | Availability | Access |
|---------|--------------|--------|
| tmQM | Public (CC BY 4.0) | https://github.com/bbskjelstad/tmqm |
| COD | Public Domain | https://www.crystallography.net/cod/ |
| CoordRep processed splits | Public | Included in repository (`data/coordrep/`) |
| Model checkpoint | Public | Included in repository (`checkpoints/pretrain_v3/`) |
| Fig. 5 CSV data | Public | Included in repository (`outputs/fig5/reproducible_data/`) |

### S10.6 Computational Requirements

| Task | Time | GPU Memory | GPU Required? |
|------|------|------------|---------------|
| Full pretraining (50 epochs) | ~8 hours | ~40 GB | Yes (A100 recommended) |
| Tool A evaluation | ~5 minutes | ~8 GB | Yes |
| Tool B evaluation | ~10 minutes | ~8 GB | Yes |
| Fig. 5 from CSV (plot only) | ~30 seconds | — | No (CPU only) |
| SI figure generation | ~2 minutes | — | No (CPU only) |

---

## Supplementary Tables

### Table S1. CShM Reference Polyhedra Coordinates
(See Section S3.4)

### Table S2. Full Donor Priority List
(See Section S4.2)

### Table S3. Complete Vocabulary
(Available as `tokenizers/coordrep_vocab.txt`)

### Table S4. Per-Metal Evaluation Results
(Available as `outputs/fig5/reproducible_data/panel_b_delta_by_metal_cn.csv`)

### Table S5. Confidence-Coverage Data Points
(Available as `outputs/fig5/reproducible_data/panel_c_confidence_coverage.csv`)

---

## Supplementary Figures

### Figure S1. Training Curves
(Loss, perplexity, and syntax validity vs. epoch)

### Figure S2. Data Distribution
(Metal, CN, and donor atom distributions in train/val/test)

### Figure S3. CShM Distribution
(For CN=4-8, comparing tmQM and COD)

### Figure S4. Extended ΔTop-1 Analysis
(All 179 (metal, CN) groups)

### Figure S5. Correction Case Studies
(30 examples where model corrects CondFreq)

---

*End of Supplementary Information*
