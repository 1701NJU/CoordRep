# Comparison scope and protocol

## Input routes and reference independence

The molecular route reads the original optimized XYZ inputs associated with
the published 44-Fe comparison. The native-source route reads the fixed 40
crystallographic cases through a local CCDC environment. The separate
conditional-CIF route reads the deposited source files for BOWCAV and
MAXNUY01, retaining compatible disorder alternatives and physical periodic
donor images. Its 18 emitted local candidates include the 16 qualified
reference alternatives and the two unscored third-orientation candidates.

Literature comparison targets are assembled from explicit donor assignments
and drawings in the original papers and their supporting material. A target
is scored as a complete set only when the reference specifies the complete
primary non-H donor set for the relevant structure/state. A coordination
number alone, a partially labelled drawing, or an automatically generated MOL
file does not qualify as an independently complete donor-set reference.
Experimental assignments are not silently transferred to a different optimized
geometry. Source-derived execution checks and independent chemical-reference
comparisons are therefore kept separate.

For XYZ inputs the frozen contact policy includes the implemented donor-element
domain, the 1.3 times covalent-radius-sum distance cap and the simultaneous
all-heavy-atom Gabriel screening. No per-case threshold fitting is used. The
native crystallographic route retains the source bond/contact assignments; it
is not the same input construction as the XYZ route.

## Table S16: optimized molecular donor sets

All 21 listed targets have complete primary-heavy-donor references for the
tested optimized state. Original case numbers and optimized-coordinate
references follow S1. All 103 assigned donors are retained, with no additional
donor in these 21 complete sets. The full execution register retains all
44 structures and 47 Fe centres, irrespective of reference qualification.
`molecular_public_reference_indices.json` supplies the zero-based Fe and donor
indices for these 21 external published XYZ inputs, with input SHA-256 and
the original SI page/figure. The element-specific assignments follow the
complete drawn donor inventories, not a distance cutoff or an automatically
generated MOL file. No native-CSD atom mappings are included.

Case 28 is reported separately rather than counted as a failed program or a
successful same-state chemical match. The experimental Fe complex in S2
provides an FeN2S2, CN=4 assignment; the tested optimized structure presents a
CN=5 contact candidate under the frozen XYZ policy. The states have different
Fe-N-O geometries. The chemical assignment of the additional Fe-O contact in
the optimized state is unresolved by the experimental reference alone. This
accounts for the 40/41 single-Fe CN agreement without changing the
21/21 complete optimized-state donor-set result.

## Table S17: crystallographic donor sets

Twelve entries provide 19 complete local reference sets. All 107 independently
assigned nonmetal donor contacts are retained. Seventeen sets are exact.
The two differences are the two Cu sites of NUXNIG (S3): the paper identifies
CuN2O2 primary spheres, whereas the source-supported record also contains two
longer Cu-O contacts per Cu. The four additional contacts are retained and
reported; they are not relabelled as exact matches or treated as missing
primary bonds. This distinction separates a primary chemical sphere from
additional recorded contacts.

Six directed Au-Tl relations in TEGVAD (S4) are tested separately from the
107 nonmetal donor contacts. Directed incidences are not a count of six
independent undirected Au-Tl bonds. Shared-donor, haptic and metal-metal
descriptions are retained as distinct relation types.

## Table S18: disorder and periodic images

BOWCAV (S5) supplies six split-site Zn pairs, with two local ZnN2I2 alternatives
per pair: 12 alternatives and 48 donor images. MAXNUY01 (S6) supplies two Na
sites under each of the first two BF4 orientations: four NaF6 alternatives
and 24 donor images. The 16 complete local sets and all 72 physical donor
images agree, including donor images reached through a nonzero translation.
A donor image denotes a physical local contact, not a unique wrapped atom
orbit; the count 72 does not mean 72 nonzero-translation edges.

The third BF4 orientation generates two additional local Na candidates that
remain in the output. It lacks the equally complete independently assigned
local donor reference required for scoring here. Its candidates are not
removed, not counted among the 16 qualified references, and not presented as
reference-validated NaF6 alternatives. These are conditional local alternatives,
not a determination of the globally correlated crystal disorder ensemble.

## Reproducibility and stability endpoints

Fresh-process execution reproduces all 44 molecular outputs (47 local Fe
records), the 40 native crystallographic outputs and the two conditional-CIF
outputs. The latter contain 12 and six local records, respectively. Reproduction
of an output is not by itself chemical validation of its donor assignments.

The molecular stress tests contain 470 centre-level trials across the 47 Fe
centres, plus 2,000 additional centre-level trials for the two Fe centres of
case 13. All preserve their tested donor-set, graph and local-record endpoints.
These tests establish local-record stability, not chemical valence correctness
or a whole-assembly identity for the two-Fe structure.

For the two conditional-CIF cases, 24 actual adapter executions yield 216
local-record comparisons; all pass. Executions and local comparisons are
different units. The test family includes coordinate-frame changes, donor/atom
reordering, periodic re-expression, and removal or poisoning of deposited
geometry hints. The report does not claim that these 24 executions are
24 different compounds or additional independent chemical determinations.

## How to read the register

`processing_and_fresh_replay=COMPLETE` means an input was encoded and its
archived output reproduced. `complete_local_reference_qualified` identifies a
complete scored reference endpoint. The unconditional and conditional columns
must be interpreted separately; a case can have useful partial literature
information without a complete qualified set. Every case remains visible.
No whole-CSD chemical-accuracy rate, unseen-test claim or human-blind expert
study is inferred from these results.
