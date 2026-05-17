    # Supplementary Note: CoordRep-v2 Beta Extensions

    ## Scope Statement

    CoordRep v1 is restricted to **mononuclear, atom-resolved η1 coordination
    snapshots**.  The v2-beta extensions described here demonstrate chemically
    motivated grammar extensibility toward multinuclear and haptic/π coordination.

    These prototype tests demonstrate chemically motivated grammar extensibility,
    not full production support for all multinuclear, periodic, or organometallic
    structures.

    ## Metal-Centered Record Graph (CoordRep-Multi-v2beta)

    For multinuclear complexes, CoordRep-v2beta introduces a **metal-centered
    record graph** consisting of:

    - **Metals block**: each metal center with element, oxidation state, d-count,
      CNsite, CNatom, and local shape.
    - **MetalGraph block**: explicit metal–metal edges with relation type
      (bridged, direct_MM_bond, contact), distance, and bond assignment.
    - **LocalSphere block**: per-metal coordination environment.
    - **Bridges block**: bridging ligands/donors with target metals and μ-value.
    - **Ligands block**: each ligand appears once regardless of bridging.
    - **ID block**: hierarchical L0–L3 identity keys plus per-metal local IDs.

    ### Example: Cu2(mu2-Cl)2(phen)2 chloro-bridged

    ```
    [Metals:
  M1=[Cu|ox:+2|d:d9|CNsite:4|CNatom:4];
  M2=[Cu|ox:+2|d:d9|CNsite:4|CNatom:4]]
[MetalGraph:
  {M1--M2|relation:bridged|dMM:3.10|MM_bond:no}]
[LocalSphere:
  M1=<ShapeBest:SP|Class:|Delta:0|V:>{S1,S2,S3,S4};
  M2=<ShapeBest:SP|Class:|Delta:0|V:>{S1,S2,S5,S6}]
[Bridges:
  {site:S1|ligand:L1|donor:ClCl1|targets:M1,M2|mu:2|mode:mu2-Cl};
  {site:S2|ligand:L1|donor:ClCl2|targets:M1,M2|mu:2|mode:mu2-Cl}]
[Ligands:
  L1=[Cl-];
  L2=c1ccnc2cccnc12;
  L3=c1ccnc2cccnc12]
[ID:
  L0_GlobalState=e9da24b45036cf29;
  L1_GlobalShape=d9678341f25fb8df;
  L2_MetalGraphTopo=31719ab392aa9df3;
  L3_Connectivity=cc65aff80fca431c;
  LocalID[M1]=a868d5b91e42;
  LocalID[M2]=bc796eed5405]
    ```

    ## Coordination-Site Object (CoordRep-Haptic-v2beta)

    For π/haptic ligands, CoordRep-v2beta introduces a **coordination-site
    object** that represents:

    - **atom**: conventional η1 donor
    - **pi_fragment**: contiguous conjugated atom set (η2–η6)
    - **centroid**: geometric center of a π system

    Each site carries η-value, μ-value, target metals, and canonical atom ordering.

    ### Example: Ferrocene Fe(eta5-Cp)2

    ```
    [Metal:
  M1=[Fe|ox:+2|d:d6|CNsite:2|eta_sum:10]]
[Sites:
  S1={type:pi_fragment|ligand:L1|atoms:CC1,CC2,CC3,CC4,CC5|eta:5|centroid:Cp1|mode:Cp};
  S2={type:pi_fragment|ligand:L2|atoms:CC10,CC6,CC7,CC8,CC9|eta:5|centroid:Cp2|mode:Cp}]
[Hapticity:
  {S1->M1|eta:5|mode:eta5-Cp};
  {S2->M1|eta:5|mode:eta5-Cp}]
[Ligands:
  L1=[cH-]1cccc1;
  L2=[cH-]1cccc1]
[Geometry:
  CNsite=2;
  eta_sum=10;
  site_level_shape=linear]
[ID:
  L0_HapticState=cf3d47e88915d796;
  L1_HapticShape=cff1b68b0dc711a0;
  L2_SiteTopo=37b2504f427685a9;
  L3_Connectivity=3ec24b7f02a0c93a]
    ```

    ## Validation Metrics

    | Metric | Multinuclear | Haptic |
    |---|---|---|
    | Curated cases | 55 | 33 |
    | Parse + roundtrip rate | 100.0% | 100.0% |
    | Invariance rate | 100.0% | 100.0% |
    | Random audit (estimated) | 72% | 68% |
    | Estimated coverage gain | +228.1% over v1 baseline |

    ## Limitations and Failure Modes

    - Ambiguous oxidation states in mixed-valence systems
    - Complex bridging topologies (μ4+, rare connectivity)
    - Mixed haptic and bridging on same ligand
    - Polymeric / extended structures beyond molecular clusters
    - Unusual hapticity (η7+, η8 in actinides)
    - Lanthanide/actinide edge cases

    ## Statement

    This is **not** part of the main v1 benchmark.  CoordRep v1 claims and
    validation remain restricted to mononuclear η1 coordination records.
    These extensions are provided as supplementary evidence of grammatical
    extensibility.
