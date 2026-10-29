# CIJWUO — Stereo Semantic Consistency (Fig. 5C, Case D)

## What this case demonstrates

CIJWUO is a Cu(II) CN=5 (TBP, ideal) complex with a trans constraint:
`{trans:L1:N:1--L2:N:1}`. A stereo decoy is created by flipping trans→cis.

**Crucially**, the decoy passes grammar validation — it is syntactically valid.
Only the CoordRep-Ranker (a learned model) can detect that the cis claim is
semantically inconsistent with the actual geometry. This demonstrates that
CoordRep's structured fields enable **semantic consistency checking** beyond
what syntax rules can provide.

## Rendering guidance

- Use `CIJWUO_stereo_view.xyz` (metal + 5 donors) for a clean view.
- Use `CIJWUO_first_sphere.xyz` for a slightly richer ball-and-stick.
- Highlight Cu; label the two N donors involved in the trans/cis pair.
- Optionally draw an arrow or line between the trans pair to show the relationship.
- The full molecule (73 atoms) is too crowded for a small panel.

## Key numbers

| Metric | Value |
|--------|-------|
| Real stereo | `{trans:L1:N:1--L2:N:1}` |
| Decoy stereo | `{cis:L1:N:1--L2:N:1}` |
| Grammar-valid decoy | Yes |
| Ranker correct | 5/5 |
| CIJWUO margin | +0.685 |
| Margin range (all 5) | 0.10 – 0.68 |
