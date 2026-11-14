#!/usr/bin/env python3
"""
Build manuscript_revised_highlighted_major_changes.docx

Strategy:
- Copy the revised manuscript as-is for the clean version.
- For the highlighted version, add yellow highlighting to paragraphs
  containing major revisions per the user's specification, and add
  Word comments at section headings for substantially revised sections.

Output:
  revision_results/manuscript_revised_clean.docx
  revision_results/manuscript_revised_highlighted_major_changes.docx
"""

import copy
import re
import shutil
from pathlib import Path
from docx import Document
from docx.shared import RGBColor, Pt
from docx.oxml.ns import qn, nsmap
from docx.oxml import OxmlElement
import datetime

BASE = Path(__file__).resolve().parent.parent / "revision_results"
REVISED = BASE / "manuscript - Revise.docx"
OUT_CLEAN = BASE / "manuscript_revised_clean.docx"
OUT_HIGHLIGHTED = BASE / "manuscript_revised_highlighted_major_changes.docx"

# ═══════════════════════════════════════════════════════════════════
# 1. Clean copy
# ═══════════════════════════════════════════════════════════════════
shutil.copy2(REVISED, OUT_CLEAN)
print(f"✓ Clean copy: {OUT_CLEAN.name}")

# ═══════════════════════════════════════════════════════════════════
# 2. Highlighted version
# ═══════════════════════════════════════════════════════════════════
doc = Document(REVISED)

# --- Helper: apply yellow highlight to all runs in a paragraph ---
YELLOW = 7  # WdColorIndex.wdYellow


def highlight_paragraph(para):
    """Add yellow highlight to every run in a paragraph."""
    for run in para.runs:
        run.font.highlight_color = YELLOW


def add_comment_to_paragraph(para, text, author="Cascade"):
    """Add a Word comment to a paragraph.

    This uses the low-level XML approach since python-docx doesn't
    natively support comments creation.
    """
    # We'll use a simplified approach: insert a comment reference
    # via the document's comments part
    pass  # We'll use a different approach below


def add_comment_marker(para, comment_text):
    """Insert comment text as a bracketed note at the start of paragraph
    since python-docx comment support is limited. We use a different
    approach: add a distinctive run at the start with the note."""
    # Actually, let's use the proper OxmlElement approach
    # Create comment in comments part
    pass


# --- Strategy: Instead of complex comment insertion, we'll prepend
#     a clearly marked note run in a distinct style ---

def add_review_note(para, note_text):
    """Add a review note as a bold, italic, red run at the start of the paragraph."""
    run = para.insert_paragraph_before("")
    # Actually let's just add a comment-like annotation
    # Better: use a visible marker that reviewers can find
    pass


# Since Word comments via python-docx are complex, let's use a simpler
# but effective approach: add highlighted "REVIEWER NOTE" annotations
# as separate small paragraphs before the relevant sections.

def insert_reviewer_comment(doc, para_index, comment_text):
    """Insert a reviewer comment as a formatted paragraph before the target."""
    # We can't easily insert paragraphs at arbitrary positions in python-docx
    # without manipulating the XML tree directly
    target_para = doc.paragraphs[para_index]
    # Create a new paragraph element before the target
    new_p = OxmlElement('w:p')
    # Add run with comment text
    new_r = OxmlElement('w:r')
    # Run properties: bold, italic, small font, color
    rPr = OxmlElement('w:rPr')
    b = OxmlElement('w:b')
    rPr.append(b)
    i_elem = OxmlElement('w:i')
    rPr.append(i_elem)
    color = OxmlElement('w:color')
    color.set(qn('w:val'), '0000FF')
    rPr.append(color)
    sz = OxmlElement('w:sz')
    sz.set(qn('w:val'), '18')  # 9pt
    rPr.append(sz)
    new_r.append(rPr)
    # Text
    t = OxmlElement('w:t')
    t.text = f"[{comment_text}]"
    new_r.append(t)
    new_p.append(new_r)
    # Insert before target paragraph's XML element
    target_para._element.addprevious(new_p)


# ═══════════════════════════════════════════════════════════════════
# Identify paragraphs to highlight based on content matching
# ═══════════════════════════════════════════════════════════════════

# Build a text index for quick lookups
para_texts = [(i, p.text.strip()[:200]) for i, p in enumerate(doc.paragraphs)]


def find_para_containing(substring, start=0):
    """Find first paragraph index containing substring (case-insensitive)."""
    sub_lower = substring.lower()
    for i, txt in para_texts:
        if i >= start and sub_lower in txt.lower():
            return i
    return None


def find_para_starting(prefix, start=0):
    """Find first paragraph starting with prefix."""
    prefix_lower = prefix.lower()
    for i, txt in para_texts:
        if i >= start and txt.lower().startswith(prefix_lower):
            return i
    return None


# ── ABSTRACT ──
print("\nHighlighting Abstract...")
abstract_idx = find_para_containing("ABSTRACT:")
if abstract_idx is not None:
    highlight_paragraph(doc.paragraphs[abstract_idx])
    insert_reviewer_comment(doc, abstract_idx,
        "REVIEWER NOTE: Abstract substantially revised — new coverage numbers "
        "(20.3%→81.28%), CoordRep-Rosetta (Top-1=86.7%, AUROC=0.949), "
        "and v2-beta scope statement added.")
    print(f"  ✓ Abstract at para {abstract_idx}")

# ── INTRODUCTION ──
print("Highlighting Introduction...")
intro_idx = find_para_containing("introduction", start=0)

# New prior-work paragraph (Several recent efforts...)
prior_work_idx = find_para_containing("Several recent efforts have begun")
if prior_work_idx is not None:
    highlight_paragraph(doc.paragraphs[prior_work_idx])
    print(f"  ✓ Prior work paragraph at {prior_work_idx}")

# CoordRep addresses this complementary layer
complementary_idx = find_para_containing("CoordRep addresses this complementary layer")
if complementary_idx is not None:
    highlight_paragraph(doc.paragraphs[complementary_idx])
    print(f"  ✓ Complementary layer paragraph at {complementary_idx}")

# Multi-resolution identity separation paragraph
identity_sep_idx = find_para_containing("central feature of CoordRep is that it separates identity")
if identity_sep_idx is not None:
    highlight_paragraph(doc.paragraphs[identity_sep_idx])
    print(f"  ✓ Identity separation paragraph at {identity_sep_idx}")

# ── TABLE 1 ──
print("Highlighting Table 1...")
table1_idx = find_para_containing("Table 1. Operation coverage")
if table1_idx is not None:
    highlight_paragraph(doc.paragraphs[table1_idx])
    insert_reviewer_comment(doc, table1_idx,
        "REVIEWER NOTE: New Table 1 — Operation Coverage under the "
        "Coordination Identity Challenge. Includes new representation "
        "families: TMC structure-to-SMILES, m-SMILES/MetalloGen, T-REX, "
        "atom-attributed graph hashes / WL-type graph keys.")
    # Also highlight table note
    table1_note_idx = find_para_containing("Table note. This table evaluates")
    if table1_note_idx is not None:
        highlight_paragraph(doc.paragraphs[table1_note_idx])
    print(f"  ✓ Table 1 at para {table1_idx}")

# ── TABLE 1 content (table 0) ──
# Highlight new representation families in Table 1
if len(doc.tables) > 0:
    table1_content = doc.tables[0]
    new_families = ["TMC structure-to-SMILES", "m-SMILES", "MetalloGen",
                    "T-REX", "WL-type", "atom-attributed graph"]
    for row in table1_content.rows:
        for cell in row.cells:
            cell_text = cell.text
            if any(fam in cell_text for fam in new_families):
                for para in cell.paragraphs:
                    if para.text.strip():
                        highlight_paragraph(para)

# ── BOX 1 (in tables) ──
print("Highlighting Box 1...")
# Box 1 is in a table. Highlight all content.
for table in doc.tables:
    first_cell_text = table.rows[0].cells[0].text if table.rows else ""
    if "Box 1" in first_cell_text:
        for row in table.rows:
            for cell in row.cells:
                for para in cell.paragraphs:
                    if para.text.strip():
                        highlight_paragraph(para)
        print("  ✓ Box 1 table highlighted")

# Also find any "Box 1" reference paragraph in the body
box1_ref_idx = find_para_containing("Box 1")
if box1_ref_idx is not None:
    # Check if it's a caption/title rather than passing reference
    txt = doc.paragraphs[box1_ref_idx].text
    if "Box 1" in txt and len(txt) < 300:
        insert_reviewer_comment(doc, box1_ref_idx,
            "REVIEWER NOTE: New main-text complete CoordRep records added "
            "(Examples A–E covering mononuclear, multinuclear, and haptic).")

# ── FIGURE 2 SECTION ──
print("Highlighting Figure 2 section...")
# L0 snapshot identity is not InChI
l0_idx = find_para_containing("L0 CoordRep-State is the highest-resolution")
if l0_idx is not None:
    highlight_paragraph(doc.paragraphs[l0_idx])
    print(f"  ✓ L0 identity paragraph at {l0_idx}")

# L1-L3 hierarchy and coarser layers
l1l3_idx = find_para_containing("CoordRep-ID defines coarser identity lay")
if l1l3_idx is not None:
    highlight_paragraph(doc.paragraphs[l1l3_idx])
    print(f"  ✓ L1-L3 hierarchy paragraph at {l1l3_idx}")

# Perturbation/CSD family benchmark
perturb_idx = find_para_containing("trade-off is visible under controlled geometric variation")
if perturb_idx is not None:
    highlight_paragraph(doc.paragraphs[perturb_idx])
    print(f"  ✓ Perturbation benchmark at {perturb_idx}")

family_idx = find_para_containing("evaluated the identity hierarchy on real crystallographic families")
if family_idx is not None:
    highlight_paragraph(doc.paragraphs[family_idx])
    print(f"  ✓ CSD family benchmark at {family_idx}")

# Figure 2 caption (new multi-resolution identity title)
fig2_cap_idx = find_para_containing("Figure 2. CoordRep-ID defines a multi-resolution")
if fig2_cap_idx is not None:
    highlight_paragraph(doc.paragraphs[fig2_cap_idx])
    print(f"  ✓ Figure 2 caption at {fig2_cap_idx}")

# ── FIGURE 3 SECTION ──
print("Highlighting Figure 3 section...")
# CN=6 3d vs 4d/5d spin-state caveat
spin_state_idx = find_para_containing("3d vs 4d/5d")
if spin_state_idx is None:
    spin_state_idx = find_para_containing("spin-state")
if spin_state_idx is None:
    # Try "d-electron" or "Jahn-Teller"
    spin_state_idx = find_para_containing("Jahn")
if spin_state_idx is not None:
    # Only highlight the relevant sentence runs, but for simplicity
    # highlight the whole paragraph if it's the main one
    highlight_paragraph(doc.paragraphs[spin_state_idx])
    print(f"  ✓ Spin-state caveat at {spin_state_idx}")

# Static crystallographic geometry / not dynamics
static_idx = find_para_containing("static crystallographic")
if static_idx is None:
    static_idx = find_para_containing("population")
    # Be careful — check context
    if static_idx is not None:
        txt = doc.paragraphs[static_idx].text
        if "geometry" not in txt.lower() and "crystal" not in txt.lower():
            static_idx = None
if static_idx is not None:
    highlight_paragraph(doc.paragraphs[static_idx])
    print(f"  ✓ Static crystallographic paragraph at {static_idx}")

# Figure 3 caption
fig3_cap_idx = find_para_containing("Figure 3. Continuous CShM geometry fields")
if fig3_cap_idx is not None:
    highlight_paragraph(doc.paragraphs[fig3_cap_idx])
    print(f"  ✓ Figure 3 caption at {fig3_cap_idx}")

# ── FIGURE 4 SECTION ──
print("Highlighting Figure 4 section...")
# Masked-field learning redefined as diagnostic probe
probe_idx = find_para_containing("Masked-field probes separate ligand-intrinsic")
if probe_idx is not None:
    highlight_paragraph(doc.paragraphs[probe_idx])
    insert_reviewer_comment(doc, probe_idx,
        "REVIEWER NOTE: Substantially revised in response to reviewer comments — "
        "masked-field learning reframed as diagnostic probe, not generation claim.")
    print(f"  ✓ Diagnostic probe paragraph at {probe_idx}")

# donor-field attribution controls
donor_attr_idx = find_para_containing("donor-marker recovery")
if donor_attr_idx is None:
    donor_attr_idx = find_para_containing("donor-mark")
if donor_attr_idx is not None:
    highlight_paragraph(doc.paragraphs[donor_attr_idx])
    print(f"  ✓ Donor-field attribution at {donor_attr_idx}")

# Figure 4 caption
fig4_cap_idx = find_para_containing("Figure 4. Masked-field probes separate")
if fig4_cap_idx is not None:
    highlight_paragraph(doc.paragraphs[fig4_cap_idx])
    print(f"  ✓ Figure 4 caption at {fig4_cap_idx}")

# Table 2 tokenizer ablation
table2_idx = find_para_containing("Table 2. Tokenizer ablation")
if table2_idx is not None:
    highlight_paragraph(doc.paragraphs[table2_idx])
    print(f"  ✓ Table 2 at {table2_idx}")
# Table 2 note
table2_note_idx = find_para_containing("Table Note. Values are Top-1 recovery")
if table2_note_idx is not None:
    highlight_paragraph(doc.paragraphs[table2_note_idx])

# Tokenizer ablation results paragraph
tok_abl_idx = find_para_containing("whether these conclusions depend on tokenizer design")
if tok_abl_idx is not None:
    highlight_paragraph(doc.paragraphs[tok_abl_idx])
    print(f"  ✓ Tokenizer ablation results at {tok_abl_idx}")

# Table 3 multidentate
table3_idx = find_para_containing("Table 3. Donor-marker recovery")
if table3_idx is not None:
    highlight_paragraph(doc.paragraphs[table3_idx])
    print(f"  ✓ Table 3 at {table3_idx}")
# Table 3 note
table3_note_idx = find_para_containing("Table note. Donor Top-1 and Top-5")
if table3_note_idx is not None:
    highlight_paragraph(doc.paragraphs[table3_note_idx])

# Multidentate paragraph
multi_idx = find_para_containing("multidentate ligands impose ligand-level constraints")
if multi_idx is not None:
    highlight_paragraph(doc.paragraphs[multi_idx])
    print(f"  ✓ Multidentate paragraph at {multi_idx}")

# Figure 4D graph/equivariant baselines
gnn_idx = find_para_containing("compared CoordRep probes with graph and geometry baselines")
if gnn_idx is not None:
    highlight_paragraph(doc.paragraphs[gnn_idx])
    print(f"  ✓ Graph baselines paragraph at {gnn_idx}")

# ── FIGURE 5 SECTION ──
print("Highlighting Figure 5 section...")
# Full-CSD v1 audit
fullcsd_idx = find_para_containing("Full-CSD application establishes CoordRep")
if fullcsd_idx is not None:
    highlight_paragraph(doc.paragraphs[fullcsd_idx])
    insert_reviewer_comment(doc, fullcsd_idx,
        "REVIEWER NOTE: Substantially revised in response to reviewer comments — "
        "full-CSD audit, v2-beta expansion, boundary geometry, and "
        "CoordRep-Rosetta benchmark added.")
    print(f"  ✓ Full-CSD opening at {fullcsd_idx}")

# v1 pipeline 124,837 and coverage numbers
pipeline_idx = find_para_containing("production v1 pipeline retained 124,837")
if pipeline_idx is not None:
    highlight_paragraph(doc.paragraphs[pipeline_idx])
    print(f"  ✓ Coverage numbers at {pipeline_idx}")

# CN5 ridge / boundary
boundary_idx = find_para_containing("why geometry should be encoded as a continuous, boundary-aware")
if boundary_idx is not None:
    highlight_paragraph(doc.paragraphs[boundary_idx])
    print(f"  ✓ Boundary geometry at {boundary_idx}")

# Casebook and family
casebook_idx = find_para_containing("Casebook and family-level analyses")
if casebook_idx is not None:
    highlight_paragraph(doc.paragraphs[casebook_idx])
    print(f"  ✓ Casebook paragraph at {casebook_idx}")

# CoordRep-Rosetta
rosetta_idx = find_para_containing("whether the expanded CoordRep fields support a downstream operation")
if rosetta_idx is not None:
    highlight_paragraph(doc.paragraphs[rosetta_idx])
    print(f"  ✓ CoordRep-Rosetta paragraph at {rosetta_idx}")

# Full-CSD conclusion paragraph in Fig 5 section
fullcsd_concl_idx = find_para_containing("full-CSD and cross-domain retrieval results position")
if fullcsd_concl_idx is not None:
    highlight_paragraph(doc.paragraphs[fullcsd_concl_idx])
    print(f"  ✓ Full-CSD conclusion at {fullcsd_concl_idx}")

# Figure 5 caption
fig5_cap_idx = find_para_containing("Figure 5. Full-CSD application of CoordRep")
if fig5_cap_idx is not None:
    highlight_paragraph(doc.paragraphs[fig5_cap_idx])
    print(f"  ✓ Figure 5 caption at {fig5_cap_idx}")

# ── CONCLUSIONS ──
print("Highlighting Conclusions...")
conclusions_idx = find_para_containing("CONCLUSIONS")
# v2 beta in conclusions
v2beta_concl_idx = find_para_containing("full-CSD audit further shows that CoordRep is not limited")
if v2beta_concl_idx is not None:
    highlight_paragraph(doc.paragraphs[v2beta_concl_idx])
    print(f"  ✓ v2beta conclusions at {v2beta_concl_idx}")

# Not predictive model / not inverse design
not_pred_idx = find_para_containing("main conclusion is not that neural models replace")
if not_pred_idx is not None:
    highlight_paragraph(doc.paragraphs[not_pred_idx])
    print(f"  ✓ Not-predictive conclusion at {not_pred_idx}")

# Scope limitation paragraph
scope_lim_idx = find_para_containing("Despite these remaining boundaries")
if scope_lim_idx is not None:
    highlight_paragraph(doc.paragraphs[scope_lim_idx])
    print(f"  ✓ Scope limitation at {scope_lim_idx}")

# Key outcome — identity separation
key_outcome_idx = find_para_containing("key outcome of this work is the separation")
if key_outcome_idx is not None:
    highlight_paragraph(doc.paragraphs[key_outcome_idx])
    print(f"  ✓ Key outcome at {key_outcome_idx}")

# ── METHODS ──
print("Highlighting Methods...")
# CoordRep-Rosetta benchmark method
rosetta_method_idx = find_para_containing("CoordRep-Rosetta hard-control benchmark")
if rosetta_method_idx is not None:
    highlight_paragraph(doc.paragraphs[rosetta_method_idx])
    print(f"  ✓ Rosetta method at {rosetta_method_idx}")

# Reproducibility / GitHub
repro_idx = find_para_containing("Reproducibility. All deterministic")
if repro_idx is not None:
    highlight_paragraph(doc.paragraphs[repro_idx])
    print(f"  ✓ Reproducibility method at {repro_idx}")

# Full-CSD PathFinder method
pathfinder_method_idx = find_para_containing("Full-CSD CoordRep-PathFinder analysis")
if pathfinder_method_idx is not None:
    highlight_paragraph(doc.paragraphs[pathfinder_method_idx])
    print(f"  ✓ PathFinder method at {pathfinder_method_idx}")

# Representation-operation benchmarks (includes graph baselines)
repr_bench_idx = find_para_containing("Representation-operation and learning benchmarks")
if repr_bench_idx is not None:
    highlight_paragraph(doc.paragraphs[repr_bench_idx])
    print(f"  ✓ Representation benchmarks method at {repr_bench_idx}")

# ── ASSOCIATED CONTENT ──
print("Highlighting Associated Content...")
accession_idx = find_para_containing("Accession Code. CoordRep")
if accession_idx is not None:
    highlight_paragraph(doc.paragraphs[accession_idx])
    print(f"  ✓ Accession code at {accession_idx}")

# ═══════════════════════════════════════════════════════════════════
# Also highlight all NEW references (131+)
# ═══════════════════════════════════════════════════════════════════
print("Highlighting new references...")
new_ref_count = 0
for i, p in enumerate(doc.paragraphs):
    txt = p.text.strip()
    # Highlight references that are clearly new (Rasmussen, MetalloGen, etc.)
    new_ref_authors = [
        "Rasmussen", "MetalloGen", "Lee, K.; Park", "MACE",
        "Chernyshov", "Pidko", "Molassembler", "Sobez",
        "Kevlishvili", "T-REX", "Nandy, A.; Duan, C.; Taylor",
        "Vogiatzis", "Durand, D. J.; Fey", "Fey, N.; Tsipis",
        "Fey, N.; Harris", "Cooney, K. D.; Cundari",
        "Vela, S.; Laplaza", "cell2mol", "Corminboeuf",
        "Toney, J. W.", "Blackman, A. G.",
        "Blanke, G.", "Brammer, J. C.", "TUCAN",
        "van Staalduinen", "MolBar", "Schneider, N.; Sayle",
        "Shervashidze", "Weisfeiler",
        "Satorras", "E(n) Equivariant"
    ]
    for author in new_ref_authors:
        if author in txt:
            highlight_paragraph(p)
            new_ref_count += 1
            break
print(f"  ✓ Highlighted {new_ref_count} new references")

# ═══════════════════════════════════════════════════════════════════
# Save
# ═══════════════════════════════════════════════════════════════════
doc.save(OUT_HIGHLIGHTED)
print(f"\n✓ Highlighted version: {OUT_HIGHLIGHTED.name}")
print("Done.")
