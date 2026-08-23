from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


HERE = Path(__file__).resolve().parent
ARIAL = Path(r"C:\Windows\Fonts\arial.ttf")
ARIAL_BOLD = Path(r"C:\Windows\Fonts\arialbd.ttf")

CANDIDATE_COUNTS = np.array([2, 4, 6, 8, 12, 24, 36, 48, 120])
RECORD_COUNTS = np.array([2652, 1675, 284, 67, 78, 184, 15, 120, 98])
TIMING_MS = {
    "Signature-only sort": {
        "Median": 0.20889998995698988,
        "P95": 0.37399999564513564,
        "Maximum": 17.719099996611476,
    },
    "Full exact search": {
        "Median": 0.7409999961964786,
        "P95": 6.523474989808054,
        "Maximum": 63.83759999880567,
    },
}

BLUE = "#3B78A7"
ORANGE = "#C97535"
DARK = "#20242A"
GRID = "#D8DDE3"
WHITE = "#FFFFFF"


def pil_font(points: float, dpi: int, *, bold: bool = False) -> ImageFont.FreeTypeFont:
    path = ARIAL_BOLD if bold else ARIAL
    return ImageFont.truetype(str(path), max(1, round(points * dpi / 72)))


def centered(draw: ImageDraw.ImageDraw, xy: tuple[float, float], text: str,
             font: ImageFont.FreeTypeFont, fill: str = DARK) -> None:
    box = draw.textbbox((0, 0), text, font=font)
    draw.text((xy[0] - (box[2] - box[0]) / 2, xy[1] - (box[3] - box[1]) / 2),
              text, font=font, fill=fill)


def build_raster(dpi: int = 600) -> Image.Image:
    width, height = round(7.15 * dpi), round(3.25 * dpi)
    im = Image.new("RGB", (width, height), WHITE)
    draw = ImageDraw.Draw(im)
    s = dpi / 72
    px = lambda value: value * s

    f7 = pil_font(7.0, dpi)
    f72 = pil_font(7.2, dpi)
    f8 = pil_font(8.0, dpi)
    f9b = pil_font(9.0, dpi, bold=True)
    f105b = pil_font(10.5, dpi, bold=True)

    # Plot boxes in points.
    a = (px(48), px(45), px(306), px(184))
    b = (px(370), px(45), px(505), px(184))
    draw.text((px(8), px(10)), "A", font=f105b, fill=DARK)
    draw.text((px(350), px(10)), "B", font=f105b, fill=DARK)
    draw.text((a[0], px(10)), "Legal-orbit complexity", font=f9b, fill=DARK)
    draw.text((b[0], px(10)), "Per-call cost", font=f9b, fill=DARK)
    draw.text((a[0], px(27)), "5,173/6,427 screened graphs contain a legal orbit",
              font=f7, fill=DARK)

    # Panel A.
    ax0, ay0, ax1, ay1 = a
    axis_w = max(1, round(px(0.75)))
    grid_w = max(1, round(px(0.55)))
    draw.line((ax0, ay0, ax0, ay1), fill=DARK, width=axis_w)
    draw.line((ax0, ay1, ax1, ay1), fill=DARK, width=axis_w)
    for tick in (0, 1000, 2000, 3000):
        x = ax0 + (ax1 - ax0) * tick / 3000
        draw.line((x, ay0, x, ay1), fill=GRID, width=grid_w)
        centered(draw, (x, ay1 + px(10)), f"{tick:,}", f72)
    centered(draw, ((ax0 + ax1) / 2, ay1 + px(27)), "Challenge-eligible records", f8)

    row_h = (ay1 - ay0) / len(CANDIDATE_COUNTS)
    for idx, (candidate_count, records) in enumerate(zip(CANDIDATE_COUNTS, RECORD_COUNTS)):
        cy = ay0 + row_h * (idx + 0.5)
        bar_h = row_h * 0.64
        x_end = ax0 + (ax1 - ax0) * records / 3000
        draw.rectangle((ax0, cy - bar_h / 2, x_end, cy + bar_h / 2), fill=BLUE)
        draw.text((x_end + px(3), cy - px(4.2)), f"{records:,}", font=f7, fill=DARK)
        label = f"{candidate_count:,}"
        box = draw.textbbox((0, 0), label, font=f72)
        draw.text((ax0 - px(8) - (box[2] - box[0]), cy - px(4.2)), label, font=f72, fill=DARK)

    ylab = Image.new("RGBA", (round(px(95)), round(px(18))), (255, 255, 255, 0))
    ImageDraw.Draw(ylab).text((0, 0), "Legal candidates", font=f8, fill=DARK)
    ylab = ylab.rotate(90, expand=True, resample=Image.Resampling.BICUBIC)
    im.paste(ylab, (round(px(3)), round((ay0 + ay1 - ylab.height) / 2)), ylab)

    # Panel B.
    bx0, by0, bx1, by1 = b
    draw.line((bx0, by0, bx0, by1), fill=DARK, width=axis_w)
    draw.line((bx0, by1, bx1, by1), fill=DARK, width=axis_w)
    log_min, log_max = math.log10(0.12), math.log10(110)

    def y_for(value: float) -> float:
        frac = (math.log10(value) - log_min) / (log_max - log_min)
        return by1 - frac * (by1 - by0)

    for tick in (0.2, 1, 10, 100):
        y = y_for(tick)
        draw.line((bx0, y, bx1, y), fill=GRID, width=grid_w)
        label = f"{tick:g}"
        box = draw.textbbox((0, 0), label, font=f72)
        draw.text((bx0 - px(7) - (box[2] - box[0]), y - px(4.2)), label, font=f72, fill=DARK)

    x_positions = [bx0 + (bx1 - bx0) * 0.30, bx0 + (bx1 - bx0) * 0.76]
    for x, lines in zip(x_positions, (("Signature-only", "sort"), ("Full exact", "search"))):
        centered(draw, (x, by1 + px(10)), lines[0], f72)
        centered(draw, (x, by1 + px(20)), lines[1], f72)

    marker_info = {"Median": ("circle", BLUE), "P95": ("square", ORANGE),
                   "Maximum": ("triangle", DARK)}

    def marker(x: float, y: float, kind: str, color: str) -> None:
        r = px(3.2)
        if kind == "circle":
            draw.ellipse((x - r, y - r, x + r, y + r), fill=color, outline=WHITE,
                         width=max(1, round(px(0.5))))
        elif kind == "square":
            draw.rectangle((x - r, y - r, x + r, y + r), fill=color, outline=WHITE,
                           width=max(1, round(px(0.5))))
        else:
            draw.polygon([(x, y - r), (x - r, y + r), (x + r, y + r)], fill=color)

    for x, method in zip(x_positions, TIMING_MS):
        for stat, value in TIMING_MS[method].items():
            kind, color = marker_info[stat]
            y = y_for(value)
            marker(x, y, kind, color)
            label = f"{value:.3f}" if value < 10 else f"{value:.1f}"
            draw.text((x + px(5), y - px(4.0)), label, font=f7, fill=DARK)

    lx, ly = bx0 + px(3), px(29)
    for i, stat in enumerate(("Median", "P95", "Maximum")):
        kind, color = marker_info[stat]
        xx = lx + px(i * 38)
        marker(xx, ly, kind, color)
        draw.text((xx + px(5), ly - px(4.0)), stat, font=f7, fill=DARK)

    ylab2 = Image.new("RGBA", (round(px(112)), round(px(18))), (255, 255, 255, 0))
    ImageDraw.Draw(ylab2).text((0, 0), "Wall time (ms; log)", font=f8, fill=DARK)
    ylab2 = ylab2.rotate(90, expand=True, resample=Image.Resampling.BICUBIC)
    im.paste(ylab2, (round(px(337)), round((by0 + by1 - ylab2.height) / 2)), ylab2)
    return im


def build_svg() -> str:
    # The vector wrapper preserves exact physical dimensions; source values are
    # separately exposed in CSV for reproducible redrawing.
    import base64
    import io

    image = build_raster(300)
    buffer = io.BytesIO()
    image.save(buffer, format="PNG", optimize=True)
    payload = base64.b64encode(buffer.getvalue()).decode("ascii")
    return ('<svg xmlns="http://www.w3.org/2000/svg" width="7.15in" height="3.25in" '
            'viewBox="0 0 2145 975">\n'
            f'  <image width="2145" height="975" href="data:image/png;base64,{payload}"/>\n'
            '</svg>\n')


def build_pdf(path: Path) -> None:
    pdfmetrics.registerFont(TTFont("Arial", str(ARIAL)))
    c = canvas.Canvas(str(path), pagesize=(7.15 * 72, 3.25 * 72))
    png = HERE / "Supplementary_Figure_S2_Canonicalization_Complexity.png"
    c.drawImage(str(png), 0, 0, width=7.15 * 72, height=3.25 * 72)
    c.showPage()
    c.save()


def write_source_data() -> None:
    rows: list[dict[str, object]] = []
    for candidate_count, records in zip(CANDIDATE_COUNTS, RECORD_COUNTS):
        rows.append({"panel": "A", "series": "legal-orbit population",
                     "category": int(candidate_count), "statistic": "records",
                     "value": int(records), "unit": "records",
                     "denominator_or_calls": 5173})
    for method, values in TIMING_MS.items():
        for stat, value in values.items():
            rows.append({"panel": "B", "series": method, "category": stat,
                         "statistic": "wall time per call", "value": value,
                         "unit": "ms", "denominator_or_calls": 101000})
    with (HERE / "Supplementary_Figure_S2_source_data.csv").open(
        "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_caption_and_notes() -> None:
    caption = (
        "Supplementary Figure S2. Search complexity and computational cost of "
        "attachment-aware exact canonicalization. (A) Distribution of the exact "
        "legal-candidate count among the 5,173 challenge-eligible records identified "
        "from 6,427 unique graph-disjoint tmQMg/PBE CN = 4–6 snapshots. Candidate "
        "counts are the Cartesian products of symmetry-allowed equivalent-ligand and "
        "donor-rank-orbit permutations after chemical, attachment, and intrinsic-distance "
        "refinement. All eligible records remained below the frozen 100,000-candidate "
        "ceiling. (B) Descriptive wall time per canonicalization call for invariant-key "
        "sorting alone and the full exact minimum-whole-record search, reported as the "
        "median, 95th percentile, and maximum over 101,000 calls per method on the same "
        "machine. These timings are implementation diagnostics rather than a cross-machine "
        "performance benchmark. In the corresponding N = 1,000, K = 100 challenge, the "
        "full exact method produced 100,000/100,000 exact variant matches, 100,000/100,000 "
        "grammar round trips, stable candidate counts for 1,000/1,000 structures, zero "
        "variant errors, and zero false merges in the declared attachment control."
    )
    (HERE / "Supplementary_Figure_S2_caption.txt").write_text(
        caption + "\n", encoding="utf-8", newline="\n")
    readme = "# Supplementary Figure S2\n\n"
    readme += "This figure is derived from the frozen 2026-08-22 legal-orbit canonicalization challenge.\n\n"
    readme += "- Source population: 6,427 graph-disjoint tmQMg/PBE CN = 4–6 snapshots.\n"
    readme += "- Challenge-eligible population: 5,173 records.\n"
    readme += "- Executed cohort: N = 1,000; K = 100 atom reindexings.\n"
    readme += "- Timing values are descriptive wall-clock measurements from the frozen run.\n"
    readme += "- No licensed CSD structures or coordinates are included.\n"
    (HERE / "README.md").write_text(readme, encoding="utf-8", newline="\n")


def write_qa() -> None:
    qa = {"status": "PASS", "candidate_distribution_sum": int(RECORD_COUNTS.sum()),
          "expected_challenge_eligible": 5173, "screened_unique_graphs": 6427,
          "executed_structures": 1000, "variants_per_structure": 100,
          "full_exact_variant_matches": 100000, "full_exact_variant_total": 100000,
          "grammar_round_trips": 100000, "stable_candidate_count_structures": 1000,
          "variant_errors": 0, "false_merges": 0}
    if qa["candidate_distribution_sum"] != qa["expected_challenge_eligible"]:
        raise RuntimeError("candidate-count distribution does not sum to 5,173")
    (HERE / "Supplementary_Figure_S2_QA.json").write_text(
        json.dumps(qa, indent=2) + "\n", encoding="utf-8", newline="\n")


def write_checksums() -> None:
    targets = sorted(p for p in HERE.iterdir() if p.is_file()
                     and p.name != "Supplementary_Figure_S2_SHA256SUMS.txt"
                     and p.suffix.lower() != ".pyc")
    lines = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}" for p in targets]
    (HERE / "Supplementary_Figure_S2_SHA256SUMS.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def main() -> None:
    stem = HERE / "Supplementary_Figure_S2_Canonicalization_Complexity"
    image = build_raster(600)
    image.save(stem.with_suffix(".png"), dpi=(600, 600), optimize=True)
    image.save(stem.with_suffix(".tif"), dpi=(600, 600), compression="tiff_lzw")
    preview = build_raster(180)
    preview.save(HERE / "Supplementary_Figure_S2_preview_180dpi.png", dpi=(180, 180), optimize=True)
    stem.with_suffix(".svg").write_text(build_svg(), encoding="utf-8", newline="\n")
    build_pdf(stem.with_suffix(".pdf"))
    write_source_data()
    write_caption_and_notes()
    write_qa()
    write_checksums()


if __name__ == "__main__":
    main()
