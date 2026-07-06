"""Render transcript-structure and protein-domain fusion diagrams as PNG."""
from __future__ import annotations

import hashlib
from typing import Optional

from .core import aa3

CATEGORY_KEYWORDS = [
    ("kinase", "#e8590c"),
    ("wd40", "#2f9e44"),
    ("beta-propeller", "#2f9e44"),
    ("help", "#0c8599"),
    ("mam domain", "#1971c2"),
]
FALLBACK_PALETTE = [
    "#495057", "#c2255c", "#5f3dc4", "#0b7285", "#e67700", "#1864ab", "#862e9c",
]
PROMOTER_WIDTH = 34.0
EXON_GAP = 12.0
EXON_MIN_WIDTH = 18.0
EXON_MAX_WIDTH = 34.0
LABEL_HALF_WIDTH_SCALE = 0.007


def color_for(name: str) -> str:
    lower = name.lower()
    for keyword, color in CATEGORY_KEYWORDS:
        if keyword in lower:
            return color
    digest = int(hashlib.md5(lower.encode()).hexdigest(), 16)
    return FALLBACK_PALETTE[digest % len(FALLBACK_PALETTE)]


def structure_segment_color(kind: str) -> str:
    if kind == "coding":
        return "#4c6ef5"
    if kind == "utr5":
        return "#a5d8ff"
    return "#d0ebff"


def _overlaps(a, b):
    return not (a["end"] < b["start"] or a["start"] > b["end"])


def _cluster_by_overlap(items):
    clusters = []
    for it in items:
        merged = [it]
        remaining = []
        for c in clusters:
            if any(_overlaps(it, m) for m in c):
                merged += c
            else:
                remaining.append(c)
        clusters = remaining + [merged]
    return clusters


def _span_overlap_frac(span, clusters):
    s, e = span
    total = e - s + 1
    covered = 0
    for c in clusters:
        cs, ce = min(m["start"] for m in c), max(m["end"] for m in c)
        os_, oe_ = max(s, cs), min(e, ce)
        if oe_ >= os_:
            covered += oe_ - os_ + 1
    return covered / total if total > 0 else 0.0


def canonicalize_domains(raw_domains, gene):
    KEEP_TYPES = {"domain", "repeat", "conserved_site"}
    items = [d for d in raw_domains
             if d["gene"] == gene and d["type"] in KEEP_TYPES and d["name"] != d["accession"]]

    kept = []
    for type_ in ("repeat", "domain", "conserved_site"):
        for c in _cluster_by_overlap([d for d in items if d["type"] == type_]):
            span = (min(m["start"] for m in c), max(m["end"] for m in c))
            if kept and _span_overlap_frac(span, kept) >= 0.5:
                continue
            kept.append(c)

    representatives = []
    for cluster in kept:
        start = min(m["start"] for m in cluster)
        end = max(m["end"] for m in cluster)
        name = min((m["name"] for m in cluster), key=len)
        statuses = {m["status"] for m in cluster}
        if statuses == {"RETAINED"}:
            status = "RETAINED"
        elif statuses == {"LOST"}:
            status = "LOST"
        else:
            status = "DISRUPTED"
        representatives.append({"name": name, "start": start, "end": end, "status": status})
    representatives.sort(key=lambda d: d["start"])
    return representatives


def _label_rows(domains, protein_len):
    groups = []
    for d in domains:
        gap_tolerance = protein_len * 0.03
        if groups and groups[-1]["name"] == d["name"] and d["start"] - groups[-1]["end"] <= gap_tolerance:
            groups[-1]["end"] = max(groups[-1]["end"], d["end"])
        else:
            groups.append({"name": d["name"], "start": d["start"], "end": d["end"]})

    row_free_at = []
    placements = []
    for g in groups:
        center = (g["start"] + g["end"]) / 2
        half_width = len(g["name"]) * protein_len * LABEL_HALF_WIDTH_SCALE
        row = 0
        while row < len(row_free_at) and row_free_at[row] > center - half_width:
            row += 1
        if row == len(row_free_at):
            row_free_at.append(center + half_width)
        else:
            row_free_at[row] = center + half_width
        placements.append((center, row, g["name"]))
    return placements


def _edge_text_position(center: float, text: str, span: float) -> tuple[float, str]:
    half_width = len(text) * span * LABEL_HALF_WIDTH_SCALE
    if center - half_width < 0:
        return 0.0, "left"
    if center + half_width > span:
        return span, "right"
    return center, "center"


def _draw_track(plt_module, ax, label, protein_len, domains, breakpoint_aa=None,
                breakpoint_label=None, junction_label=None):
    from matplotlib.patches import Rectangle

    placements = _label_rows(domains, protein_len)
    n_rows = max((row for _, row, _ in placements), default=0) + 1
    top = 1.0 + 0.28 * (n_rows - 1)

    ax.set_xlim(0, protein_len)
    ax.set_ylim(0, top + 0.45)
    ax.add_patch(Rectangle((0, 0.3), protein_len, 0.4, facecolor="#e9ecef", edgecolor="#ced4da"))
    for d in domains:
        width = max(d["end"] - d["start"], 1)
        alpha = 1.0 if d["status"] == "RETAINED" else 0.85 if d["status"] == "DISRUPTED" else 0.35
        linestyle = "dashed" if d["status"] == "DISRUPTED" else "solid"
        ax.add_patch(Rectangle((d["start"], 0.3), width, 0.4, facecolor=color_for(d["name"]),
                               edgecolor="black", linewidth=0.6, alpha=alpha, linestyle=linestyle))
    for center, row, name in placements:
        text_x, ha = _edge_text_position(center, name, protein_len)
        ax.text(text_x, 0.78 + row * 0.28, name, ha=ha, va="bottom", fontsize=9, clip_on=False)
    if breakpoint_aa is not None:
        ax.axvline(breakpoint_aa, color="#e03131", linestyle="dashed", linewidth=1.5)
        if breakpoint_label:
            text_x, ha = _edge_text_position(breakpoint_aa, breakpoint_label, protein_len)
            ax.text(text_x, top + 0.08, breakpoint_label, ha=ha, va="bottom",
                    fontsize=10, fontweight="bold", color="#e03131", clip_on=False)
        if junction_label:
            label_x = breakpoint_aa + protein_len * 0.12
            if label_x > protein_len * 0.82:
                label_x = max(protein_len * 0.18, breakpoint_aa - protein_len * 0.12)
            ax.annotate(junction_label, xy=(breakpoint_aa, 0.7),
                        xytext=(label_x, top + 0.08),
                        fontsize=10, fontweight="bold", color="#e03131",
                        arrowprops=dict(arrowstyle="-", color="#e03131", linewidth=1))
    ax.set_yticks([])
    ax.set_title(label, loc="left", fontsize=13, pad=14 + 14 * n_rows)
    for spine in ("top", "left", "right"):
        ax.spines[spine].set_visible(False)


def _exon_width(length: int) -> float:
    return max(EXON_MIN_WIDTH, min(EXON_MAX_WIDTH, 10 + (max(length, 1) ** 0.5) * 0.8))


def _layout_transcript_structure(structure: dict) -> dict:
    exons = []
    cursor = PROMOTER_WIDTH + EXON_GAP
    for exon in structure["exons"]:
        width = _exon_width(exon["length"])
        exons.append({
            "rank": exon["rank"],
            "start": cursor,
            "end": cursor + width,
            "width": width,
            "length": exon["length"],
            "segments": exon["segments"],
        })
        cursor += width + EXON_GAP
    width = (exons[-1]["end"] + EXON_GAP / 2) if exons else (PROMOTER_WIDTH + EXON_GAP)
    return {"promoter_start": 0.0, "promoter_end": PROMOTER_WIDTH, "width": width, "exons": exons}


def _transcript_breakpoint_x(partner: dict, layout: dict) -> Optional[float]:
    ctx = partner["breakpoint"]["context"]
    if ctx["region"] == "upstream":
        return (layout["promoter_start"] + layout["promoter_end"]) / 2
    if ctx["region"] == "downstream":
        return layout["width"] - EXON_GAP / 4
    if ctx["intron_rank"] is not None:
        if not (1 <= ctx["intron_rank"] < len(layout["exons"])):
            return None
        left = layout["exons"][ctx["intron_rank"] - 1]
        right = layout["exons"][ctx["intron_rank"]]
        return (left["end"] + right["start"]) / 2
    if ctx["exon_rank"] is not None:
        if not (1 <= ctx["exon_rank"] <= len(layout["exons"])):
            return None
        exon = layout["exons"][ctx["exon_rank"] - 1]
        if ctx["boundary"] == "before":
            return exon["start"]
        if ctx["boundary"] == "after":
            return exon["end"]
        if ctx["exon_offset"] is not None and ctx["exon_length"] is not None:
            frac = 0.5 if ctx["exon_length"] <= 1 else (ctx["exon_offset"] - 1) / (ctx["exon_length"] - 1)
            return exon["start"] + frac * exon["width"]
        return (exon["start"] + exon["end"]) / 2
    return None


def _transcript_breakpoint_label(partner: dict) -> str:
    bp = partner["breakpoint"]
    if bp["type"] == "genomic" and bp.get("genomic_position") is not None:
        return f"g.{bp['genomic_position']} · {bp['context']['label']}"
    return bp["context"]["label"]


def _draw_transcript_track(ax, label: str, partner: dict):
    from matplotlib.patches import Rectangle

    structure = partner.get("structure")
    if not structure:
        return

    layout = _layout_transcript_structure(structure)
    body_y = 0.38
    label_y = body_y + 0.46
    height = 0.34
    strand = "+ strand" if structure["strand"] == 1 else "- strand"

    ax.set_xlim(0, layout["width"])
    ax.set_ylim(0, 1.35)
    ax.set_yticks([])
    ax.set_title(f"{label}  ({partner['transcript']} · {strand})", loc="left", fontsize=13, pad=18)

    promoter_x, promoter_ha = _edge_text_position(
        (layout["promoter_start"] + layout["promoter_end"]) / 2,
        "promoter",
        layout["width"],
    )
    ax.text(promoter_x, label_y, "promoter", ha=promoter_ha, va="bottom", fontsize=9, clip_on=False)
    ax.add_patch(Rectangle((layout["promoter_start"], body_y),
                           layout["promoter_end"] - layout["promoter_start"], height,
                           facecolor="#fff3bf", edgecolor="#c92a2a",
                           linewidth=0.8, linestyle="dashed"))

    for i, exon in enumerate(layout["exons"]):
        if i == 0:
            ax.plot([layout["promoter_end"], exon["start"]],
                    [body_y + height / 2, body_y + height / 2],
                    color="#868e96", linewidth=1.2)
        else:
            prev = layout["exons"][i - 1]
            ax.plot([prev["end"], exon["start"]],
                    [body_y + height / 2, body_y + height / 2],
                    color="#868e96", linewidth=1.2)
        ax.text((exon["start"] + exon["end"]) / 2, label_y, str(exon["rank"]),
                ha="center", va="bottom", fontsize=9, clip_on=False)
        for seg in exon["segments"]:
            seg_start = exon["start"] + ((seg["start"] - 1) / exon["length"]) * exon["width"]
            seg_end = exon["start"] + (seg["end"] / exon["length"]) * exon["width"]
            ax.add_patch(Rectangle((seg_start, body_y), max(seg_end - seg_start, 0.8), height,
                                   facecolor=structure_segment_color(seg["kind"]), edgecolor="none"))
        ax.add_patch(Rectangle((exon["start"], body_y), exon["width"], height,
                               facecolor="none", edgecolor="#495057", linewidth=0.8))

    bp_x = _transcript_breakpoint_x(partner, layout)
    if bp_x is not None:
        ax.axvline(bp_x, color="#e03131", linestyle="dashed", linewidth=1.5)
        bp_label = _transcript_breakpoint_label(partner)
        text_x, ha = _edge_text_position(bp_x, bp_label, layout["width"])
        ax.text(text_x, 1.07, bp_label,
                ha=ha, va="bottom", fontsize=10, fontweight="bold", color="#e03131", clip_on=False)

    for spine in ("top", "left", "right"):
        ax.spines[spine].set_visible(False)


def render_domain_diagram_png(
    result: dict,
    *,
    five_uniprot: Optional[str] = None,
    three_uniprot: Optional[str] = None,
    title: Optional[str] = None,
) -> bytes:
    import io

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    iface = result["interface"]
    resolved = result["resolved"]
    five_gene, three_gene = iface["five_gene"], iface["three_gene"]
    five_len = resolved["five"]["protein_length"]
    three_len = resolved["three"]["protein_length"]

    five_domains = canonicalize_domains(iface["domains"], five_gene)
    three_domains = canonicalize_domains(iface["domains"], three_gene)

    five_last_aa = iface["five_last_aa"]
    three_first_aa = iface["three_first_aa"]
    hybrid_offset = 1 if iface["hybrid_codon"] else 0
    three_offset = five_last_aa + hybrid_offset - three_first_aa + 1

    fusion_domains = []
    for d in five_domains:
        if d["status"] == "LOST":
            continue
        fusion_domains.append({"name": d["name"], "status": d["status"],
                               "start": d["start"], "end": min(d["end"], five_last_aa)})
    for d in three_domains:
        if d["status"] == "LOST":
            continue
        start = max(d["start"], three_first_aa) + three_offset
        end = d["end"] + three_offset
        fusion_domains.append({"name": d["name"], "status": "RETAINED", "start": start, "end": end})
    fusion_domains.sort(key=lambda d: d["start"])

    tracks = []
    if resolved["five"].get("structure"):
        tracks.append(("structure", resolved["five"], f"{five_gene} transcript structure"))
    tracks.append(("protein", "five"))
    if resolved["three"].get("structure"):
        tracks.append(("structure", resolved["three"], f"{three_gene} transcript structure"))
    tracks.append(("protein", "three"))
    tracks.append(("fusion", None))

    fig_height = 2.15 * len(tracks) + 1.4
    fig, axes = plt.subplots(len(tracks), 1, figsize=(11, fig_height))
    if not isinstance(axes, (list, tuple)):
        try:
            axes = list(axes)
        except TypeError:
            axes = [axes]

    fig.suptitle(
        title or f"{five_gene}::{three_gene} — transcript structure, chimeric protein & domain retention",
        fontsize=16, fontweight="bold", x=0.02, ha="left", y=0.985)

    five_title = f"{five_gene}  (5' partner"
    if five_uniprot:
        five_title += f" · {five_uniprot}"
    five_title += f" · {five_len} aa)"
    three_title = f"{three_gene}  (3' partner"
    if three_uniprot:
        three_title += f" · {three_uniprot}"
    three_title += f" · {three_len} aa)"
    if five_last_aa > 0:
        junction = f"junction p.{aa3(iface['five_last_aa_res'])}{five_last_aa}::" \
                   f"{aa3(iface['three_first_aa_res'])}{three_first_aa}"
    else:
        junction = f"junction p.0::{aa3(iface['three_first_aa_res'])}{three_first_aa}"
    if iface["hybrid_codon"] and iface.get("junction_residue"):
        junction += f"  (hybrid {aa3(iface['junction_residue'])}{five_last_aa + 1})"
    fusion_title = (f"{five_gene}::{three_gene} fusion protein  ({iface['fusion_length']} aa · "
                    f"{iface['frame_status']})")

    ax_iter = iter(axes)
    for track in tracks:
        kind = track[0]
        ax = next(ax_iter)
        if kind == "structure":
            _draw_transcript_track(ax, track[2], track[1])
        elif kind == "protein" and track[1] == "five":
            _draw_track(plt, ax, five_title, five_len, five_domains, breakpoint_aa=five_last_aa,
                        breakpoint_label=f"breakpoint aa {five_last_aa}")
        elif kind == "protein":
            _draw_track(plt, ax, three_title, three_len, three_domains, breakpoint_aa=three_first_aa,
                        breakpoint_label=f"breakpoint aa {three_first_aa}")
        else:
            _draw_track(plt, ax, fusion_title, iface["fusion_length"], fusion_domains,
                        breakpoint_aa=five_last_aa, junction_label=junction)
            ax.set_xlabel("residue")

    plt.tight_layout(rect=(0, 0, 1, 0.97))
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150)
    plt.close(fig)
    return buf.getvalue()
