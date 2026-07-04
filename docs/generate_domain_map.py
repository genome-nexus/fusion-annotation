#!/usr/bin/env python3
"""Regenerate docs/fusion_domain_map.png from the bundled EML4::ALK fixture.

Fixes GH issue #17: the previous (hand-made, un-scripted) PNG assigned colors
to domains independently per track, so the *same* biological domain (the ALK
kinase domain) ended up orange in the ALK-only track but blue in the fusion
track. This script derives colors from a single category->color table keyed
on domain name, so a given domain always renders in the same color no matter
which track it appears in.

It also de-duplicates the raw InterPro/Pfam domain hits: Genome Nexus reports
many overlapping annotations for the same physical domain (e.g. 7 different
records all covering the ALK kinase region: Pfam PF07714, InterPro "Protein
kinase domain", "...catalytic domain" x2, "...active site", "...ATP binding
site", "...conserved site", the PANTHER-style "Receptor Tyrosine Kinase"
family, and the "Protein kinase-like domain superfamily"). Rendering all of
these as separate overlapping rectangles is illegible, so overlapping
same-type records are merged into one representative block per domain.

Run:  pip install matplotlib && python docs/generate_domain_map.py
"""
import hashlib
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from fusion_annotation import Transcript, build_exon_cds_map, build_exon_genomic_map, annotate_fusion
from fusion_annotation.providers import StaticProvider

HERE = os.path.dirname(__file__)
FIXTURE = os.path.join(HERE, os.pardir, "tests", "fixtures", "eml4_alk_fixture.json")
OUT_PNG = os.path.join(HERE, "fusion_domain_map.png")

# ---------------------------------------------------------------------------
# Domain -> color: a single lookup used for every track, so the same domain
# always renders identically wherever it appears. Curated colors for the
# domain categories that show up in this example; anything else falls back
# to a deterministic (md5-hash-based) pick from FALLBACK_PALETTE so unknown
# domains still get a stable, distinguishable color.
# ---------------------------------------------------------------------------
CATEGORY_KEYWORDS = [
    ("kinase", "#e8590c"),          # orange
    ("wd40", "#2f9e44"),            # green
    ("beta-propeller", "#2f9e44"),  # green (same family as WD40 repeats)
    ("help", "#0c8599"),            # teal
    ("mam domain", "#1971c2"),      # blue
]
FALLBACK_PALETTE = [
    "#495057", "#c2255c", "#5f3dc4", "#0b7285", "#e67700", "#1864ab", "#862e9c",
]


def color_for(name: str) -> str:
    lower = name.lower()
    for keyword, color in CATEGORY_KEYWORDS:
        if keyword in lower:
            return color
    digest = int(hashlib.md5(lower.encode()).hexdigest(), 16)
    return FALLBACK_PALETTE[digest % len(FALLBACK_PALETTE)]


def _overlaps(a, b):
    return not (a["end"] < b["start"] or a["start"] > b["end"])


def _cluster_by_overlap(items):
    """Group items that pairwise overlap into connected-component clusters."""
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
    return covered / total


def canonicalize_domains(raw_domains, gene):
    """Collapse the raw, highly redundant per-source domain hits for one gene
    into a small set of representative blocks.

    Genome Nexus returns many overlapping InterPro/Pfam records describing
    the same physical domain (different databases, different exact
    boundaries). We keep ``domain``/``repeat``/``conserved_site`` type
    records with a real curated name (dropping bare Pfam accessions used as
    a placeholder name), cluster each type separately by overlap, and — most
    specific/granular first (repeat, then domain, then conserved_site) —
    drop any later cluster that is >=50% covered by an already-kept cluster
    (a coarser "superfamily-ish" or sub-feature re-annotation of the same
    region, e.g. the tiny "ATP binding site" motif fully inside the already-
    kept kinase domain block).
    """
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


def load_result():
    fx = json.load(open(FIXTURE))
    txs = {}
    uniprot_by_gene = {}
    for key, t in fx["transcripts"].items():
        exon_cds = build_exon_cds_map(t["strand"], t["exons"], t["cds_g_start"], t["cds_g_end"])
        exon_genomic = build_exon_genomic_map(t["strand"], t["exons"])
        txs[key] = Transcript(
            gene_symbol=t["gene_symbol"], gene_id=t["gene_id"],
            transcript_id=t["transcript_id"], strand=t["strand"],
            cds=t["cds"], protein=t["protein"], uniprot=t["uniprot"],
            exon_cds=exon_cds, exon_genomic=exon_genomic,
            cds_g_start=t["cds_g_start"], cds_g_end=t["cds_g_end"], is_canonical=True)
        uniprot_by_gene[t["gene_symbol"]] = t["uniprot"]
    provider = StaticProvider(txs, domains=fx["domains"], knowledge=fx["knowledge"])
    result = annotate_fusion(provider, "EML4", "ALK", five_exon=13, three_exon=20)
    return result, uniprot_by_gene


def _label_rows(domains, protein_len):
    """Group adjacent same-name domains under one shared label, then stagger
    labels into extra rows when their estimated text width would collide
    with a neighboring label — avoids the illegible overlapping-text problem
    you get from labeling every individual rectangle independently."""
    groups = []
    for d in domains:
        gap_tolerance = protein_len * 0.03
        if groups and groups[-1]["name"] == d["name"] and d["start"] - groups[-1]["end"] <= gap_tolerance:
            groups[-1]["end"] = max(groups[-1]["end"], d["end"])
        else:
            groups.append({"name": d["name"], "start": d["start"], "end": d["end"]})

    row_free_at = []  # row index -> x position where that row is next free
    placements = []
    for g in groups:
        center = (g["start"] + g["end"]) / 2
        half_width = len(g["name"]) * protein_len * 0.0055
        row = 0
        while row < len(row_free_at) and row_free_at[row] > center - half_width:
            row += 1
        if row == len(row_free_at):
            row_free_at.append(center + half_width)
        else:
            row_free_at[row] = center + half_width
        placements.append((center, row, g["name"]))
    return placements


def draw_track(ax, label, protein_len, domains, breakpoint_aa=None, breakpoint_label=None,
               junction_label=None):
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
        ax.text(center, 0.78 + row * 0.28, name, ha="center", va="bottom", fontsize=9)
    if breakpoint_aa is not None:
        ax.axvline(breakpoint_aa, color="#e03131", linestyle="dashed", linewidth=1.5)
        ax.text(breakpoint_aa, top + 0.08, breakpoint_label, ha="center", va="bottom",
                fontsize=10, fontweight="bold", color="#e03131")
        if junction_label:
            ax.annotate(junction_label, xy=(breakpoint_aa, 0.7), xytext=(breakpoint_aa + protein_len * 0.12, top + 0.08),
                        fontsize=10, fontweight="bold", color="#e03131",
                        arrowprops=dict(arrowstyle="-", color="#e03131", linewidth=1))
    ax.set_yticks([])
    ax.set_title(label, loc="left", fontsize=13, pad=10 + 14 * n_rows)
    for spine in ("top", "left", "right"):
        ax.spines[spine].set_visible(False)


def main():
    result, uniprot_by_gene = load_result()
    iface = result["interface"]
    five_len = next(d["partner_protein_length"] for d in iface["domains"] if d["gene"] == "EML4")
    three_len = next(d["partner_protein_length"] for d in iface["domains"] if d["gene"] == "ALK")
    five_uniprot = uniprot_by_gene["EML4"]
    three_uniprot = uniprot_by_gene["ALK"]

    five_domains = canonicalize_domains(iface["domains"], "EML4")
    three_domains = canonicalize_domains(iface["domains"], "ALK")

    # Fusion-protein track: remap each partner's domains onto fusion
    # coordinates. 5' domains keep their original numbering (the fusion
    # protein starts with the 5' partner's sequence unchanged up to the
    # breakpoint). 3' domains shift by a constant offset derived from the
    # junction: fusion_pos = five_last_aa + hybrid_offset + (orig_pos -
    # three_first_aa + 1). Only domains that actually survive into the
    # fusion protein (RETAINED, or DISRUPTED clipped at the breakpoint) are
    # shown, since LOST domains are by definition entirely absent from it.
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

    fig, axes = plt.subplots(3, 1, figsize=(11, 9.5))
    fig.suptitle("EML4::ALK variant 1 (E13;A20) — chimeric protein & domain retention",
                 fontsize=16, fontweight="bold", x=0.02, ha="left", y=0.985)

    draw_track(axes[0], f"EML4  (5' partner · {five_uniprot} · {five_len} aa)", five_len,
               five_domains, breakpoint_aa=five_last_aa,
               breakpoint_label=f"breakpoint aa {five_last_aa}\n(exon 13)")
    draw_track(axes[1], f"ALK  (3' partner · {three_uniprot} · {three_len} aa)", three_len,
               three_domains, breakpoint_aa=three_first_aa,
               breakpoint_label=f"breakpoint aa {three_first_aa}\n(exon 20)")
    junction = (f"junction p.Lys{five_last_aa}::Ser{three_first_aa}"
                f" (hybrid Val{five_last_aa + 1})" if iface["hybrid_codon"] else
                f"junction p.Lys{five_last_aa}::Ser{three_first_aa}")
    draw_track(axes[2],
               f"EML4::ALK fusion protein  ({iface['fusion_length']} aa · {iface['frame_status']} · "
               f"no internal stop)" if iface["internal_stops"] == 0 else
               f"EML4::ALK fusion protein  ({iface['fusion_length']} aa · {iface['frame_status']})",
               iface["fusion_length"], fusion_domains, breakpoint_aa=five_last_aa,
               junction_label=junction)
    axes[2].set_xlabel("residue")

    plt.tight_layout(rect=(0, 0, 1, 0.97))
    plt.savefig(OUT_PNG, dpi=150)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
