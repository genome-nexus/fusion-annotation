#!/usr/bin/env python3
"""Regenerate docs/fusion_domain_map.png from the bundled EML4::ALK fixture.

Fixes GH issue #17: the previous (hand-made, un-scripted) PNG assigned colors
to domains independently per track, so the *same* biological domain (the ALK
kinase domain) ended up orange in the ALK-only track but blue in the fusion
track. Rendering now goes through fusion_annotation.domain_diagram, which
derives colors from a single category->color table keyed on domain name (so
a given domain always renders the same color no matter which track it's in)
and de-duplicates the many overlapping InterPro/Pfam records Genome Nexus
returns for the same physical domain. The MCP server (server/app.py) uses
this exact same renderer to attach a diagram to its tool response.

Run:  pip install -e .[docs] && python docs/generate_domain_map.py
"""
import json
import os

from fusion_annotation import Transcript, build_exon_cds_map, build_exon_genomic_map, annotate_fusion
from fusion_annotation.domain_diagram import render_domain_diagram_png
from fusion_annotation.providers import StaticProvider

HERE = os.path.dirname(__file__)
FIXTURE = os.path.join(HERE, os.pardir, "tests", "fixtures", "eml4_alk_fixture.json")
OUT_PNG = os.path.join(HERE, "fusion_domain_map.png")


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


def main():
    result, uniprot_by_gene = load_result()
    png = render_domain_diagram_png(
        result,
        five_uniprot=uniprot_by_gene["EML4"],
        three_uniprot=uniprot_by_gene["ALK"],
        title="EML4::ALK variant 1 (E13;A20) — chimeric protein & domain retention",
    )
    with open(OUT_PNG, "wb") as f:
        f.write(png)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
