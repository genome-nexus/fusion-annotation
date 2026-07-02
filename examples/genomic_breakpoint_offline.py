#!/usr/bin/env python3
"""Annotate EML4::ALK from GENOMIC breakpoints (offline, from the bundled fixture).

A genomic coordinate pins the transcript isoform: the tool maps it through the
exon table to an exact CDS base, so callers don't have to know exon numbering.
This example shows that the genomic breakpoints for EML4::ALK variant 1 land on
exactly the same junction as the exon-number call (E13;A20), and prints the
`resolved` block that echoes how each breakpoint was interpreted.

Coordinates are GRCh38 (both partners are on chr2):
  - EML4 3' breakpoint: chr2:42295516  (3' end of EML4 exon 13, + strand)
  - ALK  5' breakpoint: chr2:29223528  (5' start of ALK exon 20, - strand)

Run:  python examples/genomic_breakpoint_offline.py
"""
import json
import os

from fusion_annotation import (
    Transcript, build_exon_cds_map, build_exon_genomic_map, annotate_fusion,
)
from fusion_annotation.providers import StaticProvider

FIXTURE = os.path.join(os.path.dirname(__file__), os.pardir,
                       "tests", "fixtures", "eml4_alk_fixture.json")

# GRCh38 genomic breakpoints for EML4::ALK variant 1 (verified against the fixture).
EML4_BREAKPOINT = "chr2:42295516"      # also accepts 42295516 (int) or "g.42295516"
ALK_BREAKPOINT = "chr2:29223528"


def load_provider():
    fx = json.load(open(FIXTURE))
    txs = {}
    for key, t in fx["transcripts"].items():
        txs[key] = Transcript(
            gene_symbol=t["gene_symbol"], gene_id=t["gene_id"],
            transcript_id=t["transcript_id"], strand=t["strand"],
            cds=t["cds"], protein=t["protein"], uniprot=t["uniprot"],
            exon_cds=build_exon_cds_map(t["strand"], t["exons"], t["cds_g_start"], t["cds_g_end"]),
            exon_genomic=build_exon_genomic_map(t["strand"], t["exons"]),
            cds_g_start=t["cds_g_start"], cds_g_end=t["cds_g_end"], is_canonical=True)
    return StaticProvider(txs, domains=fx["domains"], knowledge=fx["knowledge"])


def main():
    provider = load_provider()

    result = annotate_fusion(
        provider, "EML4", "ALK",
        five_genomic=EML4_BREAKPOINT, three_genomic=ALK_BREAKPOINT)
    iface = result["interface"]

    print("=== EML4::ALK from genomic breakpoints ===\n")
    print(f"5' breakpoint : EML4 {EML4_BREAKPOINT}")
    print(f"3' breakpoint : ALK  {ALK_BREAKPOINT}\n")
    print("HGVS.p-like :", iface["hgvsp_like"])
    print(f"frame       : {iface['frame_status']}  "
          f"(protein {iface['fusion_length']} aa, internal stops {iface['internal_stops']})")

    print("\nresolved (echoed back by the tool):")
    for side in ("five", "three"):
        r = result["resolved"][side]
        bp = r["breakpoint"]
        print(f"  {r['gene']:5s} {r['transcript']} ({r['transcript_source']}) "
              f"-- {bp['type']} breakpoint g.{bp['genomic_position']} -> CDS coord {bp['cds_coord']}")
    if result["warnings"]:
        print("\nwarnings:")
        for w in result["warnings"]:
            print("  -", w)

    # The genomic call lands on exactly the same junction as the exon-number call.
    by_exon = annotate_fusion(provider, "EML4", "ALK", five_exon=13, three_exon=20)
    same = iface["hgvsp_like"] == by_exon["interface"]["hgvsp_like"]
    print(f"\nmatches the exon-number call (E13;A20)? {same}")


if __name__ == "__main__":
    main()
