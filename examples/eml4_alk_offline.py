#!/usr/bin/env python3
"""Annotate EML4::ALK variant 1 (E13;A20) fully offline, from the bundled fixture.

Run:  python examples/eml4_alk_offline.py
"""
import json
import os

from fusion_annotation import Transcript, build_exon_cds_map, annotate_fusion
from fusion_annotation.providers import StaticProvider

FIXTURE = os.path.join(os.path.dirname(__file__), os.pardir,
                       "tests", "fixtures", "eml4_alk_fixture.json")


def load_provider():
    fx = json.load(open(FIXTURE))
    txs = {}
    for key, t in fx["transcripts"].items():
        exon_cds = build_exon_cds_map(t["strand"], t["exons"], t["cds_g_start"], t["cds_g_end"])
        txs[key] = Transcript(
            gene_symbol=t["gene_symbol"], gene_id=t["gene_id"],
            transcript_id=t["transcript_id"], strand=t["strand"],
            cds=t["cds"], protein=t["protein"], uniprot=t["uniprot"], exon_cds=exon_cds)
    return StaticProvider(txs, domains=fx["domains"], knowledge=fx["knowledge"])


def main():
    provider = load_provider()
    result = annotate_fusion(provider, "EML4", "ALK", five_exon=13, three_exon=20)
    iface, kn = result["interface"], result["knowledge"]

    print("=== EML4::ALK variant 1 (E13;A20) ===\n")
    print("HGVS.p-like :", iface["hgvsp_like"])
    print("categorical :", iface["categorical_key"])
    print(f"frame       : {iface['frame_status']}  "
          f"(protein {iface['fusion_length']} aa, internal stops {iface['internal_stops']})")
    print(f"junction    : {iface['five_gene']} Lys{iface['five_last_aa']} :: "
          f"{iface['three_gene']} Tyr{iface['three_first_aa']}  "
          f"(hybrid codon -> {iface['junction_residue'] or '-'})")
    print("\nretained domains:")
    for d in iface["domains"]:
        if d["status"] == "RETAINED":
            print(f"  [{d['status']}] {d['name']} ({d['start']}-{d['end']})")
    print("\nknowledge:")
    print("  oncogenic :", kn["oncogenic"])
    print("  therapies :", ", ".join(kn["therapies"]))
    print("  sources   :", ", ".join(kn["sources"]))


if __name__ == "__main__":
    main()
