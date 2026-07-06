"""Offline regression tests for representative ERG/EGFR cBioPortal cases.

These cases were selected from the `msk_impact_50k_2026` structural-variant
audit and all run through the local Python package using genomic breakpoints.
The fixture pins the exact GRCh37 transcripts used during the audit so the
tests stay offline and deterministic.
"""
from __future__ import annotations

import json
import os

import pytest

from fusion_annotation import Transcript, annotate_fusion, build_exon_cds_map, translate
from fusion_annotation.providers import StaticProvider

FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "cbioportal_grch37_cases.json")


@pytest.fixture(scope="module")
def fixture_data():
    with open(FIXTURE, encoding="utf-8") as handle:
        return json.load(handle)


@pytest.fixture(scope="module")
def provider(fixture_data):
    txs: dict[str, Transcript] = {}
    for tx_id, t in fixture_data["transcripts"].items():
        exons = [{"start": start, "end": end} for start, end in t["exon_genomic"]]
        tx = Transcript(
            gene_symbol=t["gene_symbol"],
            gene_id="",
            transcript_id=tx_id,
            strand=t["strand"],
            cds=t["cds"],
            protein=translate(t["cds"]).split("*")[0],
            exon_cds=build_exon_cds_map(t["strand"], exons, t["cds_g_start"], t["cds_g_end"]),
            exon_genomic=[tuple(e) for e in t["exon_genomic"]],
            cds_g_start=t["cds_g_start"],
            cds_g_end=t["cds_g_end"],
            is_canonical=True,
        )
        txs[tx_id] = tx
        txs[t["gene_symbol"]] = tx

    out = StaticProvider(txs, domains={}, knowledge={})
    out.assembly = fixture_data["assembly"]
    return out


@pytest.mark.parametrize(
    "case_id",
    [
        "tmprss2_erg_zero_coding",
        "ewsr1_erg_inframe",
        "egfr_viii_inframe",
        "egfr_rad51_inframe",
    ],
)
def test_selected_cbioportal_cases_annotate_from_genomic_breakpoints(provider, fixture_data, case_id):
    case = next(c for c in fixture_data["cases"] if c["case_id"] == case_id)
    result = annotate_fusion(
        provider,
        case["five_gene"],
        case["three_gene"],
        five_tx=case["five_transcript"],
        three_tx=case["three_transcript"],
        five_genomic=case["five_genomic"],
        three_genomic=case["three_genomic"],
    )

    assert result["resolved"]["genome_build"] == "GRCh37"
    assert result["resolved"]["five"]["breakpoint"]["type"] == "genomic"
    assert result["resolved"]["three"]["breakpoint"]["type"] == "genomic"
    assert result["resolved"]["five"]["transcript"] == case["five_transcript"]
    assert result["resolved"]["three"]["transcript"] == case["three_transcript"]
    assert result["interface"]["frame_status"] == case["expected_frame_status"]
    assert result["interface"]["hgvsp_like"] == case["expected_hgvsp_like"]


def test_tmprss2_erg_zero_coding_case_keeps_zero_aa_five_prime_partner(provider, fixture_data):
    case = next(c for c in fixture_data["cases"] if c["case_id"] == "tmprss2_erg_zero_coding")
    result = annotate_fusion(
        provider,
        case["five_gene"],
        case["three_gene"],
        five_tx=case["five_transcript"],
        three_tx=case["three_transcript"],
        five_genomic=case["five_genomic"],
        three_genomic=case["three_genomic"],
    )
    assert result["interface"]["five_last_aa"] == 0
    assert result["interface"]["hgvsp_like"].startswith("TMPRSS2:p.0::")
    assert result["resolved"]["five"]["breakpoint"]["context"]["region"] == "intron"
    assert "5' UTR" in result["resolved"]["five"]["breakpoint"]["context"]["label"]
    assert result["resolved"]["five"]["structure"]["exons"]


def test_one_sided_egfr_case_still_raises_clear_mapping_error(provider, fixture_data):
    case = next(c for c in fixture_data["cases"] if c["case_id"] == "egfr_one_sided_unresolved")
    with pytest.raises(ValueError, match=case["expected_error"]):
        annotate_fusion(
            provider,
            case["five_gene"],
            case["three_gene"],
            five_tx=case["five_transcript"],
            three_tx=case["three_transcript"],
            five_genomic=case["five_genomic"],
            three_genomic=case["three_genomic"],
        )
