"""Offline validation of the annotation engines against the EML4::ALK v1 (E13;A20)
worked example. Every truth value here was verified against Ensembl primary data
(GRCh38, canonical transcripts) and reconstructed exactly. No network required.
"""
import json
import os
import pytest

from fusion_annotation import (
    Transcript, build_exon_cds_map, annotate_fusion, annotate_effect,
    cds_coord_at_exon_boundary,
)
from fusion_annotation.providers import StaticProvider

FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "eml4_alk_fixture.json")


@pytest.fixture(scope="module")
def provider():
    fx = json.load(open(FIXTURE))
    txs = {}
    for key, t in fx["transcripts"].items():
        exon_cds = build_exon_cds_map(t["strand"], t["exons"], t["cds_g_start"], t["cds_g_end"])
        txs[key] = Transcript(
            gene_symbol=t["gene_symbol"], gene_id=t["gene_id"],
            transcript_id=t["transcript_id"], strand=t["strand"],
            cds=t["cds"], protein=t["protein"], uniprot=t["uniprot"], exon_cds=exon_cds)
    return StaticProvider(txs, domains={k: v for k, v in fx["domains"].items()},
                          knowledge=fx["knowledge"])


@pytest.fixture(scope="module")
def result(provider):
    return annotate_fusion(provider, "EML4", "ALK", five_exon=13, three_exon=20)


# ---- transcript / CDS sanity (matches Ensembl exactly) ---------------------
def test_parent_cds_lengths(provider):
    assert provider.get_transcript("EML4").cds_len() == 2946
    assert provider.get_transcript("ALK").cds_len() == 4863


def test_breakpoint_cds_coords(provider):
    eml4 = provider.get_transcript("EML4")
    alk = provider.get_transcript("ALK")
    assert cds_coord_at_exon_boundary(eml4.exon_cds, 13, "end") == 1489
    assert cds_coord_at_exon_boundary(alk.exon_cds, 20, "start") == 3173


# ---- Layer 1: effect -------------------------------------------------------
def test_fusion_is_in_frame_no_internal_stops(result):
    iface = result["interface"]
    # fusion CDS = 496*3 + 1 (hybrid nt) + remainder of ALK = 3180 nt -> 1060 codons
    assert len(iface["fusion_protein_seq"]) * 3 + 3 == 3180  # +3 for the stop codon
    assert iface["in_frame"] is True
    assert iface["internal_stops"] == 0
    assert iface["frame_status"] == "in-frame"


def test_fusion_protein_length(result):
    # 3180 nt / 3 = 1060 codons; last is stop -> 1059 aa mature protein
    assert result["interface"]["fusion_length"] == 1059


def test_hybrid_junction_codon(result):
    # EML4 contributes 496 complete codons (..Lys496); junction codon is a hybrid
    # of 1 nt EML4 (G) + 2 nt ALK (TG) = GTG = Val.
    iface = result["interface"]
    assert iface["hybrid_codon"] is True
    assert iface["junction_residue"] == "V"       # Val
    assert iface["five_last_aa"] == 496
    assert iface["five_last_aa_res"] == "K"        # Lys496
    assert iface["three_first_aa"] == 1059         # first fully-retained ALK residue (Tyr1059)
    assert iface["three_first_aa_res"] == "Y"


# ---- Layer 2: HGVS.p-like interface ---------------------------------------
def test_hgvsp_string(result):
    hgvsp = result["interface"]["hgvsp_like"]
    assert hgvsp.startswith("EML4:p.Met1_Lys496::ALK:p.Tyr1059_")
    assert "::" in hgvsp
    assert "hybrid codon -> Val" in hgvsp


def test_categorical_key(result):
    assert result["interface"]["categorical_key"] == "EML4::ALK"


# ---- domain retention ------------------------------------------------------
def _by_status(iface, status):
    return [d for d in iface["domains"] if d["status"] == status]


def test_alk_kinase_retained(result):
    retained = " ".join(d["name"].lower() for d in _by_status(result["interface"], "RETAINED"))
    assert "kinase" in retained


def test_eml4_second_betapropeller_lost(result):
    # EML4's C-terminal (second) β-propeller sits 3' of the aa496 breakpoint -> LOST.
    lost = " ".join(d["name"].lower() for d in _by_status(result["interface"], "LOST"))
    assert "second beta-propeller" in lost


def test_eml4_help_motif_retained(result):
    # The HELP motif (aa 255-293) is 5' of the breakpoint -> RETAINED.
    retained = " ".join(d["name"].lower() for d in _by_status(result["interface"], "RETAINED"))
    assert "help motif" in retained


# ---- Layer 3: knowledge (OncoKB-like, from fixture) ------------------------
def test_knowledge_oncogenic_and_therapies(result):
    kn = result["knowledge"]
    assert kn["oncogenic"] == "Oncogenic"
    thx = " ".join(kn["therapies"]).lower()
    assert "crizotinib" in thx
