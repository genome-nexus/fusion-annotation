"""Offline tests for issue #3: genomic-coordinate breakpoints, transcript echo, and
the known-oncogenic-pair out-of-frame sanity flag. Uses the same EML4::ALK fixture
as test_eml4_alk.py, whose coordinates are verified against Ensembl GRCh38.
"""
import json
import os

import pytest

from fusion_annotation import (
    Transcript, build_exon_cds_map, build_exon_genomic_map, annotate_fusion,
    cds_coord_at_exon_boundary, cds_coord_at_genomic, parse_genomic_breakpoint,
)
from fusion_annotation.providers import StaticProvider

FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "eml4_alk_fixture.json")


@pytest.fixture(scope="module")
def provider():
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
    # Also key each transcript under its id so an explicit five_tx/three_tx resolves.
    for t in list(txs.values()):
        txs[t.transcript_id] = t
    return StaticProvider(txs, domains=fx["domains"], knowledge=fx["knowledge"])


# ---- genomic breakpoint string parsing ------------------------------------
@pytest.mark.parametrize("value,expected", [
    (117324415, 117324415),
    ("117324415", 117324415),
    ("chr6:117324415", 117324415),
    ("6:117324415", 117324415),
    ("g.117324415", 117324415),
    ("chr6:g.117324415", 117324415),
    ("chr6:g.117,324,415", 117324415),
])
def test_parse_genomic_breakpoint_forms(value, expected):
    assert parse_genomic_breakpoint(value) == expected


def test_parse_genomic_breakpoint_bad():
    with pytest.raises(ValueError):
        parse_genomic_breakpoint("not-a-position")


# ---- genomic->CDS mapping agrees with the exon-boundary path --------------
def test_genomic_maps_to_same_cds_as_exon_boundary(provider):
    """The exon-13/exon-20 boundaries have genomic positions; resolving those
    positions must land on the identical CDS coords the exon path produces."""
    eml4 = provider.get_transcript("EML4")
    alk = provider.get_transcript("ALK")

    # EML4 is +strand: exon 13's 3' end is the exon's larger genomic coordinate.
    g_five = eml4.exon_genomic[12][1]
    # ALK is -strand: exon 20's 5' start is the exon's larger genomic coordinate.
    g_three = alk.exon_genomic[19][1]

    assert cds_coord_at_genomic(eml4.strand, eml4.exon_genomic, eml4.cds_g_start,
                                eml4.cds_g_end, g_five, "end") == \
        cds_coord_at_exon_boundary(eml4.exon_cds, 13, "end")
    assert cds_coord_at_genomic(alk.strand, alk.exon_genomic, alk.cds_g_start,
                                alk.cds_g_end, g_three, "start") == \
        cds_coord_at_exon_boundary(alk.exon_cds, 20, "start")


def test_genomic_input_reproduces_exon_annotation(provider):
    eml4 = provider.get_transcript("EML4")
    alk = provider.get_transcript("ALK")
    by_exon = annotate_fusion(provider, "EML4", "ALK", five_exon=13, three_exon=20)
    by_genomic = annotate_fusion(
        provider, "EML4", "ALK",
        five_genomic=eml4.exon_genomic[12][1],
        three_genomic="chr2:g.%d" % alk.exon_genomic[19][1])
    assert by_genomic["interface"]["hgvsp_like"] == by_exon["interface"]["hgvsp_like"]
    assert by_genomic["interface"]["in_frame"] is True
    assert by_genomic["resolved"]["five"]["breakpoint"]["type"] == "genomic"
    assert by_genomic["resolved"]["three"]["breakpoint"]["type"] == "genomic"


def test_genomic_position_outside_cds_raises(provider):
    # A position 3' of the whole CDS has no downstream coding base to start from.
    # ALK is -strand, so 3' of the CDS is *below* the genomic lo bound.
    alk = provider.get_transcript("ALK")
    with pytest.raises(ValueError):
        cds_coord_at_genomic(alk.strand, alk.exon_genomic, alk.cds_g_start,
                             alk.cds_g_end, alk.cds_g_start - 10_000_000, "start")


# ---- transcript echo / pinning --------------------------------------------
def test_resolved_echoes_transcript_and_source(provider):
    r = annotate_fusion(provider, "EML4", "ALK", five_exon=13, three_exon=20)
    assert r["resolved"]["five"]["transcript"] == "ENST00000318522"
    assert r["resolved"]["five"]["transcript_source"] == "canonical"
    assert r["resolved"]["five"]["breakpoint"] == {"type": "exon", "exon": 13, "cds_coord": 1489}


def test_resolved_marks_user_specified_transcript(provider):
    # Passing the transcript id explicitly should be echoed as user-specified.
    r = annotate_fusion(provider, "EML4", "ALK", five_exon=13, three_exon=20,
                        five_tx="ENST00000318522")
    assert r["resolved"]["five"]["transcript_source"] == "user-specified"


# ---- sanity flag on a known oncogenic pair --------------------------------
def test_known_pair_out_of_frame_warns(provider):
    # E13;A29 is frameshift-truncating (see README) — a known driver pair coming
    # back out-of-frame should raise the transcript/exon re-check warning.
    r = annotate_fusion(provider, "EML4", "ALK", five_exon=13, three_exon=29)
    assert r["interface"]["frame_status"] != "in-frame"
    assert any("known oncogenic fusion pair" in w for w in r["warnings"])
    assert any("genomic breakpoints" in w for w in r["warnings"])


def test_known_pair_in_frame_no_warning(provider):
    r = annotate_fusion(provider, "EML4", "ALK", five_exon=13, three_exon=20)
    assert r["warnings"] == []


def test_out_of_frame_warning_is_context_aware_for_genomic(provider):
    # When genomic breakpoints were already supplied, the warning must not tell the
    # caller to "supply genomic breakpoints" — those already pin the isoform.
    eml4 = provider.get_transcript("EML4")
    alk = provider.get_transcript("ALK")
    r = annotate_fusion(provider, "EML4", "ALK",
                        five_genomic=eml4.exon_genomic[12][1],    # EML4 exon 13 3' end
                        three_genomic=alk.exon_genomic[28][1])    # ALK exon 29 5' start
    assert r["interface"]["frame_status"] != "in-frame"
    assert len(r["warnings"]) == 1
    warning = r["warnings"][0]
    assert "known oncogenic fusion pair" in warning
    assert "supply genomic breakpoints" not in warning
    assert "five_genomic/three_genomic" not in warning


# ---- missing breakpoint ----------------------------------------------------
def test_missing_breakpoint_raises(provider):
    with pytest.raises(ValueError):
        annotate_fusion(provider, "EML4", "ALK", five_exon=13)
