"""Unit tests for domain-diagram helpers."""
from fusion_annotation.domain_diagram import _transcript_breakpoint_x, canonicalize_domains


def test_canonicalize_domains_keeps_pfam_only_regions_when_no_named_overlap():
    raw = [
        {
            "accession": "PF02198",
            "name": "PF02198",
            "type": "domain",
            "start": 123,
            "end": 204,
            "status": "LOST",
            "gene": "ERG",
        },
        {
            "accession": "PF00178",
            "name": "PF00178",
            "type": "domain",
            "start": 317,
            "end": 399,
            "status": "RETAINED",
            "gene": "ERG",
        },
    ]

    assert canonicalize_domains(raw, "ERG") == [
        {"name": "PF02198", "start": 123, "end": 204, "status": "LOST"},
        {"name": "PF00178", "start": 317, "end": 399, "status": "RETAINED"},
    ]


def test_canonicalize_domains_prefers_curated_name_over_accession_placeholder():
    raw = [
        {
            "accession": "PF00178",
            "name": "PF00178",
            "type": "domain",
            "start": 317,
            "end": 399,
            "status": "RETAINED",
            "gene": "ERG",
        },
        {
            "accession": "IPR000837",
            "name": "ETS domain",
            "type": "domain",
            "start": 316,
            "end": 398,
            "status": "RETAINED",
            "gene": "ERG",
        },
    ]

    assert canonicalize_domains(raw, "ERG") == [
        {"name": "ETS domain", "start": 316, "end": 399, "status": "RETAINED"},
    ]


def test_transcript_breakpoint_x_returns_none_for_out_of_range_exon_rank():
    partner = {
        "breakpoint": {
            "context": {
                "region": "coding",
                "exon_rank": 0,
                "intron_rank": None,
                "boundary": None,
                "exon_offset": 1,
                "exon_length": 10,
            }
        }
    }
    layout = {
        "width": 100.0,
        "promoter_start": 0.0,
        "promoter_end": 10.0,
        "exons": [{"start": 20.0, "end": 30.0, "width": 10.0, "length": 10}],
    }

    assert _transcript_breakpoint_x(partner, layout) is None


def test_transcript_breakpoint_x_returns_none_for_out_of_range_intron_rank():
    partner = {
        "breakpoint": {
            "context": {
                "region": "intron",
                "exon_rank": None,
                "intron_rank": 2,
                "boundary": None,
                "exon_offset": None,
                "exon_length": None,
            }
        }
    }
    layout = {
        "width": 100.0,
        "promoter_start": 0.0,
        "promoter_end": 10.0,
        "exons": [{"start": 20.0, "end": 30.0, "width": 10.0, "length": 10}],
    }

    assert _transcript_breakpoint_x(partner, layout) is None
