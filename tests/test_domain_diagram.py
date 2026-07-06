"""Unit tests for transcript-structure diagram helpers."""
from fusion_annotation.domain_diagram import _transcript_breakpoint_x


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
