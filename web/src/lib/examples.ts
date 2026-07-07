import type { AnnotateParams } from "./types";

// Curated, known-good fusion presets for the "Try an example" row. Each case
// is verified against real annotation output — see docs/FRAME_VALIDATION.md,
// README.md / examples/*.py, and tests/test_cbioportal_representatives.py — so
// these are safe to run against the live deployed API without risking a
// confusing dead end.
export interface FusionExample {
  label: string;
  description: string;
  params: Partial<AnnotateParams>;
}

export const EXAMPLES: FusionExample[] = [
  {
    label: "EML4::ALK",
    description: "NSCLC driver, exon breakpoints (E13;A20)",
    params: { five_gene: "EML4", three_gene: "ALK", five_exon: "13", three_exon: "20" },
  },
  {
    label: "EML4::ALK (genomic)",
    description: "same fusion, specified by genomic breakpoints instead",
    params: {
      five_gene: "EML4",
      three_gene: "ALK",
      five_genomic: "chr2:42295516",
      three_genomic: "chr2:29223528",
      variant_type: "inversion",
    },
  },
  {
    label: "NPM1::ALK",
    description: "anaplastic large-cell lymphoma, t(2;5) (E4;A20)",
    params: { five_gene: "NPM1", three_gene: "ALK", five_exon: "4", three_exon: "20" },
  },
  {
    label: "LMNA::NTRK1",
    description: "TRK-inhibitor-responsive fusion (E2;A11)",
    params: { five_gene: "LMNA", three_gene: "NTRK1", five_exon: "2", three_exon: "11" },
  },
  {
    label: "TMPRSS2::ERG",
    description: "cBioPortal representative, promoter/5'UTR breakpoint (GRCh37 genomic)",
    params: {
      five_gene: "TMPRSS2",
      three_gene: "ERG",
      five_genomic: "chr21:42876132",
      three_genomic: "chr21:39822160",
      five_transcript: "ENST00000332149",
      three_transcript: "ENST00000398919",
      genome_build: "GRCh37",
      variant_type: "deletion",
    },
  },
  {
    label: "EWSR1::ERG",
    description: "cBioPortal representative, in-frame translocation (GRCh37 genomic)",
    params: {
      five_gene: "EWSR1",
      three_gene: "ERG",
      five_genomic: "chr22:29683558",
      three_genomic: "chr21:39759195",
      five_transcript: "ENST00000397938",
      three_transcript: "ENST00000398919",
      genome_build: "GRCh37",
      variant_type: "translocation",
    },
  },
  {
    label: "EGFRvIII-like",
    description: "cBioPortal representative, in-frame EGFR intragenic deletion (GRCh37 genomic)",
    params: {
      five_gene: "EGFR",
      three_gene: "EGFR",
      five_genomic: "chr7:55092872",
      three_genomic: "chr7:55223428",
      five_transcript: "ENST00000275493",
      three_transcript: "ENST00000275493",
      genome_build: "GRCh37",
      variant_type: "deletion",
    },
  },
  {
    label: "EGFR::RAD51",
    description: "cBioPortal representative, in-frame translocation (GRCh37 genomic)",
    params: {
      five_gene: "EGFR",
      three_gene: "RAD51",
      five_genomic: "chr7:55268799",
      three_genomic: "chr15:40996189",
      five_transcript: "ENST00000275493",
      three_transcript: "ENST00000267868",
      genome_build: "GRCh37",
      variant_type: "translocation",
    },
  },
];
