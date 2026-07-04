import type { AnnotateParams } from "./types";

// Curated, known-good fusion presets for the "Try an example" row. Each pair
// (and exon combination) is verified against real annotation output — see
// docs/FRAME_VALIDATION.md for the NPM1::ALK / LMNA::NTRK1 exon-search
// results, and README.md / examples/*.py for EML4::ALK — so these are safe
// to run against the live deployed API without risking a confusing dead end.
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
];
