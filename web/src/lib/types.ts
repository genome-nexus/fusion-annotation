// Mirrors the JSON shape returned by GET/POST /api/annotate (api/app.py),
// which is itself annotate_fusion()'s return value
// (src/fusion_annotation/core.py): {"interface", "knowledge", "resolved", "warnings"}.

export type DomainStatus = "RETAINED" | "LOST" | "DISRUPTED" | "UNKNOWN";
export type BreakpointContextRegion =
  | "upstream"
  | "utr5"
  | "coding"
  | "utr3"
  | "intron"
  | "exon_boundary"
  | "downstream"
  | "unknown";

export interface DomainCall {
  accession: string;
  name: string;
  type: string;
  start: number;
  end: number;
  status: DomainStatus;
  gene: string;
  partner_protein_length: number;
}

export interface TranscriptStructureSegment {
  kind: "utr5" | "coding" | "utr3";
  start: number;
  end: number;
}

export interface TranscriptStructureExon {
  rank: number;
  genomic_start: number;
  genomic_end: number;
  length: number;
  segments: TranscriptStructureSegment[];
}

export interface TranscriptStructure {
  strand: 1 | -1;
  chrom: string | null;
  promoter_window_bp: number;
  tss_genomic: number;
  transcript_end_genomic: number;
  transcript_length: number;
  exons: TranscriptStructureExon[];
}

export interface BreakpointContext {
  region: BreakpointContextRegion;
  label: string;
  exon_rank: number | null;
  intron_rank: number | null;
  exon_offset: number | null;
  exon_length: number | null;
  boundary: "before" | "after" | null;
}

export interface FusionInterface {
  five_gene: string;
  three_gene: string;
  five_transcript: string;
  three_transcript: string;
  five_last_aa: number | null;
  three_first_aa: number | null;
  five_last_aa_res: string | null;
  three_first_aa_res: string | null;
  in_frame: boolean | null;
  hybrid_codon: boolean;
  junction_residue: string | null;
  fusion_length: number | null;
  internal_stops: number | null;
  fusion_protein_seq: string;
  frame_status: "in-frame" | "out-of-frame" | "frameshift-truncating" | "unknown";
  domains: DomainCall[];
  breakpoint_label: string | null;
  categorical_key: string;
  hgvsp_like: string;
}

export interface FusionKnowledge {
  categorical_key: string;
  oncogenic: string | null;
  therapies: string[];
  evidence: Array<Record<string, unknown>>;
  diseases: string[];
  sources: string[];
}

export interface ResolvedPartner {
  gene: string;
  transcript: string;
  transcript_source: "user-specified" | "canonical" | "non-canonical" | "provider-default";
  breakpoint: {
    type: "genomic" | "exon" | "unknown";
    exon?: number;
    genomic_position?: number;
    cds_coord: number | null;
    context: BreakpointContext;
  };
  /** Full-length (untruncated) protein size for this partner. Used to lay
   * out the domain diagram's backbone even when the partner has zero
   * annotated domains. */
  protein_length: number;
  structure: TranscriptStructure | null;
}

export interface AnnotationResult {
  interface: FusionInterface;
  knowledge: FusionKnowledge;
  resolved: {
    genome_build: string;
    five: ResolvedPartner;
    three: ResolvedPartner;
  };
  warnings: string[];
}

export interface BatchAnnotationItemResult {
  input: AnnotateParams;
  result: AnnotationResult | null;
  error: string | null;
}

export interface BatchAnnotationResponse {
  results: BatchAnnotationItemResult[];
}

export interface GeneCurationStatus {
  enabled: boolean;
  model: string;
}

export interface GeneFusionCurationContext {
  gene: string;
  fusion: string;
  side: "five_prime" | "three_prime";
  partner_gene: string;
  five_transcript?: string | null;
  three_transcript?: string | null;
  five_exon?: string | null;
  three_exon?: string | null;
  five_genomic?: string | null;
  three_genomic?: string | null;
  five_protein_breakpoint?: string | null;
  three_protein_breakpoint?: string | null;
  retained_domains?: string[];
  lost_domains?: string[];
  disrupted_domains?: string[];
  kinase_gene?: string | null;
  kinase_gene_side?: "five_prime" | "three_prime" | null;
  kinase_domain_status?: "retained" | "lost" | "disrupted" | "unknown" | null;
  annotation_error?: string | null;
}

export interface GeneCurationGeneResult {
  gene: string;
  cancer_associated?: boolean | null;
  rationale?: string;
  supporting_pmids?: string[];
  retrieved_pmids?: string[];
  fusion_contexts?: GeneFusionCurationContext[];
  insufficient_evidence?: boolean;
  error?: string;
}

export interface GeneCurationResponse {
  genes: GeneCurationGeneResult[];
}

// The inputs a caller supplies — mirrors AnnotateRequest in api/app.py, and is
// exactly what gets encoded into the permalink's URL query string.
export interface AnnotateParams {
  five_gene: string;
  three_gene: string;
  five_exon?: string;
  three_exon?: string;
  five_genomic?: string;
  three_genomic?: string;
  five_transcript?: string;
  three_transcript?: string;
  genome_build: string;
  variant_type?: string;
  input_mode?: "exon" | "genomic";
}

export interface ApiError {
  status: number;
  detail: string;
}
