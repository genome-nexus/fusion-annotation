/** Generate external links to public databases for fusion annotation sources. */

export interface ExternalLink {
  url: string;
  title: string;
}

/** Pfam domain accession → Pfam database link. */
export function pfamLink(accession: string): ExternalLink {
  return {
    url: `https://pfam.xfam.org/family/${accession}`,
    title: `Pfam: ${accession}`,
  };
}

/** Ensembl transcript ID → Ensembl database link. */
export function ensemblTranscriptLink(transcriptId: string, genomeBuild: string = "GRCh38"): ExternalLink {
  const species = "homo_sapiens";
  const version = genomeBuild === "GRCh37" ? "?db=core" : "";
  return {
    url: `https://www.ensembl.org/${species}/Transcript/Summary?t=${transcriptId}${version}`,
    title: `Ensembl: ${transcriptId}`,
  };
}

/** Gene symbol → Ensembl gene link. */
export function ensemblGeneLink(geneSymbol: string, genomeBuild: string = "GRCh38"): ExternalLink {
  const species = "homo_sapiens";
  const version = genomeBuild === "GRCh37" ? "?db=core" : "";
  return {
    url: `https://www.ensembl.org/${species}/Search/Query?q=${geneSymbol}${version}`,
    title: `Ensembl: ${geneSymbol}`,
  };
}

/** CIViC evidence ID → CIViC evidence link. */
export function civicEvidenceLink(evidenceId: number | string): ExternalLink {
  return {
    url: `https://civicdb.org/evidence/${evidenceId}/summary`,
    title: `CIViC Evidence #${evidenceId}`,
  };
}

/** CIViC molecular profile ID → CIViC molecular profile link. */
export function civicMolecularProfileLink(mpId: number | string): ExternalLink {
  return {
    url: `https://civicdb.org/molecular-profiles/${mpId}/summary`,
    title: `CIViC Molecular Profile #${mpId}`,
  };
}

/** UCSC genomic position → UCSC Genome Browser link. */
export function ucscLink(chrom: string, pos: number, genomeBuild: string = "GRCh38"): ExternalLink {
  const db = genomeBuild === "GRCh37" ? "hg19" : "hg38";
  const region = `${chrom}:${pos - 100}-${pos + 100}`;
  return {
    url: `https://genome.ucsc.edu/cgi-bin/hgTracks?db=${db}&position=${region}`,
    title: `UCSC Genome Browser`,
  };
}

/** Genomic breakpoint string "chr7:55268799" → UCSC link. */
export function genomicBreakpointUcscLink(genomicStr: string | null, genomeBuild: string = "GRCh38"): ExternalLink | null {
  if (!genomicStr) return null;
  const match = genomicStr.match(/^(chr[\dMXY]+):(\d+)$/);
  if (!match) return null;
  const [, chrom, posStr] = match;
  const pos = parseInt(posStr, 10);
  return ucscLink(chrom, pos, genomeBuild);
}

/** OncoKB gene fusion knowledge base link. */
export function oncokbFusionLink(geneA: string, geneB: string): ExternalLink {
  return {
    url: `https://www.oncokb.org/gene/${geneA}/${geneB}`,
    title: `OncoKB: ${geneA}::${geneB}`,
  };
}

/** InterPro protein domain link. */
export function interproLink(accession: string): ExternalLink {
  return {
    url: `https://www.ebi.ac.uk/interpro/entry/${accession}`,
    title: `InterPro: ${accession}`,
  };
}
