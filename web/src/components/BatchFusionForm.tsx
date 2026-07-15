import { useState, type FormEvent } from "react";
import type { AnnotateParams } from "../lib/types";

interface Props {
  genomeBuild: string;
  onSubmit: (params: AnnotateParams[]) => void;
  onCurate?: (params: AnnotateParams[]) => void;
  loading: boolean;
  curationLoading?: boolean;
  curationEnabled?: boolean;
}

function parseBatchLine(line: string, genomeBuild: string): AnnotateParams | null {
  const trimmed = line.trim();
  if (!trimmed || trimmed.startsWith("#")) return null;

  const [fusionPart, fiveExon, threeExon] = trimmed
    .split(/[,\t]/)
    .map((part) => part.trim());
  const [fiveGene, threeGene] = fusionPart.split(/::|--/).map((part) => part.trim());
  if (!fiveGene || !threeGene) {
    throw new Error(`Could not parse fusion partners from line: ${line}`);
  }
  if (!fiveExon || !threeExon) {
    throw new Error(`Batch exon mode requires five_exon and three_exon: ${line}`);
  }
  if (!/^\d+$/.test(fiveExon) || !/^\d+$/.test(threeExon)) {
    throw new Error(`Exon numbers must be positive integers: ${line}`);
  }
  return {
    five_gene: fiveGene,
    three_gene: threeGene,
    five_exon: fiveExon,
    three_exon: threeExon,
    genome_build: genomeBuild,
    input_mode: "exon",
  };
}

export function BatchFusionForm({
  genomeBuild,
  onSubmit,
  onCurate,
  loading,
  curationLoading = false,
  curationEnabled = false,
}: Props) {
  const [text, setText] = useState("EML4::ALK,13,20\nBCR::ABL1,13,2");
  const [parseError, setParseError] = useState<string | null>(null);

  function parseInput(): AnnotateParams[] {
    const parsed = text
      .split(/\r?\n/)
      .map((line) => parseBatchLine(line, genomeBuild))
      .filter((item): item is AnnotateParams => item !== null);
    if (!parsed.length) {
      throw new Error("Add at least one fusion line before running a batch.");
    }
    setParseError(null);
    return parsed;
  }

  function handleSubmit(event: FormEvent) {
    event.preventDefault();
    try {
      onSubmit(parseInput());
    } catch (error) {
      setParseError(error instanceof Error ? error.message : String(error));
    }
  }

  function handleCurate() {
    if (!onCurate) return;
    try {
      onCurate(parseInput());
    } catch (error) {
      setParseError(error instanceof Error ? error.message : String(error));
    }
  }

  return (
    <form className="batch-form" onSubmit={handleSubmit}>
      <div className="batch-form-header">
        <div>
          <h2>Batch annotation</h2>
          <p>Run multiple exon-defined fusions against the same Genome Nexus workflow.</p>
        </div>
        <span className="batch-format">GENE1::GENE2, five_exon, three_exon</span>
      </div>
      <textarea
        value={text}
        onChange={(event) => setText(event.target.value)}
        spellCheck={false}
        rows={5}
        aria-label="Batch fusion input"
      />
      {parseError && <div className="error-box">{parseError}</div>}
      <div className="batch-actions">
        <button type="submit" disabled={loading || curationLoading} className="submit-button">
          {loading ? "Annotating batch..." : "Annotate batch"}
        </button>
        {onCurate && (
          <button
            type="button"
            disabled={!curationEnabled || loading || curationLoading}
            className="secondary-button"
            onClick={handleCurate}
          >
            {curationLoading ? "Curating genes..." : "Curate genes"}
          </button>
        )}
      </div>
    </form>
  );
}
