import { useEffect, useState, type FormEvent } from "react";
import type { AnnotateParams } from "../lib/types";

interface Props {
  initial: AnnotateParams;
  onSubmit: (params: AnnotateParams) => void;
  loading: boolean;
}

/** Lookup form for the fusion inputs. Breakpoints can be given as an exon
 * number or a genomic position per partner (mirrors the annotate_gene_fusion
 * MCP tool / AnnotateRequest schema in api/app.py). */
export function FusionForm({ initial, onSubmit, loading }: Props) {
  const [values, setValues] = useState<AnnotateParams>(initial);
  const [showAdvanced, setShowAdvanced] = useState(
    Boolean(initial.five_transcript || initial.three_transcript || initial.species !== "homo_sapiens"),
  );

  // Resync when `initial` changes after mount — e.g. browser Back/Forward
  // (App.tsx's popstate handler) updates the URL and passes a new `initial`
  // down, but this component's own state wouldn't otherwise pick it up,
  // leaving the form stale relative to the displayed result.
  useEffect(() => {
    setValues(initial);
    setShowAdvanced(
      Boolean(initial.five_transcript || initial.three_transcript || initial.species !== "homo_sapiens"),
    );
  }, [initial]);

  function update<K extends keyof AnnotateParams>(key: K, value: AnnotateParams[K]) {
    setValues((prev) => ({ ...prev, [key]: value }));
  }

  function handleSubmit(e: FormEvent) {
    e.preventDefault();
    onSubmit(values);
  }

  return (
    <form className="fusion-form" onSubmit={handleSubmit}>
      <div className="form-row">
        <label>
          5' gene
          <input
            required
            value={values.five_gene}
            onChange={(e) => update("five_gene", e.target.value)}
            placeholder="EML4"
          />
        </label>
        <label>
          3' gene
          <input
            required
            value={values.three_gene}
            onChange={(e) => update("three_gene", e.target.value)}
            placeholder="ALK"
          />
        </label>
      </div>

      <div className="form-row">
        <label>
          5' last exon
          <input
            value={values.five_exon ?? ""}
            onChange={(e) => update("five_exon", e.target.value)}
            placeholder="13"
            inputMode="numeric"
          />
        </label>
        <label>
          3' first exon
          <input
            value={values.three_exon ?? ""}
            onChange={(e) => update("three_exon", e.target.value)}
            placeholder="20"
            inputMode="numeric"
          />
        </label>
      </div>

      <div className="form-row">
        <label>
          5' genomic breakpoint <span className="hint">(instead of exon)</span>
          <input
            value={values.five_genomic ?? ""}
            onChange={(e) => update("five_genomic", e.target.value)}
            placeholder="chr2:42295516"
          />
        </label>
        <label>
          3' genomic breakpoint <span className="hint">(instead of exon)</span>
          <input
            value={values.three_genomic ?? ""}
            onChange={(e) => update("three_genomic", e.target.value)}
            placeholder="chr2:29223528"
          />
        </label>
      </div>

      <div className="form-row">
        <label>
          Genome build
          <select value={values.genome_build} onChange={(e) => update("genome_build", e.target.value)}>
            <option value="GRCh38">GRCh38</option>
            <option value="GRCh37">GRCh37</option>
          </select>
        </label>
      </div>

      <button
        type="button"
        className="advanced-toggle"
        onClick={() => setShowAdvanced((s) => !s)}
      >
        {showAdvanced ? "Hide" : "Show"} advanced options
      </button>

      {showAdvanced && (
        <div className="form-row">
          <label>
            5' transcript
            <input
              value={values.five_transcript ?? ""}
              onChange={(e) => update("five_transcript", e.target.value)}
              placeholder="ENST00000318522 (default: canonical)"
            />
          </label>
          <label>
            3' transcript
            <input
              value={values.three_transcript ?? ""}
              onChange={(e) => update("three_transcript", e.target.value)}
              placeholder="ENST00000389048 (default: canonical)"
            />
          </label>
          <label>
            Species
            <input value={values.species} onChange={(e) => update("species", e.target.value)} />
          </label>
        </div>
      )}

      <button type="submit" disabled={loading} className="submit-button">
        {loading ? "Annotating…" : "Annotate fusion"}
      </button>
    </form>
  );
}
