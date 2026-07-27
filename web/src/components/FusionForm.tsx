import { useEffect, useState, type FormEvent } from "react";
import type { AnnotateParams } from "../lib/types";
import type { DerivedInputs } from "../lib/derivedInputs";
import { derivedToParams } from "../lib/derivedInputs";

interface Props {
  initial: AnnotateParams;
  derived: DerivedInputs | null;
  onSubmit: (params: AnnotateParams) => void;
  loading: boolean;
}

const VARIANT_TYPES = [
  { value: "translocation", label: "Translocation" },
  { value: "deletion", label: "Deletion / internal deletion" },
  { value: "inversion", label: "Inversion" },
  { value: "duplication", label: "Tandem duplication" },
  { value: "insertion", label: "Insertion" },
];

/** Lookup form for the fusion inputs. Breakpoints can be given as an exon
 * number or a genomic position per partner (mirrors the annotate_gene_fusion
 * MCP tool / AnnotateRequest schema in api/app.py). */
export function FusionForm({ initial, derived, onSubmit, loading }: Props) {
  const initialMode = (initial.input_mode ??
    (initial.five_genomic || initial.three_genomic ? "genomic" : "exon")) as "exon" | "genomic";
  const [mode, setMode] = useState<"exon" | "genomic">(initialMode);
  const [values, setValues] = useState<AnnotateParams>({ ...initial, input_mode: initialMode });

  useEffect(() => {
    const newMode = (initial.input_mode ??
      (initial.five_genomic || initial.three_genomic ? "genomic" : "exon")) as "exon" | "genomic";
    setMode(newMode);
    setValues({ ...initial, input_mode: newMode });
  }, [initial]);

  function update<K extends keyof AnnotateParams>(key: K, value: AnnotateParams[K]) {
    setValues((prev) => ({ ...prev, [key]: value }));
  }

  function switchMode(newMode: "exon" | "genomic") {
    setMode(newMode);
    if (derived) {
      setValues(derivedToParams(derived, newMode, values));
    } else {
      setValues((prev) => ({ ...prev, input_mode: newMode }));
    }
  }

  function handleSubmit(e: FormEvent) {
    e.preventDefault();
    onSubmit({ ...values, input_mode: mode });
  }

  return (
    <form className="fusion-form" onSubmit={handleSubmit}>
      {/* Mode tabs */}
      <div className="input-mode-tabs" role="tablist" aria-label="Input mode">
        <button
          type="button"
          role="tab"
          aria-selected={mode === "exon"}
          className={mode === "exon" ? "mode-tab active" : "mode-tab"}
          onClick={() => switchMode("exon")}
        >
          By exon
        </button>
        <button
          type="button"
          role="tab"
          aria-selected={mode === "genomic"}
          className={mode === "genomic" ? "mode-tab active" : "mode-tab"}
          onClick={() => switchMode("genomic")}
        >
          By genomic position
        </button>
      </div>

      {/* Gene inputs — always shown */}
      <div className="form-row">
        <label>
          5′ gene
          <input
            required
            value={values.five_gene}
            onChange={(e) => update("five_gene", e.target.value)}
            placeholder="EML4"
          />
        </label>
        <label>
          3′ gene
          <input
            required
            value={values.three_gene}
            onChange={(e) => update("three_gene", e.target.value)}
            placeholder="ALK"
          />
        </label>
      </div>

      {/* Exon mode fields */}
      {mode === "exon" && (
        <div className="form-row">
          <label>
            5′ last exon
            <input
              required
              value={values.five_exon ?? ""}
              onChange={(e) => update("five_exon", e.target.value)}
              placeholder="13"
              inputMode="numeric"
            />
          </label>
          <label>
            3′ first exon
            <input
              required
              value={values.three_exon ?? ""}
              onChange={(e) => update("three_exon", e.target.value)}
              placeholder="20"
              inputMode="numeric"
            />
          </label>
        </div>
      )}

      {/* Genomic mode fields */}
      {mode === "genomic" && (
        <>
          <div className="form-row">
            <label>
              5′ genomic breakpoint
              <input
                required
                value={values.five_genomic ?? ""}
                onChange={(e) => update("five_genomic", e.target.value)}
                placeholder="chr2:42295516"
              />
            </label>
            <label>
              3′ genomic breakpoint
              <input
                required
                value={values.three_genomic ?? ""}
                onChange={(e) => update("three_genomic", e.target.value)}
                placeholder="chr2:29223528"
              />
            </label>
          </div>
          <div className="form-row">
            <label>
              Variant type <span className="hint">(optional)</span>
              <select
                value={values.variant_type ?? ""}
                onChange={(e) => update("variant_type", e.target.value || undefined)}
              >
                <option value="">— select —</option>
                {VARIANT_TYPES.map((vt) => (
                  <option key={vt.value} value={vt.value}>{vt.label}</option>
                ))}
              </select>
            </label>
          </div>
        </>
      )}

      <div className="form-row">
        <label>
          Genome build
          <select value={values.genome_build} onChange={(e) => update("genome_build", e.target.value)}>
            <option value="GRCh38">GRCh38</option>
            <option value="GRCh37">GRCh37</option>
          </select>
        </label>
      </div>

      <div className="form-row">
        <label>
          Tumor type <span className="hint">(optional)</span>
          <input
            value={values.tumor_type ?? ""}
            onChange={(e) => update("tumor_type", e.target.value || undefined)}
            placeholder="e.g. lung adenocarcinoma"
          />
        </label>
      </div>

      <div className="form-row">
        <label>
          5′ transcript <span className="hint">(optional)</span>
          <input
            value={values.five_transcript ?? ""}
            onChange={(e) => update("five_transcript", e.target.value)}
            placeholder="ENST00000318522 (default: canonical)"
          />
        </label>
        <label>
          3′ transcript <span className="hint">(optional)</span>
          <input
            value={values.three_transcript ?? ""}
            onChange={(e) => update("three_transcript", e.target.value)}
            placeholder="ENST00000389048 (default: canonical)"
          />
        </label>
      </div>

      <button type="submit" disabled={loading} className="submit-button">
        {loading ? "Annotating…" : "Annotate fusion"}
      </button>
    </form>
  );
}
