import type { AnnotationResult } from "../lib/types";
import { DomainDiagram } from "./DomainDiagram";

interface Props {
  result: AnnotationResult;
  permalink: string;
}

export function ResultView({ result, permalink }: Props) {
  const { interface: iface, knowledge, resolved, warnings } = result;

  return (
    <div className="result-view">
      <div className="permalink-row">
        <input readOnly value={permalink} onFocus={(e) => e.currentTarget.select()} />
        <button
          type="button"
          onClick={() => navigator.clipboard?.writeText(permalink)}
        >
          Copy permalink
        </button>
      </div>

      {warnings.length > 0 && (
        <ul className="warnings">
          {warnings.map((w) => (
            <li key={w}>⚠ {w}</li>
          ))}
        </ul>
      )}

      <h2 className="hgvsp">{iface.hgvsp_like}</h2>

      <dl className="summary-grid">
        <dt>Categorical fusion</dt>
        <dd>{iface.categorical_key}</dd>
        <dt>Frame</dt>
        <dd className={iface.in_frame ? "ok" : "bad"}>
          {iface.frame_status} (protein {iface.fusion_length} aa, internal stops {iface.internal_stops})
        </dd>
        <dt>Junction</dt>
        <dd>
          {iface.five_gene} {iface.five_last_aa_res}
          {iface.five_last_aa} :: {iface.three_gene} {iface.three_first_aa_res}
          {iface.three_first_aa}
          {iface.hybrid_codon && iface.junction_residue
            ? ` (hybrid codon → ${iface.junction_residue})`
            : ""}
        </dd>
        <dt>Genome build</dt>
        <dd>{resolved.genome_build}</dd>
        <dt>Transcripts used</dt>
        <dd>
          {resolved.five.transcript} ({resolved.five.transcript_source}) / {resolved.three.transcript} (
          {resolved.three.transcript_source})
        </dd>
      </dl>

      <h3>Domain retention</h3>
      <DomainDiagram
        domains={iface.domains}
        fiveGene={iface.five_gene}
        threeGene={iface.three_gene}
        fiveLastAa={iface.five_last_aa}
        threeFirstAa={iface.three_first_aa}
        fiveLength={resolved.five.protein_length}
        threeLength={resolved.three.protein_length}
      />

      <table className="domain-table">
        <thead>
          <tr>
            <th>Status</th>
            <th>Gene</th>
            <th>Domain</th>
            <th>Range</th>
          </tr>
        </thead>
        <tbody>
          {iface.domains.map((d) => (
            <tr key={`${d.gene}-${d.accession}-${d.start}-${d.end}`} className={`status-${d.status.toLowerCase()}`}>
              <td>{d.status}</td>
              <td>{d.gene}</td>
              <td>{d.name}</td>
              <td>
                {d.start}-{d.end}
              </td>
            </tr>
          ))}
        </tbody>
      </table>

      <h3>Clinical knowledge</h3>
      <dl className="summary-grid">
        <dt>Oncogenic</dt>
        <dd>{knowledge.oncogenic ?? "—"}</dd>
        <dt>Therapies</dt>
        <dd>{knowledge.therapies.length ? knowledge.therapies.join(", ") : "—"}</dd>
        <dt>Diseases</dt>
        <dd>{knowledge.diseases.length ? knowledge.diseases.join(", ") : "—"}</dd>
        <dt>Sources</dt>
        <dd>{knowledge.sources.length ? knowledge.sources.join(", ") : "—"}</dd>
      </dl>

      <p className="disclaimer">
        This is a research/informatics tool, not a diagnostic device. Results should be reviewed by a
        qualified professional before they inform patient care.
      </p>
    </div>
  );
}
