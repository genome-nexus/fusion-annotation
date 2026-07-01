# Gene Fusion Annotation System — Design & Reference Prototype

A two-engine architecture for annotating gene fusions, joined by an HGVS.p-like
protein-level interface:

1. **Effect engine (VEP-like)** — reconstructs the chimeric protein from an
   *assayed* fusion (two transcripts + breakpoints) and predicts its molecular
   consequence: reading frame, junction residue, domain retention.
2. **Knowledge engine (OncoKB-like)** — attaches curated clinical knowledge
   (oncogenicity, therapies, evidence) to a *categorical* fusion (a gene pair).
3. **Interface (HGVS.p-like)** — a normalized, computable protein-level
   representation that the effect engine emits and the knowledge engine consumes.

This document specifies the architecture and demonstrates it end-to-end on the
canonical **EML4::ALK variant 1** fusion, using live data from Ensembl,
InterPro, CIViC, and Open Targets.

---

## 1. Why this shape — and how it maps onto existing standards

The three-part vision maps cleanly onto standards that consolidated in
2021–2024, which is what makes it buildable rather than bespoke:

| Your layer | Standard it mirrors | Key concept |
|---|---|---|
| VEP-like effect | Ensembl VEP; HGVS protein consequences | chimeric protein = a delins at the junction |
| HGVS.p-like interface | HGVS `::` operator; VICC Gene Fusion Spec information model | adjoined transcripts/proteins; assayed vs. categorical |
| OncoKB-like knowledge | OncoKB / CIViC / Cat-VRS | categorical fusion concept ("ALK Fusions") |

Three load-bearing distinctions from the **VICC Gene Fusion Specification**
(fusions.cancervariants.org) and **HGVS nomenclature** drive the design:

- **Assayed vs. categorical fusion.** An *assayed* fusion is one concrete
  observation with defined breakpoints (`EML4 exon 13 :: ALK exon 20`). A
  *categorical* fusion is the abstract class a knowledgebase curates against
  (`EML4::ALK`, or even more broadly `ALK Fusions`). **The effect engine works
  on assayed fusions; the knowledge engine keys on categorical fusions.** The
  interface object carries both, and its `.categorical_key()` is the bridge.
- **The `::` junction operator.** HGVS and HGNC both adopt `GENE1::GENE2` with a
  double-colon for fusions; at the protein level the chimeric junction is
  written with `::` joining the two contributing segments. This is exactly the
  "hgvs.p-like" string you envisioned — it is a real, endorsed notation, not an
  invention.
- **Chimeric protein as delins.** HGVS treats a fusion protein as a special
  deletion-insertion: the C-terminus of the 5′ partner is "replaced" by the
  3′ partner's C-terminal portion. The junction may fall mid-codon, producing a
  **hybrid codon** whose residue belongs to neither parent — a detail any
  correct effect engine must handle.

---

## 2. Architecture

```
  ASSAYED FUSION                 INTERFACE (hgvs.p-like)            CATEGORICAL
  (transcripts + breakpoints)    normalized protein object          KNOWLEDGE
        │                              │                                │
        ▼                              ▼                                ▼
  ┌───────────────┐   emits    ┌─────────────────┐   .categorical_  ┌──────────────┐
  │ EFFECT ENGINE │──────────▶ │ FusionProtein   │─────key()──────▶ │ KNOWLEDGE    │
  │  (VEP-like)   │            │  ::-junction    │                  │ ENGINE       │
  │               │            │  frame status   │                  │ (OncoKB-like)│
  │ Ensembl CDS   │            │  domain calls   │                  │ CIViC /      │
  │ + InterPro    │            │  to_hgvsp()     │                  │ Open Targets │
  └───────────────┘            └─────────────────┘                  └──────────────┘
```

Data flows through a pluggable `DataProvider` protocol so the pure-Python core
is offline-testable and the live sources are swappable:

- **Effect engine** → Ensembl REST (`ensembl_lookup` with `expand=True` for exon
  structure, `ensembl_sequence` for CDS/protein) + InterPro
  (`get_domain_architecture`) for domains.
- **Knowledge engine** → CIViC (`civic_search_molecular_profiles` →
  `civic_search_evidence`) as the open-source OncoKB stand-in, plus Open Targets
  for target–drug tractability. (No public OncoKB connector; CIViC is the
  natural open equivalent and shares the categorical-fusion model.)

---

## 3. The interface object (the HGVS.p-like contract)

`FusionProtein` is the normalized representation. Its serializers are the API:

- **`to_hgvsp()`** → the junction string, e.g.
  `EML4:p.Met1_Lys496::ALK:p.Tyr1059_Pro1620 (junction hybrid codon → Val)`
- **`categorical_key()`** → `EML4::ALK`, the knowledge lookup key.
- **`to_dict()`** → full computable record (frame status, domain calls, sequence).

Fields the interface must carry for the two engines to interoperate:

| Field | Meaning | Consumed by |
|---|---|---|
| `five_last_aa`, `five_last_aa_res` | last residue fully from 5′ partner | junction string |
| `hybrid_codon`, `junction_residue` | mid-codon junction → chimeric residue | effect / display |
| `three_first_aa`, `three_first_aa_res` | first residue fully from 3′ partner | junction string |
| `in_frame`, `frame_status` | in-frame / out-of-frame / frameshift-truncating | oncogenicity heuristic |
| `internal_stops` | premature stop count | NMD/truncation flag |
| `domains[]` (RETAINED/LOST/DISRUPTED) | domain-level consequence | mechanism, oncogenicity |
| `categorical_key()` | gene-pair class | knowledge engine |

**Frame status is the single most predictive computed feature.** An in-frame
fusion that retains a kinase domain and a dimerization module is the classic
oncogenic driver pattern; an out-of-frame or truncating junction usually is not.

---

## 4. Worked example — EML4::ALK variant 1 (E13;A20)

All values below are **computed from primary data**, not recalled. The effect
engine reconstructs the chimeric CDS from the real Ensembl coding sequences and
the exon→CDS coordinate map; results were verified by reconstructing each parent
segment exactly.

**Inputs (Ensembl GRCh38, canonical transcripts):**

| Partner | Gene | Transcript | Strand | Protein | UniProt |
|---|---|---|---|---|---|
| 5′ | EML4 (ENSG00000143924) | ENST00000318522 | + | 981 aa | Q9HC35 |
| 3′ | ALK (ENSG00000171094) | ENST00000389048 | − | 1620 aa | Q9UM73 |

Both genes are on chr2; the fusion arises from the well-described inv(2) that
juxtaposes the + strand EML4 and − strand ALK.

**Effect engine output:**

- Breakpoint: EML4 exon 13 (CDS ends at nt 1489) joined to ALK exon 20 (CDS
  starts at nt 3173).
- Reconstructed fusion CDS = **3180 nt**, divisible by 3 → **in-frame**, with
  **zero internal stop codons**.
- Fusion protein = **1059 aa**.
- **Junction:** EML4 contributes 496 complete codons (…Lys496); the next codon
  is a **hybrid** — 1 nt from EML4 (G) + 2 nt from ALK (TG) = GTG → **Val497**;
  ALK then continues from its first fully-retained residue **Tyr1059** to
  Pro1620.

**HGVS.p-like interface string:**

```
EML4:p.Met1_Lys496::ALK:p.Tyr1059_Pro1620   (junction hybrid codon → Val)
categorical key: EML4::ALK
```

**Domain retention (InterPro):**

| Domain | Coords | Call |
|---|---|---|
| EML4 HELP motif | 255–293 | RETAINED |
| EML4 first β-propeller (partial) | 301–496 | DISRUPTED at breakpoint |
| EML4 second β-propeller + WD40 repeats | 498–864 | LOST |
| ALK extracellular (MAM / LDLa / Gly-rich) | 264–961 | LOST |
| **ALK protein kinase domain** | **1116–1392** | **RETAINED** |
| ALK ATP-binding site | 1122–1150 | RETAINED |
| ALK tyrosine-kinase active site | 1245–1257 | RETAINED |

**Mechanistic read (computed, not asserted):** the fusion **loses the entire
ALK extracellular/transmembrane region** and **retains the full ALK tyrosine
kinase domain**, now placed under the control of EML4's retained
N-terminal/coiled-coil portion. That partner-driven oligomerization is what
constitutively activates the orphaned kinase — the canonical driver mechanism,
here derived directly from frame + domain calls rather than looked up.

**Knowledge engine output (CIViC molecular profile `EML4::ALK Fusion`, id 5):**

- 10 predictive evidence items (8 sensitivity/response, 2 resistance).
- Diseases: lung non-small cell carcinoma (5), lung adenocarcinoma, mesothelioma,
  colorectal, high-grade glioma.
- Therapies cited: crizotinib (×4), alectinib, lorlatinib, entrectinib, others.
- Open Targets: 12 ALK-targeting drugs; **7 approved ALK inhibitors** — lorlatinib,
  alectinib, brigatinib, crizotinib, entrectinib, ceritinib.

The knowledge attaches at the **categorical** level (`EML4::ALK`, indeed any
`ALK` fusion for most of these drugs), which is exactly why the interface's
`categorical_key()` — not the exact breakpoint — is the join key.

---

## 5. Module

`fusion_annotator.py` implements this design in pure Python (no MCP dependency,
so the core is unit-testable offline):

- `Transcript`, `build_exon_cds_map`, `cds_coord_at_exon_boundary` — Layer-1 inputs.
- `annotate_effect(...)` → `FusionProtein` — the VEP-like engine.
- `FusionProtein.to_hgvsp()` / `.categorical_key()` / `.to_dict()` — the interface.
- `annotate_knowledge(fp, provider)` → `FusionKnowledge` — the OncoKB-like engine.
- `annotate_fusion(provider, ...)` — end-to-end orchestration.
- `DataProvider` protocol — implement `get_transcript` / `get_domains` /
  `get_fusion_knowledge` against the MCP connectors (from the `repl` tool) or any
  other backend.

The module was validated by reproducing every computed value in §4 (fusion
length 1059, zero internal stops, hybrid Val497 junction, ALK kinase RETAINED),
with each parent segment reconstructed exactly from the fetched sequences.

---

## 6. Roadmap — extending toward production

**Effect engine**
- Handle intronic/mid-exon genomic breakpoints (not just exon boundaries) by
  mapping genomic coordinates → CDS via the exon table already built.
- Non-canonical/all-transcript enumeration; report per-transcript frame status.
- NMD prediction for out-of-frame/PTC-bearing junctions (55-nt-from-last-EJC rule).
- 5′-partner retention of dimerization motifs (coiled-coil detection) as an
  explicit oncogenicity feature, not just a domain-retention side effect.

**Interface / nomenclature**
- Full alignment to the VICC Gene Fusion Specification JSON schema (categorical
  vs. assayed objects, regulatory-element and multi-partner cases) for
  interoperable exchange.
- Emit GA4GH VRS / Cat-VRS identifiers so records are globally addressable.
- Round-trip validation against `hgvs`/`cool-seq-tool`-style libraries.

**Knowledge engine**
- Add OncoKB proper when a connector/API key is available; keep CIViC + Open
  Targets as the open fallback. Normalize evidence levels across sources
  (OncoKB levels ↔ CIViC A–E ↔ AMP/ASCO/CAP tiers).
- Reciprocal/partner-agnostic lookup (e.g. any `ALK` fusion vs. `EML4::ALK`
  specifically), matching how knowledgebases scope their assertions.
- Fold in ClinGen/ClinVar and trial matching (ClinicalTrials.gov connector).

**Delivery**
- Wrap `annotate_fusion` as an MCP tool so it plugs into your agentic cBioPortal
  workflows directly — input a fusion call, get back the interface object +
  knowledge in one step.
