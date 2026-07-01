# fusion-annotation

Standards-aligned **gene-fusion annotation**: a protein-effect engine (VEP-like) and
a knowledge engine (OncoKB-like), joined by an **HGVS.p-like protein-level interface**.

The core engine has **zero runtime dependencies** (Python stdlib only) and is fully
offline-testable. Live annotation sources (Ensembl, InterPro, CIViC, Open Targets)
are reached through a pluggable `DataProvider` — including an MCP-backed provider for
agentic and [Genome Nexus](https://www.genomenexus.org/) / cBioPortal workflows.

## Why

Fusion callers (STAR-Fusion, Arriba, FusionInspector, …) tell you *that* a fusion
exists. Oncogenicity scorers (OncoFuse, FusionPath, …) tell you *how likely it is to
be a driver*. Neither emits a normalized, HGVS-style protein-level description of the
**chimeric product**, nor a clean contract between the *effect* of a specific breakpoint
and the *knowledge* about a categorical gene pair. This package fills that gap:

| Layer | Analogy | Output |
|-------|---------|--------|
| 1 · Effect | VEP | chimeric protein: frame, junction, hybrid codon, retained/lost domains |
| 2 · Interface | HGVS.p / `::` | `EML4:p.Met1_Lys496::ALK:p.Tyr1059_Pro1620` + categorical key |
| 3 · Knowledge | OncoKB | oncogenicity, therapies, evidence for the gene-pair fusion |

The interface layer follows the [VICC Gene Fusion Specification](https://fusions.cancervariants.org/)
distinction between an **assayed** fusion (a specific breakpoint, annotated by Layer 1)
and a **categorical** fusion (a gene pair, keyed by Layer 3), and uses the HGVS `::`
adjoined-protein operator for the junction string.

## Install

```bash
pip install -e .            # core only, no dependencies
pip install -e ".[test]"    # + pytest
```

## Quickstart (offline)

```python
from fusion_annotation import Transcript, build_exon_cds_map, annotate_fusion
from fusion_annotation.providers import StaticProvider

# ... seed a StaticProvider (see examples/eml4_alk_offline.py) ...
result = annotate_fusion(provider, "EML4", "ALK", five_exon=13, three_exon=20)
print(result["interface"]["hgvsp_like"])
# EML4:p.Met1_Lys496::ALK:p.Tyr1059_Pro1620  (junction hybrid codon -> Val)
```

Run the bundled worked example:

```bash
python examples/eml4_alk_offline.py
```

## Live annotation via MCP

```python
# In a Claude Science repl cell, host.mcp is the upstream connector callable.
from fusion_annotation.providers import MCPDataProvider
from fusion_annotation import annotate_fusion

provider = MCPDataProvider(host.mcp)          # Ensembl + InterPro + CIViC
result = annotate_fusion(provider, "EML4", "ALK", five_exon=13, three_exon=20)
```

`MCPDataProvider` expects a callable with signature `mcp(server, method, **kwargs)` and
uses these upstream servers:

- **genomes** (Ensembl REST): `ensembl_lookup`, `ensembl_sequence`, `ensembl_xrefs`
- **protein-annotation** (InterPro): `get_domain_architecture`
- **clinical-genomics** (CIViC / Open Targets): `civic_search_molecular_profiles`, `civic_search_evidence`

## As an MCP tool

`fusion_annotation.mcp_tool` exposes a single `annotate_gene_fusion` tool
(`TOOL_SCHEMA` + `annotate_fusion_tool()` backend) that you can register with any MCP
server framework, so an LLM agent can annotate a fusion in one call. See the module
docstring for a FastMCP example.

## The EML4::ALK worked example

`EML4::ALK` variant 1 (E13;A20) — EML4 exon 13 joined to ALK exon 20 — is the canonical
NSCLC driver. This package reproduces, from primary data:

- fusion CDS = 3180 nt → **in-frame**, **zero internal stops**, 1059 aa protein
- junction: EML4 contributes 496 complete codons (…Lys496); the junction codon is a
  **hybrid** (1 nt EML4 + 2 nt ALK = GTG = **Val**); ALK continues from Tyr1059
- **ALK kinase domain (1116–1392) fully retained**; EML4 second β-propeller lost →
  the mechanism: loss of ALK's extracellular/TM region + EML4-driven oligomerization
  of an intact kinase = constitutive activation

See [docs/DESIGN.md](docs/DESIGN.md) for the full architecture and rationale.

![EML4::ALK chimeric protein and domain retention](docs/fusion_domain_map.png)

## Tests

```bash
pytest            # 11 offline assertions against the EML4::ALK truth values
```

## Status & roadmap

v0.1 handles exon-boundary breakpoints on canonical transcripts. Planned: genomic
(non-exon-boundary) breakpoints, all-transcript enumeration, NMD prediction,
full VICC GFS JSON schema + GA4GH VRS/Cat-VRS identifiers, `hgvs` library round-trip
validation, and a proper OncoKB backend. See docs/DESIGN.md §6.

## License

Apache-2.0 (see [LICENSE](LICENSE))
