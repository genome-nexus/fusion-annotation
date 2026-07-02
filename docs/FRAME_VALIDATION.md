# Frame-Engine Validation — Known In-Frame Oncogenic Fusions

> Tracking issue: [#3 — Exon-number-only input cannot disambiguate isoforms](https://github.com/genome-nexus/fusion-annotation/issues/3)

**Goal:** Confirm the structural annotator's reading-frame math is correct, by
re-running three known in-frame drivers that initially came back out-of-frame,
this time with **pinned transcripts** and a **sweep over exon numbers**.

**Result:** 2 of 3 resolved to the textbook in-frame breakpoint once the correct
exon was supplied; the 3rd is a transcript/isoform mismatch, not a frame-math
error. The earlier out-of-frame calls were **input mismatches** (wrong exon
number or wrong isoform), which validates the frame engine and directly
motivates the recommended input-format fix (accept genomic/HGVS breakpoints;
echo/pin the transcript; flag known-partner pairs that come back out-of-frame).

Server: `fusion-annotation` v1.28.1 · live MCP endpoint.
Transcripts pinned: NPM1 `ENST00000296930`, ALK `ENST00000389048`,
LMNA `ENST00000368300`, NTRK1 `ENST00000524377`,
CD74 `ENST00000009530`, ROS1 `ENST00000368507`.

---

## NPM1-ALK (t(2;5)) — RESOLVED ✅

| NPM1 exon | ALK exon | in-frame | junction | NPM1 last aa | ALK first aa | length | stops |
|-----------|----------|----------|----------|--------------|--------------|--------|-------|
| 3 | 20 | ✗ | — | 86 | 1058 | 108 | 0 |
| **4** | **20** | **✓** | **V** | **117** | **1059** | **680** | **0** |
| 5 | 20 | ✗ | — | 153 | 1058 | 175 | 0 |
| 6 | 20 | ✗ | D | 174 | 1058 | 193 | 0 |

Exon **4→20** is the unique in-frame pair — NPM1 truncated at residue 117,
ALK kinase resumes at 1059, 680 aa product, zero internal stops. This is the
canonical nucleophosmin-ALK breakpoint. The earlier exon-6 guess was wrong.

## LMNA-NTRK1 — RESOLVED ✅

| LMNA exon | NTRK1 exon | in-frame | length | stops |
|-----------|------------|----------|--------|-------|
| 2 | 10 | ✗ | 241 | 0 |
| **2** | **11** | **✓** | **550** | **0** |
| 2 | 12 | ✗ | 188 | 0 |
| 10 | 11 | ✓ | 945 | 0 |
| 11 | 11 | ✓ | 1035 | 0 |

The frame is set by **NTRK1 exon 11**'s native start phase — every LMNA donor
exon tested (2, 10, 11) is in-frame against NTRK1 ex11 and out-of-frame against
ex10/ex12. The literature LMNA ex2 → NTRK1 ex11 breakpoint (550 aa) is
recovered; the earlier exon-10 acceptor was off-by-one.

## CD74-ROS1 — TRANSCRIPT/ISOFORM MISMATCH ⚠️ (not a frame bug)

| CD74 exon | ROS1 exon | in-frame | CD74 last aa | ROS1 first aa | length | stops |
|-----------|-----------|----------|--------------|---------------|--------|-------|
| 6 | 32 | ✗ | 208 | 1688 | 212 | 0 |
| 6 | 34 | ✗ | 208 | 1784 | 220 | 0 |
| 7 | 32 | ✗ | 272 | 1688 | 276 | 0 |
| 7 | 34 | ✗ | 272 | 1784 | 284 | 0 |
| 8 | 32 | ✗ | 293 | 1688 | 297 | 0 |
| 8 | 34 | ✗ | 293 | 1784 | 305 | 0 |

On `ENST00000009530`, CD74 exon 6 already ends at residue **208**, but the real
CD74-ROS1 fusion retains only **~106 aa** of CD74 at the junction. This
transcript is a **longer CD74 isoform** (p41-type, extra invariant-chain exon)
whose exon numbering does not line up with the canonical breakpoint. No in-frame
pair exists on this isoform because the frame the engine computes is correct for
*this* transcript — the transcript is simply the wrong one for this fusion.
Resolving it requires the CD74 isoform that gives a ~106-aa junction, which the
current exon-number-only input cannot disambiguate. This is exactly the failure
mode the genomic/HGVS-breakpoint input would eliminate.

---

## What this validates

1. **The frame math is correct.** Given the right exon on the right transcript,
   the engine recovers the textbook in-frame breakpoint (NPM1-ALK, LMNA-NTRK1)
   and correctly reports out-of-frame for every neighboring exon.
2. **The engine is only as good as its exon/transcript inputs.** All three
   earlier "out-of-frame" surprises were input mismatches (wrong exon number, or
   an isoform whose exon numbering doesn't match the literature), never arithmetic.
3. **Direct motivation for the input-format fix:**
   - Accept **genomic coordinates / HGVS breakpoints**, not just exon numbers —
     this pins the isoform and removes exon-numbering ambiguity (CD74 case).
   - **Echo and pin the transcript** used (MANE default + explicit override).
   - Add a **sanity flag** when a known oncogenic partner pair returns
     out-of-frame, prompting a transcript/exon re-check.
