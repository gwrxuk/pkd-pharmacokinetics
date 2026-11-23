# PKD — Pharmacokinetics + Drug-Repurposing Research Tool (C++17)

> **Research / education only. Not clinical advice.** Outputs are computational hypotheses.

PKD reads drug monograph PDFs, extracts pharmacokinetic (PK) parameters and molecular
targets, simulates plasma concentration–time profiles, and ranks **new-use hypotheses
(drug repurposing)** for existing medicines against a disease–target knowledge base.
It does **not** design new chemical compounds.

```
PDF ──► text ──► monograph parser ──► DrugRecord ──► PK engine (RK4 ODE)
 (built-in       (units, ranges,       F, ka, CL, V,     Cmax, Tmax, AUC,
  FlateDecode     derivations,         targets+direction, steady state,
  or pdftotext)   target mining)       label indications  coverage
                                              │
                                              ▼
                            repurposing scorer vs data/diseases.csv
                     (direction-aware target match + pathway Jaccard
                       × PK exposure − contraindication penalty)
```

## Build & run

Requires CMake ≥ 3.16, a C++17 compiler and zlib. `pdftotext` (poppler) is optional
and is used as a fallback for complex PDFs (CID fonts).

```sh
cmake -S . -B build && cmake --build build -j
./build/pkd_tests                                   # 42 unit tests
./build/mkpdf samples/monographs samples/pdfs       # regenerate sample PDFs
./build/pkd discover samples/pdfs/sildenafil.pdf    # ingest a PDF, simulate, rank new uses
./build/pkd simulate Theophylline --doses 10 --csv theo.csv
./build/pkd rank Propranolol
./build/pkd validate                                # full validation report
```

## How it works

**PK model.** Linear compartment model `gut →ka→ central(V) ⇄Q⇄ peripheral(V2)`, elimination
CL, integrated with fixed-step RK4 across dose events. Checked against the analytic
Bateman solution (error < 1e-5), AUC = F·D/CL, and Css,avg = F·D/(CL·τ).

**Parser.** Finds `Label: value` fields, converts units (mL/min, L/min, L/kg ×70 kg,
ng/mL, mcg/mL, minutes, ranges → midpoint), and derives what's missing:
V = CL·t½/ln2, CL = V·ln2/t½, ka solved from Tmax, CL/F handled as apparent clearance.
Every derivation is reported as a note.

**Target mining.** Synonyms in `data/targets.csv` are matched inside *Mechanism of Action*
(longest match wins). Direction (inhibit/block → down, activate/agonist → up) comes from the
nearest action word in the same sentence.

**Score** for drug *d* and disease *s*:

```
target   = 1 − Π(1 − w)  over disease targets hit in the right direction
pathway  = Jaccard(pathways(d), pathways(s))
coverage = Css,avg / (Css,avg + Ceff)          (neutral if Ceff unknown)
conflict = 1 − Π(1 − w)  over targets hit in the WRONG direction
score    = (target + 0.3·pathway) · (0.5 + 0.5·coverage) − conflict
```

A negative score with a conflict is reported as a **safety flag**, for example propranolol
(β2 blockade) in asthma.

## Validation (from `pkd validate`; full output in `samples/validation_report.txt`)

| Check | Result |
|---|---|
| PK Cmax vs literature (12 drugs) | 12/12 within 2-fold, GMFE 1.28 |
| Indication recovery (14 drugs × 31 diseases) | AUROC 0.997, R@3 1.00 |
| Permutation null (1000 shuffles) | null AUROC 0.52, p = 0.001 |
| Known contraindications | 2/2 flagged |
| Retrospective PDF rediscovery | 8/9 later uses in top-3 (e.g. sildenafil → PAH, methotrexate → RA/psoriasis, propranolol → tremor/hemangioma) |

**Read these honestly:**
- The disease→target table was written with knowledge of how these drugs work, so the
  indication recovery and rediscovery numbers are **optimistic**. They check internal
  consistency, not predictive power. For a real test, swap in disease associations from an
  independent source (Open Targets, DisGeNET) and a time-split.
- Predicted Cmax is biased low (fold ≈ 0.6–0.7) because a one-compartment model spreads
  the dose over the full Vss immediately. Use V1/Q/V2 when you have them.
- Values in `data/` and `samples/` are **approximate** literature figures for teaching, and
  the sample PDFs are synthetic. Verify against primary sources before relying on them.
- Real FDA/EMA labels vary a lot in layout. The parser works best with `Label: value`
  style text, and every derived or assumed value is reported.

## Layout
```
include/pkd/   headers         src/        library + CLI (main.cpp)
data/          knowledge base  samples/    monograph .txt, generated .pdf, truth set, report
tools/mkpdf    text → PDF      tests/      unit tests
```
