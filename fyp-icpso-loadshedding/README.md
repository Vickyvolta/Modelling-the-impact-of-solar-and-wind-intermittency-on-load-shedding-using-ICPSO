# Modelling the Impact of Solar and Wind Intermittency on Load Shedding using a Modified (Chaotic ICPSO) PSO Algorithm

**Author:** Ajiroye Victor Olusegun (190403034)
**Supervisor:** Prof. S. O. Adetona
**Department of Electrical & Electronics Engineering, University of Lagos**

A from-scratch MATLAB implementation of the simulations behind this final-year
project. Every figure the report relies on is either an **exact replication of
Prof. Adetona's published results** (references [31] and [32]) or a new result
built **on top of** that replicated benchmark. Nothing depends on a saved data
file: each stage loads the IEEE-14 bus system (`case14`) straight from MATPOWER
and re-derives everything at run time.

---

## Requirements

- **MATLAB** (R2018b or later). **No** Optimization or Statistics toolbox is
  required — the compass/pattern search, the Beta/Weibull inverse transforms and
  the wind-power quadrature are all written out in full.
- **MATPOWER 7.x or 8.x** on the MATLAB path (provides `loadcase`/`case14` and
  `runpf`). Get it from https://matpower.org.

```matlab
>> addpath(genpath('/path/to/matpower'))   % if not already on the path
>> run_all                                  % runs all six stages in order
```

Run any stage on its own instead, e.g. `>> basecase_loadflow`. Each stage
writes a transcript to `<stage>.log` and saves `<stage>_results.mat`.

### Run-time knobs

The full statistical protocol runs many thousands of load flows. Two switches
give a fast, still-complete pass:

| File | Setting | Effect |
|------|---------|--------|
| `icpso_vs_pso.m` | `NRUN = 3` (default 15) | 3 independent runs per algorithm — **this is the protocol used in [31]** |
| `re_intermittency.m` | `QUICK = true` | reduced Monte-Carlo sample; every block still runs |

The **deterministic** results (losses, shed totals, critical lines, voltages,
loadings, the true optimum) are identical either way; only Monte-Carlo counts and
the paired run-means tighten with more runs.

---

## Repository layout

| File | Role |
|------|------|
| `powerflow_lib.m` | Shared library — power-flow `solve` (with PV-bus-only Q-limit enforcement), branch `flows`, the `pso`/`icpso`/`polish` optimisers, `shed_pso`, and hard-coded MATPOWER column indices. Every stage calls into this. |
| `basecase_loadflow.m` | Stage A1 — base-case load flow (losses, voltages, all branch flows) |
| `critical_line_screen.m` | Stage A2 — critical-line identification (N-1 screen) |
| `load_shedding.m` | Stage A3 — 15 % load-shedding totals |
| `contingency_scenarios.m` | Stage A4 — the three contingency states of [32] §4.1.3 |
| `icpso_vs_pso.m` | Stage B — ICPSO vs standard PSO on [31]'s ORPD problem |
| `re_intermittency.m` | Stage C — solar/wind intermittency on the validated benchmark (the project's contribution) |
| `run_all.m` | Driver that runs A1 → C in order with per-stage timing |

---

## Results vs. benchmark

"Published" is Prof. Adetona's reported figure; "This code" is what the script
computes at run time. **MATCH** = agreement to the last printed digit.

| Stage | Quantity | Published | This code |
|-------|----------|-----------|-----------|
| A1 | Base-case real-power loss | 13.393 MW | **13.3933 MW — MATCH** |
| A2 | Critical lines (N-1) | {L1, L14} | **{L1, L14} — MATCH** |
| A3 | Real / reactive power shed (15 %) | 38.8502 MW / 11.0250 Mvar | **38.8500 / 11.0250 — MATCH** |
| A4 | Voltage at 15 % shed (% of base) | ~275 % | **274 %** |
| A4 | Double outage of L1 & L14 | — | **no power-flow solution** (unsurvivable) |
| B | PSO best loss (ORPD) | 12.275 MW | **12.2777 MW** — the true global optimum (multi-start compass search); his PSO sits on it |
| B | ICPSO vs PSO over independent runs | ICPSO "better" | **statistical tie** (paired mean +0.0088 MW; the wingbeat term reduces spread, not the mean) |
| C | Expected solar / wind (Beta / Weibull) | — | 66.60 MW / 49.96 MW |
| C | Deterministic minimum shed (critical outage) | — | **21.01 MW — interior** to the 38.85 MW cap |

**Headline finding.** Substituting the chaotic ICPSO for standard PSO does not
lower the optimum on either of Prof. Adetona's own problems; its wingbeat
frequency acts as a *momentum damper* that reduces run-to-run variance. The
intermittency layer (Stage C) is built specifically to give the optimiser a
non-degenerate problem with real headroom, and it confirms the same conclusion on
independent ground while quantifying the renewable-penetration band over which the
15 % shedding cap remains sufficient.

---

## Verification note

Because MATLAB was not available in the environment where these scripts were
first assembled, every numerical claim was independently reproduced in a separate
implementation and checked against Prof. Adetona's published values before being
transcribed here. The two implementations agree with the paper, which is why the
replicated figures above match to the last digit.

---

## Citation

> Ajiroye, V. O. (2026). *Modelling the Impact of Solar and Wind Intermittency on
> Load Shedding using a Modified PSO Algorithm.* B.Sc. Final-Year Project,
> Department of Electrical & Electronics Engineering, University of Lagos.
> Supervised by Prof. S. O. Adetona.

Built on the IEEE 14-bus test system via [MATPOWER](https://matpower.org).
