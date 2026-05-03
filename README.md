# Optimizing an Investment Portfolio under Solvency II

Bachelor thesis (Econometrics & Operations Research, Maastricht University) that maximizes the expected one-period surplus of a life insurer's investment portfolio, subject to the risk constraints imposed by the Solvency II standard formula. The problem is formulated as a linear program and solved in R.

**Author:** Florian Koenigstorfer

**Supervisor:** Prof. Dr. Antoon Pelsser

**Submitted:** March 22, 2017

## Problem & solution

Solvency II imposes pre-defined shocks (interest rate, equity, mortality, etc.) on an insurer's balance sheet, and the insurer must remain solvent after each shock. This thesis turns those rules into LP constraints and maximizes the sample-average expected surplus at `T = 1`.

The objective function uses the sample counterpart of expected surplus computed across `N` simulated economic scenarios. Each scenario contributes one constraint per shock, plus a budget constraint and non-negativity constraints (no short-selling, no credit lines).

## Highlights

- LP formulation of the surplus-maximization problem with Solvency II shocks as constraints.
- Toy model (3 bonds, 3 scenarios, 2 interest-rate shocks) used to illustrate the setup and check binding constraints.
- Stability analysis showing that increasing the number of simulated economic scenarios (`N`) from 3 to 200 stabilizes the optimal allocation.
- Mortality shock scenario showing that a budget increase from 1000 to 3000 restores feasibility.
- Mixed bond/equity portfolio analysis showing how the bond/stock split shifts as the budget grows.

The script runs end-to-end and prints/plots the results referenced in the thesis (toy model setup, optimal allocations, binding constraints, mortality-shock results, and the bond/stock-vs-budget plots).

## Usage

The full analysis lives in `CodeThesis.R`. Sourcing the file runs the entire pipeline; the script is organized into 14 numbered sections that mirror the thesis. Sections 1–6 are reusable building blocks (define them once); sections 7–14 are the four model variants from §4–§5 of the thesis plus the summary printout.

A browser-viewable version of the full pipeline is also provided as `CodeThesis.html` output, so the results can be inspected without re-running R.

### Building blocks (sections 1–6)

- **§1 Term-structure & cashflow primitives** — `draw_interest_term_structure`, `shock_term_structure`, `build_bond_cashflows`, `discounted_sum_bond`, `discounted_sum_liabilities`.
- **§2 Investment universe** — `build_investment_universe`, `calc_initial_values`.
- **§3 Scenario engine** — `shock_liability_value`, `shock_bond_value`, `shock_equity_value` (per-shock dispatch helpers), and `generate_scenarios` (the main scenario loop, replaces the original `NStatesJShocks`).
- **§4 LP construction** — `calc_expected_surplus`, `build_constraint_inputs`, `build_lp_model`.
- **§5 Diagnostics** — `find_binding_constraints`, `check_allocation_feasibility`, `report_solver_status`.
- **§6 End-to-end driver** — `solve_portfolio` (runs scenarios → objective → constraints → LP → solve in one call) and `summarise_solution` (extracts allocation, expected return, % return).

### Analysis blocks (sections 7–14)

- **§7  Thesis §4.1–4.2 — Toy model.** `toy_model`, `toy_summary$allocation`, `toy_summary$expected_return`, `toy_summary$percentage_return`.
- **§8  Thesis §4.3 — Binding constraints.** `binding_toy` (1 = binding, 0 = slack; final `K + M` entries are the non-negativity bounds).
- **§9  Thesis §4.4 — Old allocation under freshly-drawn constraints.** `toy_model_redraw`, `toy_redraw_summary`, `feasibility_old_in_new` (indicator vector of failing constraints), plus `return_new_alloc_old_obj` / `return_old_alloc_new_obj` cross-evaluations.
- **§10 Thesis §5.1 — Stability with `N = 200`.** `big_model_a`, `big_model_b`, `big_summary_a`, `big_summary_b`.
- **§11 Thesis §5.2 — Mortality shock.** `mort_model_low_budget$status` (= 2, infeasible at budget 1000), then `mort_status_hi`, `mort_alloc_hi`, `mort_return_hi`, `mort_pct_hi` at budget 3000.
- **§12 Thesis §5.3 — Equities + equity shock.** `eq_model`, `eq_summary` (3 bonds + 2 equity indices, mortality shock excluded).
- **§13 Thesis §5.4 — Budget sensitivity.** Sweeps the budget from 800 upward in 2 % increments and records `allocation_grid`, `budget_seq`, `bond_share_seq`, `equity_share_seq`.
- **§14 Summary printout** — prints results for each section and renders the two budget-vs-share charts via `qplot`.

To inspect a specific result interactively, open `CodeThesis.R` in RStudio, run the file from the top through the end of the §6 driver, then run the section of interest. The `cfg_toy` / `cfg_big` / `cfg_mort` / `cfg_eq` config lists at the start of each analysis block make it easy to tweak parameters without touching the building blocks.

## Project structure

```
.
├── CodeThesis.R         # Full analysis script (building blocks + all thesis sections)
├── BScThesisEOR.pdf     # Thesis writeup
├── install_packages.R   # Pinned R dependency installer
└── README.md            # This file
```

## Caveats

- **No external data.** Interest term structures are drawn from a uniform distribution and equity returns from a log-normal distribution; the script is fully self-contained. A commented-out `# set.seed(1)` line sits at the top of §0 — uncomment it if you need bit-for-bit reproducibility. (The Quarto wrapper sets the seed explicitly so its rendered output is stable.)
- **`lpSolve::lp()` builds and solves in a single call**, returning a regular R object whose `$solution` and `$objval` fields hold the optimum. The `solve_portfolio` driver returns this object in its `$lp` field — read results via `model$lp$solution` and `model$lp$objval`. The objective is negated inside `build_lp_model` (the solver is called with `direction = "min"` on `-obj * x`) so that the maximisation problem is solved as a minimisation; downstream code such as `summarise_solution` un-negates `$objval` before reporting.
- **`qplot()`** is used at the end of §14 for the budget-sensitivity plots. It still works in modern ggplot2 but is soft-deprecated; the pinned `ggplot2 3.5.1` keeps the original output clean.
- **Linux audio:** `beepr::beep()` calls inside `report_solver_status` will silently fail (with a warning) if no audio backend is on PATH. The analysis runs fine without sound.
- **`install_packages.R`** uses `remotes::install_version()` to pin `lpSolve`, `beepr`, and `ggplot2` to the versions used to produce the rendered HTML. Update the pins only after re-running and re-rendering the full pipeline.

## Errata
**Mortality shock direction is inconsistent with the liability formula.** - During a later model review, I identified a mismatch between the thesis prose and the liability cashflow specification. The implemented liability model is appropriate for an annuity/pension-style book, where higher mortality reduces future payments; the thesis text described a death-benefit interpretation, where higher mortality would increase claims. I preserve the original thesis result for reproducibility and document the corrected interpretation here.
discounted_sum_liabilities models the liability cashflow at year t as L0 * exp(-lambda * t) — a survivorship curve, where lambda plays the role of force of mortality and the insurer pays surviving policyholders. Under this model the PV of liabilities is monotonically decreasing in lambda: higher mortality means fewer survivors and smaller liabilities. This is the right shape for an annuity / pension book.
The thesis prose in §5.2, however, describes the shock as a death-benefit story: "an increase in the mortality rate of the insured. The effect of this increase in the mortality is an increase in the number of claims to the insurer." That framing is incompatible with the survivorship formula actually used in the code.
The §5.2 result still holds — the model is infeasible at budget 1000 and feasible at budget 3000 — but the mechanism the thesis gives for it ("more claims tighten the constraint") doesn't match what the code computes. The infeasibility is a real LP-corner effect rather than a constraint-tightening effect.
The cleanest fix would be to rebrand the §5.2 shock as a longevity shock (set mort_shock_size = +0.15 so lambda_new = 0.85 * lambda; people live longer; PV of liabilities rises; constraint genuinely tightens). Solvency II's Life Risk module has a dedicated longevity sub-module precisely for annuity books, so this framing is well-motivated. The repository has not been modified to apply this fix; this note records the issue and the fix so a reader of the code is not misled.