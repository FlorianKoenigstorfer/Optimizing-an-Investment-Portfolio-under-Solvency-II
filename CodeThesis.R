###############################################################################
# Optimising an Investment Portfolio under Solvency II
#
# Companion R script for the bachelor thesis (1.7) "Optimising an Investment
# Portfolio under Solvency II" by Florian Königstorfer (I6080626).
#
# Overview
# --------
# The model maximises an insurer's expected surplus at T = 1, subject to
#   (a) a Solvency-II–style risk constraint that requires asset value to
#       dominate liability value under each of J stress scenarios in each of
#       N simulated economic states, and
#   (b) a budget constraint on the cost of the portfolio at T = 0.
# Decision variables are the units held of K bonds and M Type-1 equities.
#
# Pipeline
#   1. Build the investment universe (random bond cashflows + equity slots).
#   2. Simulate N economic states; in each state apply up to four shocks:
#        1 = interest rate UP        2 = interest rate DOWN
#        3 = mortality decrease      4 = equity decrease (-39%)
#   3. Compute expected surplus (objective) and the LP constraints.
#   4. Solve via lpSolve; report binding constraints and feasibility.
#   5. Repeat the pipeline for the four model variants in §4–§5 of the thesis.
#
# Setup
#   Run install_packages.R once to obtain the pinned package versions, then
#     source("CodeThesis.R")
###############################################################################


# ============================================================================
# 0. Setup -------------------------------------------------------------------
# ============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(lpSolve)      # LP construction and solver (lp() builds + solves)
  library(beepr)        # Audible solver-status notifications
  library(ggplot2)      # qplot() for the budget-sensitivity charts in §5.4
})

# To exactly reproduce a single run, uncomment:
set.seed(1)

# Shock identifiers used throughout the scenario engine.
SHOCK_INTEREST_UP   <- 1L
SHOCK_INTEREST_DOWN <- 2L
SHOCK_MORTALITY     <- 3L
SHOCK_EQUITY        <- 4L


# ============================================================================
# 1. Term-structure & cashflow primitives ------------------------------------
# ============================================================================

#' Draw a linear risk-free term structure
#'
#' Endpoints `a` and `b` are drawn independently from U(lower, upper); rates
#' are linearly interpolated between them across `horizon` years. The returned
#' vector has `horizon + 1` elements (years 0..horizon), matching the
#' original `seq(a, b, by = (b - a) / horizon)` behaviour.
draw_interest_term_structure <- function(lower, upper, horizon) {
  a <- runif(1, lower, upper)
  b <- runif(1, lower, upper)
  step <- (b - a) / horizon
  seq(from = a, to = b, by = step)
}

#' Apply a parallel shift to a term structure
shock_term_structure <- function(term_structure, magnitude,
                                 direction = c("up", "down")) {
  direction <- match.arg(direction)
  shift <- if (direction == "up") magnitude else -magnitude
  term_structure + shift
}

#' Build a level-coupon bond cashflow vector
build_bond_cashflows <- function(maturity, coupon, face_value) {
  cf <- rep(coupon, times = maturity)
  cf[maturity] <- cf[maturity] + face_value
  cf
}

#' Discounted-cashflow value of a bond
#'
#' Mirrors the thesis' `DiscountedSummation`: when valuing at T = 1, the
#' first cashflow is rolled forward one period (CF1 * (1 + i1)) and
#' subsequent cashflows are discounted from year `year` to T = 0.
discounted_sum_bond <- function(term_structure, cashflows, valuation_time) {
  start <- 1L + valuation_time
  total <- 0
  if (valuation_time == 1L) {
    total <- cashflows[1] * (1 + term_structure[1])
  }
  if (start <= length(cashflows)) {
    for (year in start:length(cashflows)) {
      i  <- as.numeric(term_structure[year])
      cf <- as.numeric(cashflows[year])
      total <- total + cf / ((1 + i) ^ year)
    }
  }
  total
}

#' Discounted-cashflow value of liabilities under exponential mortality
#'
#' Liability cashflow at year t is L0 * exp(-lambda * t); discounted using
#' the supplied term structure.
discounted_sum_liabilities <- function(term_structure, L0,
                                       valuation_time, lambda) {
  start <- 1L + valuation_time
  total <- 0
  if (start <= length(term_structure)) {
    for (year in start:length(term_structure)) {
      i <- as.numeric(term_structure[year])
      total <- total + L0 * exp(-lambda * year) / ((1 + i) ^ year)
    }
  }
  total
}


# ============================================================================
# 2. Investment universe -----------------------------------------------------
# ============================================================================

#' Build the (K + M) x horizon cashflow matrix
#'
#' Rows 1..K hold randomly-generated bond cashflows; rows (K+1)..(K+M) hold
#' the equity investment sizes in column 1 (treated as zero-coupon — equity
#' has no contractual future cashflows).
build_investment_universe <- function(K, mat_lower, mat_upper,
                                      coupon_lower, coupon_upper,
                                      fv_lower, fv_upper,
                                      M = 0, equity_sizes = numeric(0)) {
  Y <- matrix(0, nrow = K + M, ncol = mat_upper)

  if (K > 0) {
    for (k in seq_len(K)) {
      maturity <- round(runif(1, mat_lower, mat_upper))
      coupon   <- runif(1, coupon_lower, coupon_upper)
      fv       <- runif(1, fv_lower, fv_upper)
      cf       <- build_bond_cashflows(maturity, coupon, fv)
      Y[k, seq_len(maturity)] <- cf
    }
  }
  if (M > 0) {
    Y[(K + 1):(K + M), 1] <- equity_sizes
  }
  Y
}

#' T = 0 market values of the (K + M) assets
calc_initial_values <- function(K, Y, term_structure,
                                M = 0, equity_sizes = NULL) {
  if (is.null(equity_sizes)) equity_sizes <- rep(1, M)
  Y0 <- numeric(K + M)
  if (K > 0) {
    for (k in seq_len(K)) {
      Y0[k] <- discounted_sum_bond(term_structure, Y[k, ], valuation_time = 0)
    }
  }
  if (M > 0) {
    Y0[(K + 1):(K + M)] <- equity_sizes
  }
  Y0
}


# ============================================================================
# 3. Scenario engine ---------------------------------------------------------
# ============================================================================

# For each of the four shock types, three small helpers describe how that
# shock acts on liabilities, bonds, and equities respectively. This dispatch
# table replaces the deeply-nested if/else trees of the original script.

shock_liability_value <- function(shock, ctx) {
  switch(as.character(shock),
    "1" = discounted_sum_liabilities(ctx$ts_up,   ctx$L0, ctx$dt, ctx$lambda),
    "2" = discounted_sum_liabilities(ctx$ts_down, ctx$L0, ctx$dt, ctx$lambda),
    "3" = discounted_sum_liabilities(ctx$ts,      ctx$L0, ctx$dt, ctx$lambda_new),
    "4" = discounted_sum_liabilities(ctx$ts,      ctx$L0, ctx$dt, ctx$lambda),
    NA_real_)
}

shock_bond_value <- function(shock, bond_cf, ctx) {
  switch(as.character(shock),
    "1" = discounted_sum_bond(ctx$ts_up,   bond_cf, ctx$dt),
    "2" = discounted_sum_bond(ctx$ts_down, bond_cf, ctx$dt),
    "3" = discounted_sum_bond(ctx$ts,      bond_cf, ctx$dt),  # bonds unaffected
    "4" = discounted_sum_bond(ctx$ts,      bond_cf, ctx$dt),  # bonds unaffected
    NA_real_)
}

shock_equity_value <- function(shock, base, mean_bond_pchange) {
  switch(as.character(shock),
    "1" = base,                                     # interest UP: no equity reaction
    "2" = base * (1 + mean_bond_pchange * 0.5),     # interest DOWN: half-correlation w/ bonds
    "3" = base,                                     # mortality:    no equity reaction
    "4" = base * 0.61,                              # Type-1 equity: -39% per Solvency II
    NA_real_)
}

#' Generate scenarios and stress them
#'
#' @param K,M           number of bonds / Type-1 equities
#' @param Y             (K + M) x horizon cashflow matrix from
#'                      `build_investment_universe`
#' @param N             number of states to simulate
#' @param J             number of shocks to apply (1..4)
#' @param exclude_shocks integer vector of shock IDs to leave at zero
#'
#' @return list with
#'   - economic_states  (K + M) x N matrix of unstressed asset values at T = 1
#'   - shocked_assets   (K + M) x N x J array of stressed asset values
#'   - shocked_liab     N x J matrix of stressed liability values
generate_scenarios <- function(K, Y, N, lower, upper, horizon, J,
                               i_magnitude, lambda, L0, mort_shock_size,
                               M = 0, exclude_shocks = integer(0)) {

  if (K + M == 0L) {
    warning("Investment universe is empty.")
    return(list(economic_states = matrix(0, 0, N),
                shocked_assets  = array(0, c(0, N, max(J, 1))),
                shocked_liab    = matrix(0, N, max(J, 1))))
  }
  # Out-of-range J: mirror the original's silent no-op behaviour.
  if (J < 1L || J > 4L) {
    return(list(economic_states = matrix(0, K + M, N),
                shocked_assets  = array(0, c(K + M, N, max(J, 1))),
                shocked_liab    = matrix(0, N, max(J, 1))))
  }

  active_shocks <- setdiff(seq_len(J), exclude_shocks)
  lambda_new    <- lambda * (1 - mort_shock_size)
  valuation_time <- 1L

  ESMat          <- matrix(0, nrow = K + M, ncol = N)
  shocked_assets <- array(0, dim = c(K + M, N, J))
  shocked_liab   <- matrix(0, nrow = N, ncol = J)

  for (i in seq_len(N)) {
    # --- Sample the state ------------------------------------------------
    ts <- draw_interest_term_structure(lower, upper, horizon)
    eq_returns <- if (M > 0) rlnorm(M, meanlog = 2, sdlog = 0.4) / 100 + 1
                  else numeric(0)

    ctx <- list(
      ts         = ts,
      ts_up      = shock_term_structure(ts, i_magnitude, "up"),
      ts_down    = shock_term_structure(ts, i_magnitude, "down"),
      L0         = L0,
      dt         = valuation_time,
      lambda     = lambda,
      lambda_new = lambda_new
    )

    # --- Unstressed asset values at T = 1 --------------------------------
    if (K > 0) {
      for (k in seq_len(K)) {
        ESMat[k, i] <- discounted_sum_bond(ts, Y[k, ], valuation_time)
      }
    }
    if (M > 0) {
      for (m in seq_len(M)) {
        ESMat[K + m, i] <- Y[K + m, 1] * eq_returns[m]
      }
    }

    # --- Stressed liabilities --------------------------------------------
    for (j in active_shocks) {
      shocked_liab[i, j] <- shock_liability_value(j, ctx)
    }

    # --- Stressed bonds; track % change for equity correlation -----------
    pchange <- numeric(K)  # bond %-change under interest-rate-DOWN shock
    if (K > 0) {
      for (k in seq_len(K)) {
        for (j in active_shocks) {
          shocked_assets[k, i, j] <- shock_bond_value(j, Y[k, ], ctx)
        }
        if (SHOCK_INTEREST_DOWN %in% active_shocks && ESMat[k, i] != 0) {
          pchange[k] <- (shocked_assets[k, i, SHOCK_INTEREST_DOWN] -
                          ESMat[k, i]) / ESMat[k, i]
        }
      }
    }
    mean_pchange <- if (K > 0) mean(pchange) else 0

    # --- Stressed equities -----------------------------------------------
    if (M > 0) {
      for (m in seq_len(M)) {
        base <- ESMat[K + m, i]
        for (j in active_shocks) {
          shocked_assets[K + m, i, j] <- shock_equity_value(j, base, mean_pchange)
        }
      }
    }
  }

  list(economic_states = ESMat,
       shocked_assets  = shocked_assets,
       shocked_liab    = shocked_liab)
}


# ============================================================================
# 4. LP construction --------------------------------------------------------
# ============================================================================

#' Sample expected surplus per asset (objective coefficients)
calc_expected_surplus <- function(economic_states, N) {
  rowSums(economic_states, na.rm = TRUE) / N
}

#' Pack the (N * J + 1) constraints expected by lpSolve
#'
#' Risk constraints (rows 1..N*J): asset value >= liability value for each
#' (state, shock) pair. Budget constraint (final row): cost of portfolio
#' at T = 0 must not exceed `max_budget`.
build_constraint_inputs <- function(N, K, J, Y0,
                                    shocked_assets, shocked_liab,
                                    max_budget, M = 0) {
  n_constr <- N * J + 1L
  n_vars   <- K + M

  C   <- matrix(0, nrow = n_constr, ncol = n_vars)
  rhs <- numeric(n_constr)
  dir <- rep(">=", n_constr); dir[n_constr] <- "<="

  # Risk constraints, block by shock.
  for (j in seq_len(J)) {
    for (i in seq_len(N)) {
      C[(j - 1L) * N + i, ]  <- shocked_assets[, i, j]
      rhs[(j - 1L) * N + i]  <- shocked_liab[i, j]
    }
  }
  # Budget constraint at T = 0.
  C[n_constr, ] <- Y0
  rhs[n_constr] <- max_budget

  list(mat = C, dir = dir, rhs = rhs)
}

#' Build and solve the LP via lpSolve
#'
#' lpSolve's `lp()` constructs and solves the model in a single call. It
#' assumes all decision variables are non-negative by default — exactly our
#' "no short-selling, no credit lines" requirement — so no explicit bounds
#' are needed. The objective is negated so that minimising `-obj * x` is
#' equivalent to maximising `obj * x`; downstream code un-negates `$objval`.
#'
#' Returns the solved `lp` object, with fields including:
#'   $status   solver status code (0 = optimal, 2 = infeasible, ...)
#'   $objval   optimal objective value (negated; see above)
#'   $solution optimal decision-variable values
build_lp_model <- function(objective, const_mat, const_dir, const_rhs,
                           N, J, K, M = 0) {
  lp(direction    = "min",
     objective.in = -1 * objective,
     const.mat    = const_mat,
     const.dir    = const_dir,
     const.rhs    = const_rhs)
}


# ============================================================================
# 5. Diagnostics -------------------------------------------------------------
# ============================================================================

#' Mark which constraints — including non-negativity bounds — bind at the
#' optimum. Returns a 0/1 vector of length N*J + 1 + (K + M).
find_binding_constraints <- function(const_mat, const_rhs, N, J, allocation) {
  n_explicit <- N * J + 1L
  binding    <- integer(n_explicit + length(allocation))
  lhs        <- as.numeric(const_mat %*% allocation)

  # FIX: the original used `for (i in 1:N*J)`, which evaluates as
  # `(1:N) * J` and silently skipped most constraints. Use seq_len(N*J).
  for (i in seq_len(N * J)) {
    if (round(lhs[i], 4) == round(const_rhs[i], 4)) binding[i] <- 1L
  }
  if (round(lhs[n_explicit], 0) == round(const_rhs[n_explicit], 0)) {
    binding[n_explicit] <- 1L
  }
  for (k in seq_along(allocation)) {
    if (allocation[k] == 0) binding[n_explicit + k] <- 1L
  }
  binding
}

#' Does an existing allocation still satisfy a (new) constraint set?
#'
#' Returns an indicator vector of length N*J + 1 with 1 in positions that
#' fail. Logs a one-line summary.
check_allocation_feasibility <- function(allocation, N, J,
                                         const_mat, const_rhs) {
  fail <- integer(N * J + 1L)
  lhs  <- as.numeric(const_mat %*% allocation)

  ok <- TRUE
  # FIX: original `for (i in 1:N*J)` operator-precedence bug; see above.
  for (i in seq_len(N * J)) {
    if (round(lhs[i], 4) < round(const_rhs[i], 4)) {
      ok <- FALSE
      fail[i] <- 1L
    }
  }
  bi <- N * J + 1L
  if (round(lhs[bi], 0) > round(const_rhs[bi], 0)) {
    ok <- FALSE
    fail[bi] <- 1L
  }
  message(if (ok) "Allocation satisfies the new constraint set."
                else "Allocation FAILS at least one new constraint.")
  fail
}

#' Translate an lpSolve status code, with audible cue
report_solver_status <- function(code) {
  msgs <- c(
    "0" = "Solved to optimality.",
    "1" = "Sub-optimal solution returned.",
    "2" = "Model is infeasible.",
    "3" = "Model is unbounded.",
    "4" = "Model is degenerate.",
    "5" = "Numerical failure encountered."
  )
  text <- msgs[as.character(code)]
  # FIX: original referenced `Result4` (a typo) in its catch-all branch.
  if (is.na(text)) text <- sprintf("Unknown solver status (code %s).", code)
  message(text)
  if (identical(code, 0L) || identical(code, 0)) beep(8) else beep(7)
  invisible(text)
}


# ============================================================================
# 6. End-to-end driver -------------------------------------------------------
# ============================================================================

#' One call that runs the full pipeline (scenarios -> objective ->
#' constraints -> LP -> solve) for a given configuration.
solve_portfolio <- function(K, M, Y, N, J, exclude_shocks = integer(0),
                            i_magnitude, lambda, L0, mort_shock_size,
                            ts_initial, max_budget,
                            lower, upper, horizon) {

  scen <- generate_scenarios(K = K, Y = Y, N = N,
                             lower = lower, upper = upper, horizon = horizon,
                             J = J, i_magnitude = i_magnitude,
                             lambda = lambda, L0 = L0,
                             mort_shock_size = mort_shock_size,
                             M = M, exclude_shocks = exclude_shocks)

  obj  <- calc_expected_surplus(scen$economic_states, N)
  Y0   <- calc_initial_values(K, Y, ts_initial, M)
  con  <- build_constraint_inputs(N, K, J, Y0,
                                  scen$shocked_assets, scen$shocked_liab,
                                  max_budget, M)
  lp     <- build_lp_model(obj, con$mat, con$dir, con$rhs, N, J, K, M)
  status <- lp$status   # lpSolve's lp() builds AND solves in one call

  list(
    config      = list(K = K, M = M, N = N, J = J,
                       max_budget = max_budget,
                       exclude_shocks = exclude_shocks),
    scenarios   = scen,
    objective   = obj,
    Y0          = Y0,
    constraints = con,
    lp          = lp,
    status      = status
  )
}

#' Extract solution summary statistics from a solved model
summarise_solution <- function(model) {
  if (model$status != 0L) {
    return(list(status = model$status,
                allocation = NA, expected_return = NA,
                amount_invested = NA, percentage_return = NA))
  }
  alloc            <- model$lp$solution
  expected_return  <- -1 * model$lp$objval                # un-negate
  amount_invested  <- as.numeric(model$Y0 %*% alloc)
  pct_return       <- (expected_return - amount_invested) / amount_invested
  list(status            = model$status,
       allocation        = alloc,
       expected_return   = expected_return,
       amount_invested   = amount_invested,
       percentage_return = pct_return)
}


###############################################################################
# === ANALYSIS ================================================================
#
# The remainder of this script reproduces the four model variants of the
# thesis. Each block is self-contained: it constructs its own configuration,
# solves the LP, and stores the result under a clearly-named variable.
###############################################################################


# ============================================================================
# 7. §4.1 – §4.2 : Toy Model -------------------------------------------------
# ============================================================================
# Three randomly-generated bonds, no equities; two shocks (interest UP/DOWN);
# N = 3 economic states; budget = 1000.

cfg_toy <- list(
  K = 3, M = 0,
  mat_lower = 1, mat_upper = 30,
  coupon_lower = 1, coupon_upper = 10,
  fv_lower = 1000, fv_upper = 2000,
  N = 3,
  lower = 0, upper = 0.1,    # interest-rate bounds
  J = 2, i_magnitude = 0.01,
  lambda = 1.375, L0 = 11750,
  mort_shock_size = 0.15,
  max_budget = 1000
)

# Investment universe + initial term structure are drawn here and reused
# (where appropriate) by later experiments — matches the thesis structure.
Bonds <- build_investment_universe(
  K = cfg_toy$K, mat_lower = cfg_toy$mat_lower, mat_upper = cfg_toy$mat_upper,
  coupon_lower = cfg_toy$coupon_lower, coupon_upper = cfg_toy$coupon_upper,
  fv_lower = cfg_toy$fv_lower, fv_upper = cfg_toy$fv_upper,
  M = cfg_toy$M, equity_sizes = numeric(0)
)
ts_initial <- draw_interest_term_structure(cfg_toy$lower, cfg_toy$upper,
                                           cfg_toy$mat_upper)

toy_model <- solve_portfolio(
  K = cfg_toy$K, M = cfg_toy$M, Y = Bonds,
  N = cfg_toy$N, J = cfg_toy$J,
  i_magnitude = cfg_toy$i_magnitude,
  lambda = cfg_toy$lambda, L0 = cfg_toy$L0,
  mort_shock_size = cfg_toy$mort_shock_size,
  ts_initial = ts_initial, max_budget = cfg_toy$max_budget,
  lower = cfg_toy$lower, upper = cfg_toy$upper, horizon = cfg_toy$mat_upper
)
report_solver_status(toy_model$status)
toy_summary <- summarise_solution(toy_model)


# ============================================================================
# 8. §4.3 : Binding constraints in the Toy Model -----------------------------
# ============================================================================

binding_toy <- find_binding_constraints(
  toy_model$constraints$mat, toy_model$constraints$rhs,
  cfg_toy$N, cfg_toy$J, toy_summary$allocation
)


# ============================================================================
# 9. §4.4 : Does the old allocation satisfy a freshly-drawn constraint set? --
# ============================================================================
# Re-run the same configuration (re-drawing scenarios) and check whether the
# original allocation remains feasible under the new constraints.

toy_model_redraw <- solve_portfolio(
  K = cfg_toy$K, M = cfg_toy$M, Y = Bonds,
  N = cfg_toy$N, J = cfg_toy$J,
  i_magnitude = cfg_toy$i_magnitude,
  lambda = cfg_toy$lambda, L0 = cfg_toy$L0,
  mort_shock_size = cfg_toy$mort_shock_size,
  ts_initial = ts_initial, max_budget = cfg_toy$max_budget,
  lower = cfg_toy$lower, upper = cfg_toy$upper, horizon = cfg_toy$mat_upper
)
toy_redraw_summary <- summarise_solution(toy_model_redraw)

# How does the new allocation perform under the OLD objective and vice versa?
return_new_alloc_old_obj <- sum(toy_redraw_summary$allocation * toy_model$objective)
return_old_alloc_new_obj <- sum(toy_summary$allocation       * toy_model_redraw$objective)

# Does the old allocation remain feasible under the new constraints?
feasibility_old_in_new <- check_allocation_feasibility(
  toy_summary$allocation, cfg_toy$N, cfg_toy$J,
  toy_model_redraw$constraints$mat, toy_model_redraw$constraints$rhs
)


# ============================================================================
# 10. §5.1 : Stability with N = 200 ------------------------------------------
# ============================================================================
# Hypothesis: a larger sample of economic states should produce more stable
# expected-surplus estimates, and therefore more stable optimal allocations.

cfg_big <- modifyList(cfg_toy, list(N = 200))

big_model_a <- solve_portfolio(
  K = cfg_big$K, M = cfg_big$M, Y = Bonds,
  N = cfg_big$N, J = cfg_big$J,
  i_magnitude = cfg_big$i_magnitude,
  lambda = cfg_big$lambda, L0 = cfg_big$L0,
  mort_shock_size = cfg_big$mort_shock_size,
  ts_initial = ts_initial, max_budget = cfg_big$max_budget,
  lower = cfg_big$lower, upper = cfg_big$upper, horizon = cfg_big$mat_upper
)
big_model_b <- solve_portfolio(
  K = cfg_big$K, M = cfg_big$M, Y = Bonds,
  N = cfg_big$N, J = cfg_big$J,
  i_magnitude = cfg_big$i_magnitude,
  lambda = cfg_big$lambda, L0 = cfg_big$L0,
  mort_shock_size = cfg_big$mort_shock_size,
  ts_initial = ts_initial, max_budget = cfg_big$max_budget,
  lower = cfg_big$lower, upper = cfg_big$upper, horizon = cfg_big$mat_upper
)
big_summary_a <- summarise_solution(big_model_a)
big_summary_b <- summarise_solution(big_model_b)


# ============================================================================
# 11. §5.2 : Adding the mortality shock --------------------------------------
# ============================================================================
# J is increased to 3 (interest UP / DOWN / mortality). The mortality shock
# enters with negative size (mortality DECREASE -> liabilities increase).
# At budget = 1000 the model is infeasible; raising the budget to 3000
# restores feasibility.

cfg_mort <- modifyList(cfg_big, list(
  J = 3,
  mort_shock_size = -0.15
))
ts_initial_mort <- draw_interest_term_structure(cfg_mort$lower, cfg_mort$upper,
                                                cfg_mort$mat_upper)

mort_model_low_budget <- solve_portfolio(
  K = cfg_mort$K, M = cfg_mort$M, Y = Bonds,
  N = cfg_mort$N, J = cfg_mort$J,
  i_magnitude = cfg_mort$i_magnitude,
  lambda = cfg_mort$lambda, L0 = cfg_mort$L0,
  mort_shock_size = cfg_mort$mort_shock_size,
  ts_initial = ts_initial_mort, max_budget = 1000,
  lower = cfg_mort$lower, upper = cfg_mort$upper, horizon = cfg_mort$mat_upper
)

# Reuse the same scenario draw for the higher-budget run by rebuilding only
# the constraint inputs and the LP — replicates the original §5.2 logic.
mort_constraints_hi <- build_constraint_inputs(
  cfg_mort$N, cfg_mort$K, cfg_mort$J, mort_model_low_budget$Y0,
  mort_model_low_budget$scenarios$shocked_assets,
  mort_model_low_budget$scenarios$shocked_liab,
  max_budget = 3000, M = cfg_mort$M
)
mort_lp_hi <- build_lp_model(
  mort_model_low_budget$objective,
  mort_constraints_hi$mat, mort_constraints_hi$dir, mort_constraints_hi$rhs,
  cfg_mort$N, cfg_mort$J, cfg_mort$K, cfg_mort$M
)
mort_status_hi <- mort_lp_hi$status

mort_alloc_hi    <- mort_lp_hi$solution
mort_return_hi   <- -1 * mort_lp_hi$objval
mort_invested_hi <- as.numeric(mort_model_low_budget$Y0 %*% mort_alloc_hi)
mort_pct_hi      <- (mort_return_hi - mort_invested_hi) / mort_invested_hi


# ============================================================================
# 12. §5.3 : Adding equities and the equity shock ----------------------------
# ============================================================================
# Add two Type-1 equity indices and the J = 4 equity shock. Mortality shock
# is excluded from this experiment (matches original `ExcludeShocks5`).

cfg_eq <- list(
  K = 3, M = 2,
  mat_upper = cfg_toy$mat_upper,
  N = 200, J = 4,
  lower = 0, upper = 0.1,
  i_magnitude = 0.01, mort_shock_size = 0.2,
  lambda = 1.375, L0 = 11750,
  max_budget = 1500,
  exclude_shocks = c(SHOCK_MORTALITY)   # original used c(0,3); only 3 is meaningful
)
equity_sizes_eq <- rep(1, cfg_eq$M)

# Reuse bond cashflows from §4 and append equity rows — same construction
# as in the original (§5.3) script.
universe_eq <- matrix(0, nrow = cfg_eq$K + cfg_eq$M, ncol = cfg_eq$mat_upper)
universe_eq[seq_len(cfg_eq$K), ] <- Bonds[seq_len(cfg_eq$K), ]
universe_eq[(cfg_eq$K + 1):(cfg_eq$K + cfg_eq$M), 1] <- equity_sizes_eq

ts_initial_eq <- draw_interest_term_structure(cfg_eq$lower, cfg_eq$upper,
                                              cfg_eq$mat_upper)

eq_model <- solve_portfolio(
  K = cfg_eq$K, M = cfg_eq$M, Y = universe_eq,
  N = cfg_eq$N, J = cfg_eq$J, exclude_shocks = cfg_eq$exclude_shocks,
  i_magnitude = cfg_eq$i_magnitude,
  lambda = cfg_eq$lambda, L0 = cfg_eq$L0,
  mort_shock_size = cfg_eq$mort_shock_size,
  ts_initial = ts_initial_eq, max_budget = cfg_eq$max_budget,
  lower = cfg_eq$lower, upper = cfg_eq$upper, horizon = cfg_eq$mat_upper
)
eq_summary <- summarise_solution(eq_model)


# ============================================================================
# 13. §5.4 : Sensitivity of the optimal allocation to the budget -------------
# ============================================================================
# Sweep the budget from 800 upward in 2 % increments and record the share
# allocated to bonds vs. equities at each step.

n_budget_steps <- 20L
allocation_grid <- matrix(0, nrow = n_budget_steps, ncol = cfg_eq$K + cfg_eq$M)
bond_share_seq   <- numeric(0)
equity_share_seq <- numeric(0)
budget_seq       <- numeric(0)

for (k in seq_len(n_budget_steps)) {
  budget_k <- 800 * (1.02 ^ (k - 1L))
  message(sprintf("Budget step %2d : max budget = %.2f", k, budget_k))

  # Each step uses freshly-drawn equity rows and a fresh initial term
  # structure (matching the original loop). Bond rows are reused.
  universe_k <- matrix(0, nrow = cfg_eq$K + cfg_eq$M, ncol = cfg_eq$mat_upper)
  universe_k[seq_len(cfg_eq$K), ] <- Bonds[seq_len(cfg_eq$K), ]
  universe_k[(cfg_eq$K + 1):(cfg_eq$K + cfg_eq$M), 1] <- equity_sizes_eq

  ts_initial_k <- draw_interest_term_structure(cfg_eq$lower, cfg_eq$upper,
                                               cfg_eq$mat_upper)

  model_k <- solve_portfolio(
    K = cfg_eq$K, M = cfg_eq$M, Y = universe_k,
    N = cfg_eq$N, J = cfg_eq$J, exclude_shocks = cfg_eq$exclude_shocks,
    i_magnitude = cfg_eq$i_magnitude,
    lambda = cfg_eq$lambda, L0 = cfg_eq$L0,
    mort_shock_size = cfg_eq$mort_shock_size,
    ts_initial = ts_initial_k, max_budget = budget_k,
    lower = cfg_eq$lower, upper = cfg_eq$upper, horizon = cfg_eq$mat_upper
  )

  if (model_k$status == 0L) {
    alloc_k <- model_k$lp$solution
    allocation_grid[k, ] <- alloc_k

    bonds_idx  <- seq_len(cfg_eq$K)
    equity_idx <- (cfg_eq$K + 1L):(cfg_eq$K + cfg_eq$M)

    bond_share_seq   <- c(bond_share_seq,
                          (model_k$Y0[bonds_idx]  %*% alloc_k[bonds_idx])  / budget_k)
    equity_share_seq <- c(equity_share_seq,
                          (model_k$Y0[equity_idx] %*% alloc_k[equity_idx]) / budget_k)
    budget_seq       <- c(budget_seq, budget_k)
  }
}


# ============================================================================
# 14. Summary printout -------------------------------------------------------
# ============================================================================

cat("\n=== §4: Toy Model =========================================\n")
print(toy_model$lp)
cat("Optimal allocation : "); print(toy_summary$allocation)
cat("Expected return    : ", toy_summary$expected_return,    "\n")
cat("Pct. return        : ", toy_summary$percentage_return,  "\n")

cat("\n=== §4.3: Binding constraints (toy model) ================\n")
cat("1 = binding, 0 = slack. Final (K + M) entries are the\n",
    "non-negativity constraints on the decision variables.\n", sep = "")
print(binding_toy)

cat("\n=== §4.4: Old allocation under new constraints ===========\n")
cat("Allocation (original) : "); print(toy_summary$allocation)
cat("Allocation (re-drawn) : "); print(toy_redraw_summary$allocation)
cat("Expected return       : ", toy_redraw_summary$expected_return, "\n")
cat("Pct. return           : ", toy_redraw_summary$percentage_return, "\n")
cat("Old allocation feasibility under new constraints: \n")
print(feasibility_old_in_new)

cat("\n=== §5.1: Stability with N = 200 =========================\n")
cat("Allocation A : "); print(big_summary_a$allocation)
cat("Allocation B : "); print(big_summary_b$allocation)
cat("Returns      : ", big_summary_a$expected_return, " / ",
                       big_summary_b$expected_return, "\n")
cat("Pct. returns : ", big_summary_a$percentage_return, " / ",
                       big_summary_b$percentage_return, "\n")

cat("\n=== §5.2: Mortality shock ================================\n")
cat("Status at budget = 1000 : ", mort_model_low_budget$status,
    "  (2 = infeasible)\n")
cat("Status at budget = 3000 : ", mort_status_hi, "\n")
cat("Allocation     : "); print(mort_alloc_hi)
cat("Expected return: ", mort_return_hi, "\n")
cat("Pct. return    : ", mort_pct_hi,    "\n")

cat("\n=== §5.3: Equities + equity shock ========================\n")
cat("Allocation     : "); print(eq_summary$allocation)
cat("Expected return: ", eq_summary$expected_return,   "\n")
cat("Pct. return    : ", eq_summary$percentage_return, "\n")

cat("\n=== §5.4: Budget sensitivity =============================\n")
print(allocation_grid)
cat("Budget grid   : "); print(budget_seq)
cat("Bond shares   : "); print(as.numeric(bond_share_seq))
cat("Equity shares : "); print(as.numeric(equity_share_seq))

# Sensitivity charts (matching the thesis figures).
print(qplot(budget_seq, as.numeric(bond_share_seq),
            xlab = "Budget", ylab = "Bond share of portfolio"))
print(qplot(budget_seq, as.numeric(equity_share_seq),
            xlab = "Budget", ylab = "Equity share of portfolio"))
