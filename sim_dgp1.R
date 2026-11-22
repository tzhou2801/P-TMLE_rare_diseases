# =============================================================================
# sim_dgp1.R
#
# Simulation study, data-generating process (DGP) 1: Table 2 and Figure 1 of
#
#   Zhou T, Gruber S, Lee H, Lee W, Nie L, van der Laan M (2026). Improving
#   Efficiency and Power in Small Rare-Disease Trials Through Prespecified,
#   Data-Adaptive Covariate Adjustment. arXiv:2611.15071.
#
# Section, equation, table and figure numbers refer to arXiv:2611.15071.
#
# Each replication simulates a randomized trial with n = 50 participants
# (Section 6.1) and applies the five estimators of Table 1 in the effect
# setting (tau = -1.3) and in the null setting (tau = 0). Within a replication
# the two settings share the covariates, treatment assignments, errors and
# cross-validation folds; their outcomes differ by tau * A.
#
# Run from the folder that should receive the output:
#     Rscript sim_dgp1.R
#
# Requires the CRAN packages glmnet, future, future.apply and ggplot2. The
# results in the paper were produced with R 4.2.1, glmnet 4.1-4, future 1.33.0,
# future.apply 1.11.0 and ggplot2 3.4.2.
#
# Output, written to the working directory (<time> = start time of the run):
#     sim_dgp1_raw_<time>.csv      one row per replication, setting and estimator
#     sim_dgp1_table2_<time>.csv   the rows of Table 2
#     sim_dgp1_figure1_<time>.png  Figure 1
# The rows of Table 2 are also printed.
#
# Run time: about 8 hours for 1,000 replications with 6 workers.
# =============================================================================


# ---- settings ---------------------------------------------------------------
reps      <- 1:1000     # replications to run; replication i always uses random-number stream i
seed      <- 20260622L  # seed from which the streams of all replications are derived
n_workers <- NA         # parallel R workers; NA = all available cores but one


# ---- design (Sections 6.1 and 6.2) ------------------------------------------
N       <- 50L                          # participants per trial
G0      <- 2/3                          # known treatment mechanism g0(1 | W) of the 2:1 allocation
NFOLDS  <- 25L                          # V = 25 cross-validation folds of 2 participants
K_INNER <- 24L                          # folds for cross-validation within a set of training rows
TAU     <- c(effect = -1.3, null = 0)   # treatment effect psi_0 = tau in the two settings
COVARS  <- c("W1", "W2", "W3", "W4", "W5")


# ---- data-generating process (Section 6.1) -----------------------------------
# Covariates, eq. (18): W1 age, W2 sex, W3 years since first symptom (a
# Gamma(shape 3.2, scale 2.4) draw rounded to one decimal), W4 baseline miglustat
# use, W5 baseline R4DNPCCSS score.
# Treatment: 2:1 allocation by random permutation within each stratum of W4.
# Error: N(0, 1.3^2).
GAMMA_SHAPE <- 3.2
GAMMA_SCALE <- 2.4
SIGMA       <- 1.3

# The covariates, treatment assignments and errors of one trial. The random-number
# draws are made in this order; the only other draw of a replication is the fold
# assignment in run_one_rep().
sim_base <- function(n = N) {
  W1 <- 2L + rbinom(n, size = 17, prob = 0.52)                          # age
  W2 <- rbinom(n, size = 1, prob = 0.53)                                # sex
  W3 <- round(rgamma(n, shape = GAMMA_SHAPE, scale = GAMMA_SCALE), 1)   # years since first symptom
  W4 <- rbinom(n, size = 1, prob = 0.8)                                 # baseline miglustat use
  W5 <- sample(0:20, size = n, replace = TRUE)                          # baseline R4DNPCCSS score
  A  <- integer(n)
  for (s in c(0L, 1L)) {              # 2:1 allocation within each stratum of W4
    idx <- which(W4 == s); m <- length(idx)
    if (m == 0L) next
    n_trt <- round(2 * m / 3)
    lab   <- c(rep(1L, n_trt), rep(0L, m - n_trt))
    A[idx] <- lab[sample.int(m)]
  }
  eps <- rnorm(n, mean = 0, sd = SIGMA)
  data.frame(W1 = W1, W2 = W2, W3 = W3, W4 = W4, W5 = W5, A = A, eps = eps)
}

# Outcome under DGP 1, eq. (19): Y = mu0(W) + tau * A + error, with
#   mu0(W) = 1.9 - 2 (W4 - 0.8) + 0.2 (W3 - 7.68) - 0.3 (W5 - 10)
#          = 4.964 + 0.2 W3 - 2 W4 - 0.3 W5,
# the mean change from baseline under placebo.
MU0    <- 1.9
B4     <- -2.0
B3     <- 0.20
B5     <- -0.3
CENTER <- c(W3 = 7.68, W4 = 0.8, W5 = 10)
mu0 <- function(base)
  MU0 + B4 * (base$W4 - CENTER[["W4"]]) + B3 * (base$W3 - CENTER[["W3"]]) + B5 * (base$W5 - CENTER[["W5"]])
make_outcome <- function(base, tau) mu0(base) + tau * base$A + base$eps


# ---- outcome-regression library (Table 1; Section 6.2) ----------------------
# 33 candidates: the linear working models containing A and one of the 32 subsets
# of {W1, ..., W5}, from the unadjusted model to the model with all five
# covariates, and a LASSO regression on A and W1, ..., W5 with the intercept and A
# unpenalized. Every candidate is linear in A, so its Qbar(1, W) - Qbar(0, W) is
# its coefficient of A.

# covariate-set label: "unadjusted", "W4+W5", ...
set_label <- function(v) if (is.null(v) || length(v) == 0) "unadjusted" else paste(v, collapse = "+")
powerset <- function(x) {
  out <- list(character(0))
  for (k in seq_along(x)) out <- c(out, combn(x, k, simplify = FALSE))
  out
}
Q_LIBRARY <- c(
  lapply(powerset(COVARS), function(s) list(learner = "ols", vars = s, label = set_label(s))),
  list(list(learner = "lasso", vars = COVARS, label = "LASSO")))
NQ <- length(Q_LIBRARY)

# index of the linear working model with exactly this covariate set
q_index <- function(vars) {
  hit <- which(vapply(Q_LIBRARY, function(q) q$learner == "ols" && setequal(q$vars, vars), logical(1)))
  stopifnot(length(hit) == 1L)
  hit
}

# Least squares on a design matrix that contains the intercept and A. Aliased
# coefficients are set to 0.
fit_Q_ols <- function(X, Y) {
  b <- lm.fit(X, Y)$coefficients
  b[is.na(b)] <- 0
  b
}

# LASSO (X: columns A, W1, ..., W5) with the intercept and A unpenalized. The
# penalty is cv.glmnet's lambda.min for the given folds: rep_len(1:K_INNER, n)
# within the training rows of a fold, and the V outer folds on all n
# participants. Returns the coefficients of (1, A, W1, ..., W5) and the retained
# covariates.
fit_Q_lasso_cv <- function(X, Y, foldid) {
  pf <- c(0, rep(1, ncol(X) - 1L))
  cv <- suppressWarnings(glmnet::cv.glmnet(X, Y, alpha = 1, penalty.factor = pf, foldid = foldid,
                                           thresh = 1e-10))
  b  <- as.numeric(coef(cv, s = "lambda.min"))
  list(beta = b, retained = colnames(X)[-1L][b[-(1:2)] != 0])
}


# ---- treatment-mechanism library (Table 1; Section 6.2) ---------------------
# The known g0 and logistic LASSO regressions of A on the main terms W1, ..., W5
# and the hinge functions (W1 - m1)+, (W3 - m3)+ and (W5 - m5)+, where m1, m3 and
# m5 are the sample medians of W1, W3 and W5 in the simulated trial. The
# penalties are NLAMG values spaced evenly on the log scale from lambda_max (all
# coefficients zero) of the rows being fitted down to lambda_max / LADDER_FLOOR.
# Fitted probabilities are truncated to G_BOUNDS.
NLAMG        <- 8L
LADDER_FLOOR <- 100
G_BOUNDS     <- c(0.025, 0.975)
G_COEF_TOL   <- 1e-8   # a dictionary column is retained by a fit if |coefficient| > G_COEF_TOL

# n x 8 dictionary: W1, (W1 - m1)+, W2, W3, (W3 - m3)+, W4, W5, (W5 - m5)+
make_dictionary <- function(W, knots)
  cbind(W1 = W$W1, h_W1 = pmax(W$W1 - knots[["W1"]], 0),
        W2 = W$W2,
        W3 = W$W3, h_W3 = pmax(W$W3 - knots[["W3"]], 0),
        W4 = W$W4,
        W5 = W$W5, h_W5 = pmax(W$W5 - knots[["W5"]], 0))

lambda_ladder <- function(lmax, n_rung = NLAMG, floor_ratio = LADDER_FLOOR)
  exp(seq(log(lmax), log(lmax / floor_ratio), length.out = n_rung))

# Logistic LASSO of A on the dictionary at the penalty ladder of these rows.
# Predictions and coefficients are always taken at the ladder (s = ladder), which
# also covers a glmnet path that stops before the last penalty.
fit_g_lasso <- function(Xd, A) {
  path <- suppressWarnings(glmnet::glmnet(Xd, A, family = "binomial", alpha = 1))
  lad  <- lambda_ladder(path$lambda[1])
  fit  <- suppressWarnings(glmnet::glmnet(Xd, A, family = "binomial", alpha = 1, lambda = lad))
  list(fit = fit, ladder = lad)
}
truncate_g <- function(g, bounds = G_BOUNDS) pmin(pmax(g, bounds[1]), bounds[2])
# n x NLAMG matrix of fitted probabilities g(1 | W), one column per penalty
predict_g <- function(gl, Xnew, truncate = TRUE) {
  p <- suppressWarnings(predict(gl$fit, newx = Xnew, s = gl$ladder, type = "response"))
  if (truncate) truncate_g(p) else p
}
# dictionary columns retained at penalty r of the ladder
g_retained <- function(gl, r, tol = G_COEF_TOL) {
  cf <- as.matrix(suppressWarnings(coef(gl$fit, s = gl$ladder[r])))[-1L, 1L]
  names(cf)[abs(cf) > tol]
}


# ---- TMLE+SL (Sections 3 and 5; Algorithm 1) ---------------------------------
# cv_grid() fits every outcome candidate and every treatment-mechanism candidate
# on the training rows of each fold and evaluates, on the validation rows, the
# two parts of the efficient influence curve at the known g0,
#   DY = H_g0(A, W) {Y - Qbar*(A, W)}  and  DW = Qbar*(1, W) - Qbar*(0, W) - psi_v.
# Here Qbar* is the candidate's training-row fit after the linear fluctuation
# (eq. 8) with clever covariate H_g(A, W) = A / g(1 | W) - (1 - A) / g(0 | W), for
# g = g0 or an estimated g, and psi_v is the TMLE of the pair on the training
# rows (the mean of its Qbar*(1, W) - Qbar*(0, W) over those rows). The CV-eIC^2
# risk (eq. 17) of a pair is the mean of (DY + DW)^2 over all n participants.
# With g = g0, Qbar*(1, W) - Qbar*(0, W) is the same for every participant, so DW = 0.
#
# d: data frame with W1, ..., W5, A, Y; folds: fold of each row (values 1..V);
# X_g: the dictionary rows of d. Column j of the grid is the treatment-mechanism
# candidate: j = 1 the known g0, j = 1 + r the logistic LASSO at penalty r.
# store = TRUE also returns the pieces that proj_var_cf() needs.
cv_grid <- function(d, folds, X_g, g0 = G0, store = TRUE) {
  n <- nrow(d); A <- d$A; Y <- d$Y
  Xq <- cbind(1, A, as.matrix(d[, COVARS, drop = FALSE]))
  colnames(Xq) <- c("(Intercept)", "A", COVARS)
  nG <- 1L + NLAMG
  DY <- DW <- array(NA_real_, c(n, NQ, nG))
  D0_va <- matrix(NA_real_, n, NQ)
  fold_ids <- sort(unique(folds))
  D0_tr <- g_keep <- vector("list", max(fold_ids))
  Hg0 <- A / g0 - (1 - A) / (1 - g0)

  for (v in fold_ids) {
    itr <- which(folds != v); iva <- which(folds == v)
    n_tr <- length(itr); n_va <- length(iva)

    # outcome candidates: fitted values on the training and validation rows
    QAW_tr <- matrix(NA_real_, n_tr, NQ); QAW_va <- matrix(NA_real_, n_va, NQ)
    betaA  <- numeric(NQ)
    for (q in seq_len(NQ)) {
      cand <- Q_LIBRARY[[q]]
      if (cand$learner == "ols") {
        cols <- c("(Intercept)", "A", cand$vars)
        b <- fit_Q_ols(Xq[itr, cols, drop = FALSE], Y[itr])
        QAW_tr[, q] <- Xq[itr, cols, drop = FALSE] %*% b
        QAW_va[, q] <- Xq[iva, cols, drop = FALSE] %*% b
        betaA[q]    <- b[["A"]]
      } else {
        fl <- fit_Q_lasso_cv(Xq[itr, -1L, drop = FALSE], Y[itr],
                             foldid = rep_len(seq_len(K_INNER), n_tr))
        QAW_tr[, q] <- Xq[itr, , drop = FALSE] %*% fl$beta
        QAW_va[, q] <- Xq[iva, , drop = FALSE] %*% fl$beta
        betaA[q]    <- fl$beta[2L]
      }
    }

    # treatment-mechanism candidates: g0 in column 1, the logistic LASSO ladder in columns 2..nG
    G_tr <- matrix(g0, n_tr, nG); G_va <- matrix(g0, n_va, nG)
    gl <- fit_g_lasso(X_g[itr, , drop = FALSE], A[itr])
    G_tr[, 1L + seq_len(NLAMG)] <- predict_g(gl, X_g[itr, , drop = FALSE])
    G_va[, 1L + seq_len(NLAMG)] <- predict_g(gl, X_g[iva, , drop = FALSE])

    # linear fluctuation of every (outcome, treatment-mechanism) pair, fitted on the
    # training rows in closed form, and the influence-curve parts at g0
    R_tr <- Y[itr] - QAW_tr
    for (j in seq_len(nG)) {
      Htr <- A[itr] / G_tr[, j] - (1 - A[itr]) / (1 - G_tr[, j])
      Hva <- A[iva] / G_va[, j] - (1 - A[iva]) / (1 - G_va[, j])
      eps <- drop(crossprod(Htr, R_tr)) / sum(Htr^2)          # one fluctuation parameter per outcome candidate
      QA_va <- QAW_va + outer(Hva, eps)
      blip  <- matrix(betaA, n_va, NQ, byrow = TRUE) + outer(1 / (G_va[, j] * (1 - G_va[, j])), eps)
      psi_v <- colMeans(matrix(betaA, n_tr, NQ, byrow = TRUE) + outer(1 / (G_tr[, j] * (1 - G_tr[, j])), eps))
      DY[iva, , j] <- Hg0[iva] * (Y[iva] - QA_va)
      DW[iva, , j] <- sweep(blip, 2L, psi_v)
    }

    if (store) {
      # influence curve at g0 of each initial (untargeted) training-row fit, on the
      # validation and on the training rows, and the dictionary columns retained
      # by the training-row treatment-mechanism fit at each penalty
      D0_va[iva, ] <- Hg0[iva] * (Y[iva] - QAW_va)
      D0_tr[[v]]   <- Hg0[itr] * (Y[itr] - QAW_tr)
      g_keep[[v]]  <- lapply(seq_len(NLAMG), function(r) g_retained(gl, r))
    }
  }
  risk <- apply((DY + DW)^2, c(2L, 3L), mean)      # NQ x nG CV-eIC^2 risks
  out <- list(risk = risk, n = n, folds = folds)
  if (store) { out$DY <- DY; out$DW <- DW; out$D0_va <- D0_va; out$D0_tr <- D0_tr; out$g_keep <- g_keep }
  out
}

# Discrete Super Learner selection by the CV-eIC^2 risk (Algorithm 1, Step 1).
# TMLE+SL (g0): the outcome candidate with the smallest risk at g0 (Step 1a).
# TMLE+SL (gn): the same outcome candidate; then, given it, the treatment-
# mechanism candidate with the smallest risk (Step 1b). An estimated g replaces g0
# only if it lowers the risk by more than the relative tolerance TIE_TOL.
# which.min() takes the first minimum, so exact ties go to the candidate listed first.
TIE_TOL <- 1e-6
select_Q <- function(grid) {
  r <- grid$risk[, 1L]
  r[!is.finite(r)] <- Inf
  which.min(r)
}
select_Qg <- function(grid, tie_tol = TIE_TOL) {
  r <- grid$risk
  r[!is.finite(r)] <- Inf
  q <- which.min(r[, 1L])
  j <- which.min(r[q, ])
  if (j > 1L && !(r[q, j] < r[q, 1L] * (1 - tie_tol))) j <- 1L
  list(q = q, j = j, rung = if (j == 1L) NA_integer_ else j - 1L)
}

# Steps 2 and 3 of Algorithm 1 on all n participants: fit outcome candidate q,
# target it with the linear fluctuation (eq. 8), with the clever covariate at the
# known g0 (rung = NULL) or at the logistic LASSO fit at penalty `rung`
# (truncated), and average Qbar*(1, W) - Qbar*(0, W) over the participants.
# `sel` is the covariate set of the fit: that of the working model, or the
# covariates the LASSO retained.
tmle_point <- function(d, folds, q, rung = NULL, X_g = NULL, g0 = G0) {
  n <- nrow(d); A <- d$A; Y <- d$Y
  Xq <- cbind(1, A, as.matrix(d[, COVARS, drop = FALSE]))
  colnames(Xq) <- c("(Intercept)", "A", COVARS)
  cand <- Q_LIBRARY[[q]]
  if (cand$learner == "ols") {
    cols <- c("(Intercept)", "A", cand$vars)
    b    <- fit_Q_ols(Xq[, cols, drop = FALSE], Y)
    QAW  <- drop(Xq[, cols, drop = FALSE] %*% b); betaA <- b[["A"]]
    sel  <- cand$label
  } else {
    fl   <- fit_Q_lasso_cv(Xq[, -1L, drop = FALSE], Y, foldid = folds)
    QAW  <- drop(Xq %*% fl$beta); betaA <- fl$beta[2L]
    sel  <- set_label(fl$retained)
  }
  if (is.null(rung)) {
    g <- rep(g0, n)
  } else {
    gl    <- fit_g_lasso(X_g, A)
    g_raw <- predict_g(gl, X_g, truncate = FALSE)[, rung]
    g     <- truncate_g(g_raw)
  }
  H   <- A / g - (1 - A) / (1 - g)
  eps <- sum(H * (Y - QAW)) / sum(H^2)
  psi <- mean(betaA + eps / (g * (1 - g)))       # mean of Qbar*(1, W) - Qbar*(0, W)
  list(psi = psi, sel = sel, learner = cand$learner)
}

# Cross-validated influence-curve variance (Algorithm 1, Step 4; Section 5.2):
# var() of the n out-of-fold influence-curve values DY + DW of the pair (q, j),
# divided by n.
cv_ic_var <- function(grid, q, j = 1L) var(grid$DY[, q, j] + grid$DW[, q, j]) / grid$n

# Variance of TMLE+SL (gn), which accounts for the estimation of the treatment
# mechanism (Section 6.2): the variance of a cross-fitted projection of the
# influence curve. For each outer fold v:
#   1. TMLE+SL (gn) is selected again on the training rows of fold v, with the
#      cross-validation grid of those rows (inner_grids[[v]]);
#   2. the influence curve at g0 of the selected outcome candidate's initial fit
#      on the training rows of fold v, H_g0(A, W) {Y - Qbar(A, W)}, is regressed
#      over those training rows on the score S = (A - g0) (1, W_keep) of the
#      selected logistic model, where W_keep are the dictionary columns that the
#      training-row fit retains at the selected penalty (none, i.e. the
#      intercept score alone, when g0 is selected);
#   3. the residual of this regression is evaluated on the validation rows of fold v.
# The variance estimate is var() of the n cross-fitted residuals, divided by n.
proj_var_cf <- function(grid, inner_grids, d, X_g, g0 = G0) {
  n <- grid$n; folds <- grid$folds; A <- d$A
  r <- rep(NA_real_, n)
  for (v in seq_along(inner_grids)) {
    itr <- which(folds != v); iva <- which(folds == v)
    w   <- select_Qg(inner_grids[[v]])
    keep_v <- if (w$j == 1L) character(0) else grid$g_keep[[v]][[w$rung]]
    D0_tr <- grid$D0_tr[[v]][, w$q]
    D0_va <- grid$D0_va[iva, w$q]
    S_tr  <- (A[itr] - g0) * cbind(1, X_g[itr, keep_v, drop = FALSE])
    S_va  <- (A[iva] - g0) * cbind(1, X_g[iva, keep_v, drop = FALSE])
    cf <- lm.fit(S_tr, D0_tr)$coefficients
    cf[is.na(cf)] <- 0
    r[iva] <- D0_va - drop(S_va %*% cf)
  }
  var(r) / n
}

# Wald confidence interval and two-sided test at level alpha, with a Student t
# reference distribution with n - 2 degrees of freedom (Section 6.2).
inference_tdist <- function(truth, psi, var, n, alpha = 0.05) {
  se     <- sqrt(var)
  cutoff <- qt(alpha / 2, df = n - 2, lower.tail = FALSE)
  tstat  <- psi / se
  data.frame(
    psi.hat = psi, se = se,
    cover   = (psi - cutoff * se <= truth) & (truth <= psi + cutoff * se),
    reject  = 2 * pt(abs(tstat), df = n - 2, lower.tail = FALSE) < alpha
  )
}


# ---- the five estimators (Table 1) -------------------------------------------
est_row <- function(label, pt, var, tau, n = N) {
  inf <- inference_tdist(tau, pt$psi, var, n)
  data.frame(estimator = label, psi_hat = inf$psi.hat, se = inf$se, cover = inf$cover, reject = inf$reject,
             sel = pt$sel, learner = pt$learner, error = NA_character_, stringsAsFactors = FALSE)
}
# an estimator that fails in a replication gives a row of NA with the error message
safely <- function(expr, label) tryCatch(expr, error = function(e)
  data.frame(estimator = label, psi_hat = NA_real_, se = NA_real_, cover = NA, reject = NA,
             sel = NA_character_, learner = NA_character_, error = conditionMessage(e), stringsAsFactors = FALSE))

# 1-3. Unadjusted, ANCOVA(W4, W5), ANCOVA(W3, W4, W5): one linear working model,
# targeted at g0 (which leaves the estimate equal to its coefficient of A, up to
# rounding), with the cross-validated influence-curve variance.
est_fixed <- function(d, grid, vars, label, tau) {
  q <- q_index(vars)
  est_row(label, tmle_point(d, grid$folds, q), cv_ic_var(grid, q, 1L), tau)
}
# 4. TMLE+SL (g0): discrete Super Learner over the 33 outcome candidates at g0.
est_tmle_sl_g0 <- function(d, grid, label, tau) {
  q <- select_Q(grid)
  est_row(label, tmle_point(d, grid$folds, q), cv_ic_var(grid, q, 1L), tau)
}
# 5. TMLE+SL (gn): the outcome candidate of TMLE+SL (g0) and the discrete Super
# Learner over g0 and the logistic LASSO ladder, with the projection variance.
est_tmle_sl_gn <- function(d, grid, inner_grids, X_g, label, tau) {
  w  <- select_Qg(grid)
  pt <- if (w$j == 1L) tmle_point(d, grid$folds, w$q)
        else           tmle_point(d, grid$folds, w$q, rung = w$rung, X_g = X_g)
  est_row(label, pt, proj_var_cf(grid, inner_grids, d, X_g), tau)
}


# ---- one replication ---------------------------------------------------------
# Draws the trial and its folds, then computes all estimators in the effect
# setting and in the null setting. The estimators draw no random numbers.
run_one_rep <- function(i) {
  base  <- sim_base(N)
  folds <- sample(rep(seq_len(NFOLDS), length.out = N))
  # hinge knots: the sample medians (type-7 quantiles at 1/2) of W1, W3 and W5
  knots <- lapply(c(W1 = "W1", W3 = "W3", W5 = "W5"),
                  function(v) as.numeric(quantile(base[[v]], 0.5, type = 7, names = FALSE)))
  X_g <- make_dictionary(base, knots)

  one_setting <- function(st) {
    tau <- TAU[[st]]
    d   <- base[, c(COVARS, "A")]
    d$Y <- make_outcome(base, tau)
    grid  <- cv_grid(d, folds, X_g, store = TRUE)
    # the cross-validation grid within the training rows of each fold, for proj_var_cf()
    inner <- lapply(seq_len(NFOLDS), function(v) {
      itr <- which(folds != v)
      cv_grid(d[itr, , drop = FALSE], rep_len(seq_len(K_INNER), length(itr)), X_g[itr, , drop = FALSE],
              store = FALSE)
    })
    rows <- rbind(
      safely(est_fixed(d, grid, character(0),        "Unadjusted",       tau), "Unadjusted"),
      safely(est_fixed(d, grid, c("W4", "W5"),       "ANCOVA(W4,W5)",    tau), "ANCOVA(W4,W5)"),
      safely(est_fixed(d, grid, c("W3", "W4", "W5"), "ANCOVA(W3,W4,W5)", tau), "ANCOVA(W3,W4,W5)"),
      safely(est_tmle_sl_g0(d, grid,                 "TMLE+SL(g0)",      tau), "TMLE+SL(g0)"),
      safely(est_tmle_sl_gn(d, grid, inner, X_g,     "TMLE+SL(gn)",      tau), "TMLE+SL(gn)"))
    data.frame(rep = as.integer(i), setting = st, tau = tau, rows, stringsAsFactors = FALSE)
  }
  rbind(one_setting("effect"), one_setting("null"))
}


# ---- replications in parallel ------------------------------------------------
# Replication i uses its own L'Ecuyer-CMRG random-number stream: stream i of the
# sequence that future.apply generates for future.seed = TRUE after set.seed(seed).
# Its results therefore do not depend on the number of workers or on which other
# replications are run.
replicate_seeds <- function(n, seed) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  RNGkind("L'Ecuyer-CMRG")
  s <- get(".Random.seed", envir = globalenv())
  out <- vector("list", n)
  for (i in seq_len(n)) {
    out[[i]] <- parallel::nextRNGSubStream(s)
    s <- parallel::nextRNGStream(s)
  }
  RNGkind("Mersenne-Twister")
  out
}

run_simulation <- function(reps, seed, n_workers) {
  stopifnot(length(reps) >= 1L, all(reps >= 1L), !anyDuplicated(reps))
  oplan <- if (n_workers <= 1L) future::plan(future::sequential)
           else if (future::supportsMulticore()) future::plan(future::multicore, workers = n_workers)
           else future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(oplan), add = TRUE)
  seeds <- replicate_seeds(max(reps), seed)[reps]
  out <- future.apply::future_lapply(reps, run_one_rep, future.seed = seeds, future.packages = "glmnet")
  do.call(rbind, out)
}


# ---- Table 2 -------------------------------------------------------------------
# Bias, Monte Carlo SD, MSE (= bias^2 + SD^2), coverage of the 95% confidence
# interval and power from the effect setting; Type I error rate from the null
# setting (two-sided tests at the 5% level).
table_rows <- function(raw) {
  ok <- raw[is.na(raw$error), ]
  one <- function(e) {
    eff  <- ok[ok$setting == "effect" & ok$estimator == e, ]
    nul  <- ok[ok$setting == "null" & ok$estimator == e, ]
    bias <- mean(eff$psi_hat) - TAU[["effect"]]
    data.frame(estimator = e, bias = bias, sd = sd(eff$psi_hat), mse = bias^2 + var(eff$psi_hat),
               coverage = mean(eff$cover), power = mean(eff$reject), type1_error = mean(nul$reject),
               n_effect = nrow(eff), n_null = nrow(nul), stringsAsFactors = FALSE)
  }
  do.call(rbind, lapply(unique(raw$estimator), one))
}
print_table <- function(tab, title) {
  f3 <- function(x) formatC(x, digits = 3, format = "f")
  cat("\n", title, "\n", sep = "")
  print(data.frame(Estimator = tab$estimator, Bias = f3(tab$bias), SD = f3(tab$sd), MSE = f3(tab$mse),
                   Coverage = f3(tab$coverage), Power = f3(tab$power), `Type I error rate` = f3(tab$type1_error),
                   check.names = FALSE), row.names = FALSE)
}


# ---- Figure 1 ------------------------------------------------------------------
# Covariates in the outcome working model selected by TMLE+SL (g0), effect
# setting: the proportion of replications in which each covariate set was
# selected. A LASSO selection is counted by the covariates it retained; sets
# selected in fewer than 2% of the replications are omitted; dark bars mark the
# sets that include all three prognostic covariates W3, W4 and W5.
set_labeller <- function(x) {
  x <- as.character(x)
  is_l <- grepl("^LASSO: ", x); core <- sub("^LASSO: ", "", x)
  core <- ifelse(core == "unadjusted", "'unadjusted'", gsub("W([1-5])", "W[\\1]", gsub("+", "*'+'*", core, fixed = TRUE)))
  parse(text = ifelse(is_l, paste0("'LASSO: '*", core), core))
}
make_figure1 <- function(raw, file) {
  selection_theme <- theme_classic(base_size = 11) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
          axis.text.y = element_text(colour = "black"), panel.spacing.x = grid::unit(1.25, "lines"))
  R1 <- length(unique(raw$rep))
  s1 <- raw[raw$setting == "effect" & raw$estimator == "TMLE+SL(g0)" & is.na(raw$error), ]
  d1 <- as.data.frame(table(set = s1$sel), stringsAsFactors = FALSE)
  d1$rate <- d1$Freq / nrow(s1)
  d1 <- d1[d1$rate >= 0.02, ]
  d1$set <- factor(d1$set, levels = d1$set[order(d1$rate)])
  d1$prog <- factor(ifelse(vapply(strsplit(as.character(d1$set), "+", fixed = TRUE),
                                  function(v) all(c("W3", "W4", "W5") %in% v), logical(1)), "all", "some"),
                    levels = c("all", "some"))
  p1 <- ggplot(d1, aes(set, rate, fill = prog)) + geom_col(width = 0.7) + coord_flip() +
    scale_fill_manual(values = c(all = "grey25", some = "grey70"), drop = FALSE, name = NULL,
                      labels = c(all = "Includes all three prognostic covariates", some = "Omits a prognostic covariate")) +
    scale_x_discrete(labels = set_labeller) +
    scale_y_continuous(limits = c(0, max(0.5, ceiling(max(d1$rate) * 10) / 10)), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Covariates in the selected working model", y = sprintf("Proportion of %s Monte Carlo replications", format(R1, big.mark = ","))) +
    selection_theme + theme(legend.position = "bottom")
  ggsave(file, p1, width = 6.5, height = 3.9, dpi = 300, bg = "white")
  invisible(d1)
}


# ---- run the simulation ----------------------------------------------------------
for (pkg in c("glmnet", "future", "future.apply", "ggplot2"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Please install the R package '", pkg, "'.")
if (is.na(n_workers)) n_workers <- max(1L, future::availableCores() - 1L)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
cat(sprintf("DGP 1: %d replications with %d workers, started %s\n", length(reps), n_workers, format(Sys.time())))
t0  <- Sys.time()
raw <- run_simulation(reps, seed, n_workers)
cat(sprintf("finished after %.2f hours; %d estimator errors\n",
            as.numeric(difftime(Sys.time(), t0, units = "hours")), sum(!is.na(raw$error))))
write.csv(raw, sprintf("sim_dgp1_raw_%s.csv", stamp), row.names = FALSE)
tab <- table_rows(raw)
write.csv(tab, sprintf("sim_dgp1_table2_%s.csv", stamp), row.names = FALSE)
print_table(tab, "Table 2 (DGP 1)")
suppressPackageStartupMessages(library(ggplot2))
make_figure1(raw, sprintf("sim_dgp1_figure1_%s.png", stamp))
cat(sprintf("output written to %s\n", getwd()))
