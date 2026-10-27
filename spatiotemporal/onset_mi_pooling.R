# =============================================================================
# onset_mi_pooling.R — the multiple-imputation pooling estimators
# BDBV 2026 DRC · Spatiotemporal Invasion Forecast Suite
#
# Factored OUT of 48_onset_mi_pool.R so the arithmetic can be unit-tested
# (tests/test_onset_mi.R) without loading a sweep: 48 is a script that executes on source
# and needs M completed arms on disk, so its helpers were unreachable from a test.
#
# Two estimators are in use across this suite and they are NOT interchangeable:
#   * the DRAW MIXTURE (exact; used for the current forecast, where posterior draws exist) —
#     that lives in 21_bayesian_renewal.R::.mix_and_summarise_draws() and is reused, not
#     reimplemented, so one estimator pools the GT mixture and the imputation mixture alike;
#   * RUBIN'S RULES (used wherever the artifact is already a summary with a within-arm
#     variance) — that is what this file provides.
#
# REFERENCE: Rubin DB (1987) Multiple Imputation for Nonresponse in Surveys, §3.1 — the
# combining rules, the reference df nu = (M - 1) / lambda^2 and the fmi formula, which is
# everything implemented here. The Barnard & Rubin (1999, Biometrika 86:948-955) small-sample
# df is NOT used: it needs a complete-data df, which these summaries (bootstrap bands, AUCs,
# Monte-Carlo cascade means) do not have. An earlier header cited it as the df in use.
# =============================================================================

#' Rubin's rules for one estimand.
#'
#' Q_bar = mean(Q_m);  B = var(Q_m);  T = U_bar + (1 + 1/M) B.
#'
#' THE INTERVAL IS ON THE t DISTRIBUTION, NOT THE NORMAL. Rubin's reference distribution is
#' t with nu = (M - 1) / lambda^2 degrees of freedom. At M = 20 with an unknown within
#' component (lambda = 1) that is nu = 19 and t_{19, 0.95} = 1.729 against z_{0.95} = 1.645 —
#' so a normal interval is 5% too NARROW, on every number this sweep publishes. When B is
#' small relative to U_bar, nu -> Inf and the t collapses to the normal, which is why the
#' normal looks adequate until it is not.
#'
#' TWO DIFFERENT QUANTITIES ARE RETURNED AND THEY MUST NOT BE CONFLATED:
#'   lambda = (1 + 1/M) B / T   the PROPORTION OF TOTAL VARIANCE attributable to the
#'                              missing data. This is what a variance-decomposition figure
#'                              should show.
#'   fmi    = lambda + (1 - lambda) * 2/(nu + 3)   the FRACTION OF MISSING INFORMATION,
#'                              i.e. lambda with Rubin's small-sample df correction. It is
#'                              never below lambda, and equals it at lambda = 0 and at
#'                              lambda = 1 (where nu = M - 1 but the factor 1 - lambda is 0).
#'                              With nu = (M - 1)/lambda^2 the gap is
#'                              2 lambda^2 (1 - lambda) / (M - 1 + 3 lambda^2); at M = 20 its
#'                              maximum is 0.0146, at lambda = 0.652. Small here, but the two
#'                              are different statistics and are labelled as such.
#' Earlier revisions of this file returned lambda under the name `fmi`. Publishing the one
#' labelled as the other is a misstatement of a reported statistic, which is why both are
#' now returned and the caller chooses explicitly.
#'
#' @param q per-arm point estimates. Non-finite entries are dropped, and the corresponding
#'   entries of `u` are dropped with them so the pairing cannot slip.
#' @param u per-arm WITHIN variances, same length as `q`. NULL or all-NA is treated as ZERO
#'   and flagged in `u_known`: T then collapses to (1 + 1/M)B and lambda to 1. That is the
#'   honest reading for a quantity with no published within-arm interval — it says "we cannot
#'   separate the two components" — not a measurement that the within component is nil.
#' @return list(M, n_u, q, B, U, T, se, riv, lambda, fmi, df, lo, hi, u_known); the interval
#'   is 90%, matching this suite's band convention (90%/60%, never 95%/50%).
mi_rubin <- function(q, u = NULL) {
  q <- as.numeric(q)
  if (!is.null(u)) {
    u <- as.numeric(u)
    # A LENGTH MISMATCH IS FATAL, not something to recycle around. `u[ok]` on a vector of
    # the wrong length silently mis-pairs each arm's within variance with another arm's
    # estimate, which no downstream check could detect.
    if (length(u) != length(q))
      stop(sprintf("[mi_rubin] q has %d element(s) but u has %d; they must be parallel.",
                   length(q), length(u)), call. = FALSE)
  }
  ok <- is.finite(q)
  if (!is.null(u)) u <- u[ok]
  q <- q[ok]
  m <- length(q)
  .empty <- list(M = 0L, n_u = 0L, q = NA_real_, B = NA_real_, U = NA_real_, T = NA_real_,
                 se = NA_real_, riv = NA_real_, lambda = NA_real_, fmi = NA_real_,
                 df = NA_real_, lo = NA_real_, hi = NA_real_, u_known = FALSE)
  if (!m) return(.empty)

  qbar <- mean(q)
  # var() is NA at m = 1; a single completion carries NO between-imputation information, and
  # B = 0 is the correct statement of that (T then equals U_bar, i.e. the single-draw answer).
  B <- if (m > 1L) stats::var(q) else 0
  .ufin <- if (is.null(u)) logical(0) else is.finite(u)
  n_u <- sum(.ufin)
  ukn <- n_u > 0L
  Ubar <- if (ukn) mean(u[.ufin]) else 0
  Tv <- Ubar + (1 + 1/m) * B

  # lambda, then df from lambda (nu = (M-1)/lambda^2 — algebraically identical to Rubin's
  # (M-1)(1 + 1/riv)^2 but numerically safe at riv = 0 and riv = Inf, where the product form
  # is 0*Inf). At m = 1 there is no between-imputation df at all: nu = Inf, i.e. the ordinary
  # single-draw normal interval, which is exactly what one completion earns.
  lambda <- if (is.finite(Tv) && Tv > 0) (1 + 1/m) * B / Tv else NA_real_
  riv    <- if (Ubar > 0) (1 + 1/m) * B / Ubar else if (B > 0) Inf else 0
  df <- if (m <= 1L || !is.finite(lambda) || lambda <= 0) Inf else (m - 1) / lambda^2
  # fmi = lambda + (1-lambda) * 2/(nu+3). Reduces to lambda as nu -> Inf and to 1 at
  # lambda = 1, both correct.
  fmi <- if (!is.finite(lambda)) NA_real_
         else if (is.infinite(df)) lambda
         else lambda + (1 - lambda) * 2 / (df + 3)

  se <- sqrt(max(Tv, 0))
  tcrit <- mi_tcrit(df)          # ONE definition of the critical value, shared with callers
  list(M = m, n_u = n_u, q = qbar, B = B, U = Ubar, T = Tv, se = se,
       riv = riv, lambda = lambda, fmi = fmi, df = df,
       lo = qbar - tcrit * se, hi = qbar + tcrit * se, u_known = ukn)
}

#' The critical value mi_rubin() used, so a caller drawing a SECOND band (the within-arm
#' component alone, for a within-vs-total comparison) uses the same one.
#'
#' Drawing a t-based total beside a z-based within component would make the comparison
#' partly an artefact of two different reference distributions rather than of the
#' between-imputation variance it exists to show.
#' VECTORISED, because its callers use it inside mutate()/transmute() over a whole column of
#' per-row degrees of freedom. A scalar-only `if (!is.finite(df))` aborted there with
#' "the condition has length > 1".
mi_tcrit <- function(df) {
  df <- as.numeric(df)
  # ifelse() evaluates both arms over the full vector; qt(0.95, Inf) is finite in R
  # (it returns the normal quantile), so the unused branch cannot error.
  ifelse(is.finite(df) & df > 0, stats::qt(0.95, df = df), stats::qnorm(0.95))
}

#' Recover a within-arm VARIANCE from a published central percentile-bootstrap interval.
#'
#' evaluate_invasion() (16_invasion_eval.R: quantile(v, c(0.05, 0.95))) and
#' compute_detection_curve() (20_forecast_detail.R: .detc_boot_ci) publish 5%/95% bootstrap
#' quantiles — a 90% interval — not standard errors, so U is recovered as
#' ((hi - lo) / (2 z_{(1+level)/2}))^2. THIS IS A NORMAL APPROXIMATION to the bootstrap
#' spread, and z (not t) is correct here: it inverts an interval that was itself formed as a
#' percentile pair, not as estimate +/- t*se. It touches only U_bar; the between component B —
#' which is what the sweep exists to measure — is exact.
#' @param level the central coverage of the interval being inverted (0.90 for every interval
#'   this suite publishes). Passing the wrong level scales U by (z_wrong / z_right)^2.
mi_u_from_ci <- function(lo, hi, level = 0.90) {
  stopifnot(is.numeric(level), length(level) == 1L, level > 0, level < 1)
  w <- (as.numeric(hi) - as.numeric(lo)) / (2 * stats::qnorm((1 + level) / 2))
  ifelse(is.finite(w) & w >= 0, w^2, NA_real_)
}

#' The metric STEMS of a frame that carry a published within-arm interval.
#'
#' A stem qualifies when the frame has the stem itself AND both `<stem>_lo` and `<stem>_hi`.
#' Derived from the columns rather than listed, because evaluate_invasion() gains interval
#' columns over time (auc_pr_skill_date_lo/hi were added after auc_pr and log_score): a fixed
#' list silently pooled the new bounds as metrics of their own and gave their estimand U = 0.
#' @param nms column names.
#' @return character vector of stems (possibly empty), in column order.
mi_ci_stems <- function(nms) {
  nms <- as.character(nms)
  lo <- nms[grepl("_lo$", nms)]
  stems <- sub("_lo$", "", lo)
  stems[stems %in% nms & paste0(stems, "_hi") %in% nms]
}

#' Per-row MI weights: equal total weight per completion WITHIN each slice.
#'
#' Arms score different numbers of rows (the at-risk support moves with the imputed labels),
#' so a consumer averaging the stacked rows of one slice, e.g. within (method, horizon), would
#' weight each completion by how many rows it contributed there. The weight is 1 / (M * n_scored) with n_scored counted
#' WITHIN (arm, slice), so every arm's scored rows sum to 1/M in every slice. Normalising over
#' an arm's rows across ALL slices instead gives each arm equal weight only in the stack as a
#' whole, not in the slice a consumer actually averages.
#' @param scored logical, the rows a consumer scores (others get NA).
#' @param arm per-row completion id.
#' @param slice a list of per-row slice keys (e.g. list(method, horizon)).
#' @param M number of completions.
mi_slice_weights <- function(scored, arm, slice, M) {
  scored <- scored %in% TRUE
  stopifnot(length(arm) == length(scored), is.list(slice),
            all(vapply(slice, length, integer(1)) == length(scored)),
            length(M) == 1L, is.finite(M), M >= 1)
  # paste(), not interaction(): an NA key then forms its own group instead of being dropped
  # by split() and leaving its rows with an unnormalised count.
  g <- do.call(paste, c(list(arm), slice, sep = "\r"))
  n <- stats::ave(as.numeric(scored), g, FUN = sum)
  ifelse(scored, 1 / (M * pmax(n, 1)), NA_real_)
}

#' The probability scale of a draw mixture, from each completion's DEPLOYED delta.
#'
#' Each arm rescales its own draws by its own deployed recalibration factor before mixing.
#' An arm without one (47 only WARNS when the delta table fails) would contribute RAW draws to
#' a mixture labelled recalibrated. So: every arm resolves a finite positive delta for every
#' horizon -> "recalibrated-deployed-mi"; no arm carries one -> "raw-mi"; anything in between
#' is a mixture of two probability scales and is refused. The lookup mirrors
#' 21_bayesian_renewal.R::.bayes_recal_delta() (a named vector keyed "<h>" or "h<h>", or one
#' unnamed scalar shared by every horizon).
#' @param deltas list of per-arm deploy deltas (NULL allowed), named by arm tag.
#' @param horizons the horizons the draws carry.
mi_deploy_scale <- function(deltas, horizons) {
  .val <- function(d, h) {
    if (length(d) == 1L && is.null(names(d))) return(suppressWarnings(as.numeric(d)))
    nms <- names(d)
    if (is.null(nms)) return(NA_real_)
    hit <- match(c(as.character(h), paste0("h", h)), nms)
    hit <- hit[!is.na(hit)]
    if (!length(hit)) NA_real_ else suppressWarnings(as.numeric(d[[hit[1]]]))
  }
  absent <- vapply(deltas, function(d) is.null(d) || !length(d), logical(1))
  if (all(absent)) return("raw-mi")
  ok <- vapply(deltas, function(d) {
    if (is.null(d) || !length(d)) return(FALSE)
    all(vapply(horizons, function(h) {
      v <- .val(d, h)
      length(v) == 1L && is.finite(v) && v > 0
    }, logical(1)))
  }, logical(1))
  if (!all(ok)) {
    tags <- names(deltas)
    if (is.null(tags)) tags <- as.character(seq_along(deltas))
    stop(sprintf(paste0("[mi-pool] arm(s) %s carry no usable deployed recalibration factor for ",
                        "every horizon (%s) while other arms do. Their draws would enter the ",
                        "mixture on the RAW scale beside recalibrated ones, under a ",
                        "'recalibrated' label. Re-run those arms (47 only warns when its delta ",
                        "table fails)."),
                 paste(tags[!ok], collapse = ", "), paste(horizons, collapse = ", ")),
         call. = FALSE)
  }
  "recalibrated-deployed-mi"
}

#' Pool a long frame of per-arm rows: one Rubin result per key combination.
#'
#' No magrittr pipe and no bare dplyr verbs: this file is sourced by a test that does not
#' attach dplyr, and a helper that works only when its caller happens to have a package on
#' the search path is a trap rather than a dependency.
#'
#' @param df long frame with one row per (key..., arm).
#' @param keys grouping columns.
#' @param value column holding the per-arm point estimate.
#' @param u optional column holding the per-arm within variance.
mi_pool_by <- function(df, keys, value, u = NULL) {
  stopifnot(is.data.frame(df), all(keys %in% names(df)), value %in% names(df))
  if (!is.null(u) && !u %in% names(df))
    stop("[mi_pool_by] no column '", u, "' in the frame.", call. = FALSE)
  g <- dplyr::group_by(df, dplyr::across(dplyr::all_of(keys)))
  out <- dplyr::group_modify(g, function(d, k) {
    r <- mi_rubin(d[[value]], if (!is.null(u)) d[[u]] else NULL)
    tibble::tibble(q = r$q, B = r$B, U = r$U, T = r$T, se = r$se,
                   riv = r$riv, lambda = r$lambda, fmi = r$fmi, df = r$df,
                   lo = r$lo, hi = r$hi,
                   n_arms = r$M, n_u = r$n_u, u_known = r$u_known)
  })
  dplyr::ungroup(out)
}

#' Apply the "already invaded in THIS completion" rule to one arm's draws.
#'
#' The estimand is "P(this zone records its FIRST confirmed case in the forecast window)",
#' which is genuinely ZERO in a completion where the zone has already been invaded — so the
#' zone belongs in the mixture at zero, not as a missing value. Zeroing the cumulative hazard
#' is exactly that: p = 1 - exp(-0) = 0, and mu_forecast = colMeans(cum) picks up the zero.
#'
#' The rejected alternative — conditioning on at-risk-in-every-completion — silently changes
#' the population the figure describes and is selection on the outcome.
#'
#' @param drl a `.bayes_invasion_draws()` list (named by horizon).
#' @param affected zones already invaded in this completion.
mi_zero_affected <- function(drl, affected) {
  if (!length(affected)) return(drl)
  lapply(drl, function(el) {
    j <- which(el$zones %in% affected)
    if (length(j)) {
      el$cum[, j] <- 0
      el$p[, j]   <- 0
    }
    el
  })
}
