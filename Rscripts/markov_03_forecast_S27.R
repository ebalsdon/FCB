library(dplyr)
library(tibble)

ensure_col <- function(df, col) {
  if (!(col %in% colnames(df))) {
    df[, col] <- NA_real_
  }
  df
}

get_states <- function(transition_probabilities, transition_matrices) {
  if (length(transition_probabilities)) {
    return(colnames(transition_probabilities[[1]]))
  }
  if (length(transition_matrices)) {
    return(colnames(transition_matrices[[1]]))
  }
  stop("No transition matrices are available.")
}

sem_at_term <- function(cohort, term_label, term_starts, terms) {
  if (!cohort %in% names(term_starts)) return(NA_integer_)
  start_index <- match(term_starts[[cohort]], terms)
  target_index <- match(term_label, terms)
  if (is.na(start_index) || is.na(target_index)) return(NA_integer_)
  target_index - start_index + 1L
}

avg_ref_P <- function(cohort, to_sem, n_reference, transition_probabilities, states) {
  group <- substr(cohort, 1, 2)
  year <- as.integer(substr(cohort, 3, 4))
  buckets <- list()

  for (offset in seq_len(40L)) {
    candidate <- paste0(group, sprintf("%02d", year - offset))
    key <- paste0(candidate, "_", to_sem, "_", to_sem + 1L)
    if (!(key %in% names(transition_probabilities))) next

    P <- transition_probabilities[[key]]
    if (max(P, na.rm = TRUE) > 1.0001) P <- P / 100
    P <- P[states, states, drop = FALSE]
    P[is.na(P)] <- 0
    buckets[[length(buckets) + 1L]] <- P
    if (length(buckets) == n_reference) break
  }

  if (!length(buckets)) return(NULL)
  average <- Reduce(`+`, buckets) / length(buckets)
  attr(average, "reference_count") <- length(buckets)
  average
}

# Regime adjustment for the transfer-into-FCBpre policy (30 units + GPA gate, first applied F25>S26).
# If the most recent same-stage reference transition started in or after `post_from`, its
# college->FCBpre column replaces the n_reference average; non-switchers stay on the college diagonal
# (row sums preserved). Optional `uplift` raises P(FCBpre->FCBud) for stage 3>4 (scenario knob, default 0).
apply_regime <- function(P, cohort, to_sem, transition_probabilities, states,
                         term_starts, terms, post_from, uplift = 0) {
  latest <- avg_ref_P(cohort, to_sem, 1L, transition_probabilities, states)
  group <- substr(cohort, 1, 2)
  year <- as.integer(substr(cohort, 3, 4))
  ref_cohort <- NULL
  for (offset in seq_len(40L)) {
    candidate <- paste0(group, sprintf("%02d", year - offset))
    if (paste0(candidate, "_", to_sem, "_", to_sem + 1L) %in% names(transition_probabilities)) {
      ref_cohort <- candidate
      break
    }
  }
  applied <- FALSE
  if (!is.null(ref_cohort)) {
    ref_from_term <- terms[match(term_starts[[ref_cohort]], terms) + to_sem - 1L]
    if (!is.na(ref_from_term) && match(ref_from_term, terms) >= match(post_from, terms)) {
      colleges <- setdiff(states, c("FCBud", "FCBpre", "NE"))
      for (i in colleges) {
        delta <- P[i, "FCBpre"] - latest[i, "FCBpre"]
        P[i, "FCBpre"] <- latest[i, "FCBpre"]
        P[i, i] <- P[i, i] + delta
      }
      applied <- TRUE
    }
  }
  if (uplift != 0 && to_sem == 3L) {
    P["FCBpre", "FCBud"] <- P["FCBpre", "FCBud"] + uplift
    P["FCBpre", "FCBpre"] <- P["FCBpre", "FCBpre"] - uplift
  }
  attr(P, "reference_count") <- NULL
  attr(P, "regime_applied") <- applied
  P
}

f26_state_vector <- function(cohort, to_sem, transition_matrices, states) {
  key <- if (to_sem == 1L) {
    paste0(cohort, "_1_1")
  } else {
    paste0(cohort, "_", to_sem - 1L, "_", to_sem)
  }
  if (!(key %in% names(transition_matrices))) return(NULL)

  counts <- colSums(transition_matrices[[key]], na.rm = TRUE)
  counts <- counts[states]
  counts[is.na(counts)] <- 0
  matrix(counts, nrow = 1L, dimnames = list(NULL, states))
}

recompute_totals <- function(df, sum_from = "F22") {
  start <- match(sum_from, colnames(df))
  if (is.na(start)) stop(sprintf("sum_from '%s' not found in df.", sum_from))
  numeric_cols <- colnames(df)[start:ncol(df)]
  numeric_cols <- numeric_cols[vapply(df[numeric_cols], is.numeric, logical(1))]
  if (!length(numeric_cols)) stop("No numeric columns to sum.")

  if (!("Total" %in% rownames(df))) {
    total_row <- df[1, , drop = FALSE]
    total_row[] <- NA
    rownames(total_row) <- "Total"
    df <- rbind(df, total_row)
  }

  data_rows <- which(rownames(df) != "Total")
  df["Total", numeric_cols] <- colSums(df[data_rows, numeric_cols, drop = FALSE], na.rm = TRUE)
  df
}

forecast_s27 <- function(
    fcb_pre_df, fcb_ud_df,
    transition_matrices, transition_probabilities,
    term_starts, terms,
    n_reference = 3L,
    regime_aware = FALSE,
    regime_post_from = "F25",
    regime_uplift = 0
) {
  states <- get_states(transition_probabilities, transition_matrices)
  if (!all(c("FCBud", "FCBpre") %in% states)) {
    stop("Transition matrices must include FCBud and FCBpre states.")
  }

  fcb_pre_df <- ensure_col(fcb_pre_df, "S27")
  fcb_ud_df <- ensure_col(fcb_ud_df, "S27")
  fcb_pre_df$S27 <- as.numeric(fcb_pre_df$S27)
  fcb_ud_df$S27 <- as.numeric(fcb_ud_df$S27)

  cohorts <- setdiff(rownames(fcb_ud_df), "Total")
  for (cohort in cohorts) {
    to_sem <- sem_at_term(cohort, "F26", term_starts, terms)
    s27_sem <- sem_at_term(cohort, "S27", term_starts, terms)
    if (is.na(to_sem) || is.na(s27_sem)) next

    max_sem <- if (startsWith(cohort, "FT")) 16L else 12L
    if (s27_sem > max_sem) {
      if (is.na(fcb_pre_df[cohort, "S27"])) fcb_pre_df[cohort, "S27"] <- 0
      if (is.na(fcb_ud_df[cohort, "S27"])) fcb_ud_df[cohort, "S27"] <- 0
      next
    }

    xt <- f26_state_vector(cohort, to_sem, transition_matrices, states)
    if (is.null(xt)) {
      message(sprintf("Skipping %s S27 forecast: no F26 actual matrix.", cohort))
      next
    }

    P <- avg_ref_P(cohort, to_sem, n_reference, transition_probabilities, states)
    if (is.null(P)) {
      message(sprintf("Skipping %s S27 forecast: no reference probability matrices for step %d->%d.",
                      cohort, to_sem, to_sem + 1L))
      next
    }

    ref_count <- attr(P, "reference_count")
    if (isTRUE(regime_aware)) {
      P <- apply_regime(P, cohort, to_sem, transition_probabilities, states,
                        term_starts, terms, regime_post_from, regime_uplift)
    }

    x_s27 <- xt %*% P
    if (is.na(fcb_ud_df[cohort, "S27"])) {
      fcb_ud_df[cohort, "S27"] <- round(x_s27[1, "FCBud"])
    }
    if (is.na(fcb_pre_df[cohort, "S27"])) {
      fcb_pre_df[cohort, "S27"] <- round(x_s27[1, "FCBpre"])
    }
    message(sprintf("Forecasted %s S27 from actual F26 using %d reference cohort(s)%s.",
                    cohort, ref_count,
                    if (isTRUE(attr(P, "regime_applied"))) " + post-policy transfer adjustment" else ""))
  }

  list(
    fcb_pre_df = recompute_totals(fcb_pre_df),
    fcb_ud_df = recompute_totals(fcb_ud_df)
  )
}

# Forecast settings: defined in source_code.R; defaults here keep the script runnable on its own.
if (!exists("forecast_n_reference")) forecast_n_reference <- 3L
if (!exists("forecast_regime_aware")) forecast_regime_aware <- TRUE
if (!exists("forecast_regime_post_from")) forecast_regime_post_from <- "F25"
if (!exists("forecast_regime_uplift")) forecast_regime_uplift <- 0

res <- forecast_s27(
  fcb_pre_df, fcb_ud_df,
  transition_matrices, transition_probabilities,
  term_starts, terms,
  n_reference = forecast_n_reference,
  regime_aware = forecast_regime_aware,
  regime_post_from = forecast_regime_post_from,
  regime_uplift = forecast_regime_uplift
)

fcb_pre_df <- res$fcb_pre_df
fcb_ud_df <- res$fcb_ud_df

save(fcb_pre_df, fcb_ud_df,
     file = here::here("_outputs/Rdata", "cohorts_with_S27_forecast.RData"))
