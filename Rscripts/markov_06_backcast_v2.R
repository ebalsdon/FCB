# --- markov_06_backcast.R ---
# Back-casting F25 and S26 by cohort using only information available through S25.
# Produces both F25 and S26 backcasts for each cohort, saving the full state vectors
# as well as the FCBpre/FCBud summaries for each step.
#
# This version aligns the step used for transition probabilities for each forecast step,
# and supports straightforward multi-step Markov projections.

library(dplyr)
library(tibble)
library(glue)
library(gt)

# Assumes these are loaded in your R session or loaded from files:
# - fcb_pre_df
# - fcb_ud_df
# - transition_matrices
# - transition_probabilities
# - term_starts
# - terms

backcast_F25_S26_from_S25 <- function(
    fcb_pre_df, fcb_ud_df,
    transition_matrices, transition_probabilities,
    term_starts, terms,
    n_reference = 2
) {
  get_states <- function() {
    if (length(transition_probabilities)) {
      return(colnames(transition_probabilities[[1]]))
    }
    if (length(transition_matrices)) {
      return(colnames(transition_matrices[[1]]))
    }
    c("FCBud","FCBpre","AL","EDU","ENG","HHS","PSFA","SCI","UNDCL","NE")
  }
  
  # Find the transition matrix key where S25 is the "to" term for a cohort
  find_s25_to_key <- function(cohort) {
    # S25 index in terms
    s25_idx <- match("S25", terms)
    cohort_start_idx <- match(term_starts[[cohort]], terms)
    if (is.na(cohort_start_idx)) return(NULL)
    step_to_sem <- (s25_idx - cohort_start_idx) + 1
    # Look for key with to_sem == step_to_sem, and cohort name
    possible_keys <- grep(paste0("^", cohort, "_\\d+_", step_to_sem, "$"), names(transition_matrices), value = TRUE)
    if (length(possible_keys) == 0) return(NULL)
    possible_keys[[1]]
  }
  
  # Extract from_sem, to_sem from a key like FT22_3_4
  extract_step_from_key <- function(key) {
    # returns list(from=3, to=4) for FT22_3_4
    parts <- unlist(strsplit(key, "_"))
    list(from = as.integer(parts[2]), to = as.integer(parts[3]))
  }
  
  # For a given step (from_sem -> to_sem), get transition probabilities from n_reference previous cohorts
  avg_ref_P <- function(group, from_sem, to_sem, n_ref, states) {
    buckets <- list(); k <- 1
    while (length(buckets) < n_ref) {
      prev_y <- as.integer(group$year) - k
      cand <- paste0(group$prefix, sprintf("%02d", prev_y))
      key <- paste0(cand, "_", from_sem, "_", to_sem)
      if (key %in% names(transition_probabilities)) {
        P <- transition_probabilities[[key]]
        if (max(P, na.rm=TRUE) > 1.0001) P <- P / 100
        P <- P[states, states, drop=FALSE]; P[is.na(P)] <- 0
        buckets[[length(buckets) + 1]] <- P
      }
      k <- k + 1; if (k > 40) break
    }
    if (!length(buckets)) return(NULL)
    Reduce(`+`, buckets) / length(buckets)
  }
  
  # --- Main logic ---
  states <- get_states()
  cohorts <- setdiff(rownames(fcb_ud_df), "Total")
  
  # Set up results tibble with all forecasted states for F25 and S26
  results <- tibble(Cohort = cohorts)
  for (state in states) {
    results[[paste0("F25_", state, "_forecast")]] <- NA_real_
    results[[paste0("S26_", state, "_forecast")]] <- NA_real_
  }
  results$F25_pre_forecast <- NA_real_
  results$F25_ud_forecast  <- NA_real_
  results$S26_pre_forecast <- NA_real_
  results$S26_ud_forecast  <- NA_real_
  
  for (cohort in cohorts) {
    # 1. Find S25 "to" transition matrix for this cohort
    key <- find_s25_to_key(cohort)
    if (is.null(key)) next
    step <- extract_step_from_key(key)
    from_sem <- step$from
    to_sem <- step$to
    
    # 2. Extract S25 state vector (column totals)
    xt <- colSums(transition_matrices[[key]], na.rm=TRUE)
    xt <- xt[states]; xt[is.na(xt)] <- 0
    xt <- matrix(xt, nrow = 1, dimnames = list(NULL, states))
    
    # 3. Average corresponding *_X_Y transition matrices from previous cohorts for F25
    group <- list(
      prefix = substr(cohort, 1, 2),
      year   = as.integer(substr(cohort, 3, 4))
    )
    # F25 step: use next step after S25
    P_F25 <- avg_ref_P(group, from_sem + 1, to_sem + 1, n_reference, states)
    if (is.null(P_F25)) next
    
    # 4. Forecast full F25 state vector
    x_F25 <- xt %*% P_F25
    
    # Store all F25 state values
    for (state in states) {
      results[[paste0("F25_", state, "_forecast")]][results$Cohort == cohort] <- x_F25[1, state]
    }
    results$F25_pre_forecast[results$Cohort == cohort] <- x_F25[1, "FCBpre"]
    results$F25_ud_forecast [results$Cohort == cohort] <- x_F25[1, "FCBud"]
    
    # --- S26 backcast ---
    # Use F25 forecasted state vector as input
    # For S26, use next transition step: from_sem+2 -> to_sem+2
    P_S26 <- avg_ref_P(group, from_sem + 2, to_sem + 2, n_reference, states)
    if (is.null(P_S26)) next
    
    x_S26 <- x_F25 %*% P_S26
    
    # Store all S26 state values
    for (state in states) {
      results[[paste0("S26_", state, "_forecast")]][results$Cohort == cohort] <- x_S26[1, state]
    }
    results$S26_pre_forecast[results$Cohort == cohort] <- x_S26[1, "FCBpre"]
    results$S26_ud_forecast [results$Cohort == cohort] <- x_S26[1, "FCBud"]
  }
  
  results
}

# ---- Run backcasting ----
backcast_results <- backcast_F25_S26_from_S25(
  fcb_pre_df, fcb_ud_df,
  transition_matrices, transition_probabilities,
  term_starts, terms,
  n_reference = 2
)

# ---- Gather actuals for comparison ----
get_actual_safe <- function(df, cohort, term) {
  if (cohort %in% rownames(df) && term %in% colnames(df)) {
    val <- df[cohort, term]
    if (is.na(val)) NA_real_ else as.numeric(val)
  } else {
    NA_real_
  }
}

backcast_results <- backcast_results %>%
  rowwise() %>%
  mutate(
    F25_pre_actual = get_actual_safe(fcb_pre_df, Cohort, "F25"),
    F25_ud_actual  = get_actual_safe(fcb_ud_df,  Cohort, "F25"),
    S26_pre_actual = get_actual_safe(fcb_pre_df, Cohort, "S26"),
    S26_ud_actual  = get_actual_safe(fcb_ud_df,  Cohort, "S26")
  ) %>%
  ungroup() %>%
  mutate(
    across(
      c(
        F25_pre_actual, F25_ud_actual, S26_pre_actual, S26_ud_actual,
        F25_pre_forecast, F25_ud_forecast, S26_pre_forecast, S26_ud_forecast
      ),
      ~ifelse(Cohort %in% c("FT17", "TR19"), 0, .x)
    )
  ) %>%
  mutate(
    F25_pre_error  = F25_pre_forecast - F25_pre_actual,
    F25_ud_error   = F25_ud_forecast  - F25_ud_actual,
    S26_pre_error  = S26_pre_forecast - S26_pre_actual,
    S26_ud_error   = S26_ud_forecast  - S26_ud_actual
  )

# ---- Formatted comparison table ----
backcast_table <- backcast_results %>%
  mutate(
    across(
      starts_with("F25_") | starts_with("S26_"),
      ~ifelse(is.na(.x), NA, round(.x, 0))
    ),
    F25_pre_error = sprintf("%+d", F25_pre_error),
    F25_ud_error  = sprintf("%+d", F25_ud_error),
    S26_pre_error = sprintf("%+d", S26_pre_error),
    S26_ud_error  = sprintf("%+d", S26_ud_error)
  ) %>%
  select(
    Cohort,
    F25_pre_actual, F25_pre_forecast, F25_pre_error,
    F25_ud_actual,  F25_ud_forecast,  F25_ud_error,
    S26_pre_actual, S26_pre_forecast, S26_pre_error,
    S26_ud_actual,  S26_ud_forecast,  S26_ud_error
  )

# ---- Save results (optional) ----
save(backcast_results, backcast_table, file = here::here("_outputs/Rdata", "fcb_backcast_F25_S26_S25.RData"))
save.image(here::here("_outputs/Rdata/data_image1.RData"))
