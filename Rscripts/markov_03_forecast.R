forecast_with_averaging <- function(
    fcb_pre_df, fcb_ud_df,
    transition_matrices, transition_probabilities,
    term_starts, terms,
    n_reference = 1,
    new_cohort_proxy = c("skip", "prev_cohort")  # proxy S26 for FT25/TR25?
) {
  new_cohort_proxy <- match.arg(new_cohort_proxy)
  
  # ---- helpers ----
  ensure_col <- function(df, col) {
    if (!(col %in% colnames(df))) {
      df[, col] <- NA_real_
      df <- df[, unique(c(setdiff(colnames(df), col), col)), drop = FALSE]
    }
    df
  }
  
  get_states <- function() {
    if (length(transition_probabilities)) {
      return(colnames(transition_probabilities[[1]]))
    }
    if (length(transition_matrices)) {
      return(colnames(transition_matrices[[1]]))
    }
    # last-resort: canonical order you’ve been using
    c("FCBud","FCBpre","AL","EDU","ENG","HHS","PSFA","SCI","UNDCL","NE")
  }
  
  prev_cohort <- function(cohort, k = 1) {
    g <- substr(cohort, 1, 2); y <- as.integer(substr(cohort, 3, 4))
    paste0(g, sprintf("%02d", y - k))
  }
  
  sem_at_term <- function(cohort, term_label) {
    if (!cohort %in% names(term_starts)) return(NA_integer_)
    s0 <- match(term_starts[[cohort]], terms)
    t  <- match(term_label, terms)
    if (is.na(s0) || is.na(t)) return(NA_integer_)
    (t - s0) + 1
  }
  
  avg_ref_P <- function(cohort, to_sem, n_ref, states) {
    group <- substr(cohort, 1, 2)
    buckets <- list(); k <- 1
    while (length(buckets) < n_ref) {
      cand <- prev_cohort(cohort, k)
      if (substr(cand,1,2) != group) break
      key <- paste0(cand, "_", to_sem, "_", to_sem + 1)
      if (key %in% names(transition_probabilities)) {
        P <- transition_probabilities[[key]]
        if (max(P, na.rm = TRUE) > 1.0001) P <- P / 100
        # align to full state order and coerce NAs to 0
        P <- P[states, states, drop = FALSE]
        P[is.na(P)] <- 0
        buckets[[length(buckets) + 1]] <- P
      } else {
        message(sprintf("   • avg_ref_P: missing prob matrix %s", key))
      }
      k <- k + 1; if (k > 40) break
    }
    if (!length(buckets)) return(NULL)
    Reduce(`+`, buckets) / length(buckets)
  }
  
  # base S25 vector uses (to_sem-1 -> to_sem) column sums
  x_from_columnsums <- function(cohort, to_sem, states) {
    from_sem <- to_sem - 1
    if (from_sem < 1) return(NULL)
    key <- paste0(cohort, "_", from_sem, "_", to_sem)
    if (!(key %in% names(transition_matrices))) {
      message(sprintf("   • x_t: missing counts matrix %s", key))
      return(NULL)
    }
    xt <- colSums(transition_matrices[[key]], na.rm = TRUE)
    xt <- xt[states]
    xt[is.na(xt)] <- 0
    matrix(xt, nrow = 1, dimnames = list(NULL, states))
  }
  
  write_if_empty <- function(df, cohort, term, value) {
    if (cohort %in% rownames(df) && term %in% colnames(df)) {
      if (is.na(df[cohort, term])) df[cohort, term] <- as.numeric(value)
    }
    df
  }
  
  # ---- prep ----
  fcb_pre_df <- ensure_col(fcb_pre_df, "F25")
  fcb_pre_df <- ensure_col(fcb_pre_df, "S26")
  fcb_ud_df  <- ensure_col(fcb_ud_df,  "F25")
  fcb_ud_df  <- ensure_col(fcb_ud_df,  "S26")
  
  # normalize types
  fcb_pre_df[, c("F25","S26")] <- lapply(fcb_pre_df[, c("F25","S26")], as.numeric)
  fcb_ud_df[,  c("F25","S26")] <- lapply(fcb_ud_df[,  c("F25","S26")], as.numeric)
  
  states <- get_states()
  stopifnot(all(c("FCBud","FCBpre") %in% states))
  
  cohorts <- setdiff(rownames(fcb_ud_df), "Total")
  
  # ---- (1) main: cohorts with S25 base → F25, S26 ----
  for (cohort in cohorts) {
    to_sem <- sem_at_term(cohort, "S25")
    if (is.na(to_sem)) next
    
    xt <- x_from_columnsums(cohort, to_sem, states)
    if (is.null(xt)) next
    
    P  <- avg_ref_P(cohort, to_sem, n_reference, states)  # S25->F25
    if (is.null(P)) {
      message(sprintf("   • no ref P for %s step %s->%s", cohort, to_sem, to_sem+1))
      next
    }
    
    x1 <- xt %*% P         # F25
    x2 <- x1 %*% P         # S26
    
    fcb_ud_df  <- write_if_empty(fcb_ud_df,  cohort, "F25", round(x1[1, "FCBud"]))
    fcb_pre_df <- write_if_empty(fcb_pre_df, cohort, "F25", round(x1[1, "FCBpre"]))
    fcb_ud_df  <- write_if_empty(fcb_ud_df,  cohort, "S26", round(x2[1, "FCBud"]))
    fcb_pre_df <- write_if_empty(fcb_pre_df, cohort, "S26", round(x2[1, "FCBpre"]))
    
    message(sprintf("✅ Forecasted %s: base=S25 (to_sem=%d), wrote F25/S26 if empty", cohort, to_sem))
  }
  
  # ---- (2) force zeros for cohorts beyond program length ----
  for (cohort in cohorts) {
    sF25 <- sem_at_term(cohort, "F25")
    sS26 <- sem_at_term(cohort, "S26")
    group <- substr(cohort, 1, 2)
    max_sem <- if (group == "FT") 16L else 12L
    
    if (!is.na(sF25) && sF25 > max_sem) {
      if (is.na(fcb_ud_df[cohort,"F25"]))  fcb_ud_df[cohort,"F25"]  <- 0
      if (is.na(fcb_pre_df[cohort,"F25"])) fcb_pre_df[cohort,"F25"] <- 0
      message(sprintf("↳ Set %s F25 = 0 (beyond max semesters)", cohort))
    }
    if (!is.na(sS26) && sS26 > max_sem) {
      if (is.na(fcb_ud_df[cohort,"S26"]))  fcb_ud_df[cohort,"S26"]  <- 0
      if (is.na(fcb_pre_df[cohort,"S26"])) fcb_pre_df[cohort,"S26"] <- 0
      message(sprintf("↳ Set %s S26 = 0 (beyond max semesters)", cohort))
    }
  }
  
# ---- (3) proxy S26 for FT25/TR25 using prev cohort S25 -> F25 (1_2 then 2->3) ----
  if (new_cohort_proxy == "prev_cohort") {
    for (cohort in c("FT25", "TR25")) {
      if (!(cohort %in% rownames(fcb_ud_df))) next
      
      # Require manual F25 numbers present for the new cohort (what we will override with)
      if (is.na(fcb_ud_df[cohort, "F25"]) || is.na(fcb_pre_df[cohort, "F25"])) {
        message(sprintf("   • proxy: %s skipped (missing manual F25 FCBud/FCBpre)", cohort))
        next
      }
      
      # Averaged P for new cohort F25->S26 (1->2 step for FT25/TR25)
      to_sem_F25 <- sem_at_term(cohort, "F25")  # = 1
      if (is.na(to_sem_F25)) {
        message(sprintf("   • proxy: %s skipped (no sem_at_term for F25)", cohort)); next
      }
      P_1_2 <- avg_ref_P(cohort, to_sem_F25, n_reference, states)  # (1->2)
      if (is.null(P_1_2)) {
        message(sprintf("   • proxy: %s skipped (no prior P for 1->2)", cohort)); next
      }
      
      # Build PREVIOUS cohort's F25 composition via its S25 and (2->3) transition
      prev <- prev_cohort(cohort, 1)       # FT24 or TR24
      k_prev_F25 <- sem_at_term(prev, "F25")   # = 3 for FT24/TR24
      if (is.na(k_prev_F25)) {
        message(sprintf("   • proxy: %s skipped (prev cohort %s lacks F25 anchor)", cohort, prev))
        next
      }
      
      # Get prev S25 composition from counts matrix (1_2): FT24_1_2 / TR24_1_2
      key_prev_1_2 <- paste0(prev, "_", k_prev_F25-2, "_", k_prev_F25-1)  # 1_2
      if (!(key_prev_1_2 %in% names(transition_matrices))) {
        message(sprintf("   • proxy: %s skipped (missing %s)", cohort, key_prev_1_2))
        next
      }
      x_prev_S25 <- colSums(transition_matrices[[key_prev_1_2]], na.rm = TRUE)
      x_prev_S25 <- x_prev_S25[states]; x_prev_S25[is.na(x_prev_S25)] <- 0
      x_prev_S25 <- matrix(x_prev_S25, nrow = 1, dimnames = list(NULL, states))
      
      # Need averaged probabilities for (2->3) from prior cohorts (same group as prev)
      P_2_3 <- avg_ref_P(prev, k_prev_F25-1, n_reference, states)  # (2->3)
      if (is.null(P_2_3)) {
        message(sprintf("   • proxy: %s skipped (no prior P for prev %s 2->3)", cohort, prev))
        next
      }
      
      # Derive prev F25 composition: S25 %*% P(2->3)
      comp_prev_F25 <- x_prev_S25 %*% P_2_3
      comp_prev_F25 <- as.numeric(comp_prev_F25)
      names(comp_prev_F25) <- states
      
      # Override with manual new-cohort F25 FCB values
      comp_prev_F25["FCBud"]  <- fcb_ud_df [cohort, "F25"]
      comp_prev_F25["FCBpre"] <- fcb_pre_df[cohort, "F25"]
      xt_F25 <- matrix(comp_prev_F25, nrow = 1, dimnames = list(NULL, states))
      
      # One step to S26 using P(1->2) for the NEW cohort
      xS26 <- xt_F25 %*% P_1_2
      
      # Write only if empty
      if (is.na(fcb_ud_df[cohort, "S26"]))  fcb_ud_df [cohort, "S26"] <- round(xS26[1, "FCBud"])
      if (is.na(fcb_pre_df[cohort, "S26"])) fcb_pre_df[cohort, "S26"] <- round(xS26[1, "FCBpre"])
      
      message(sprintf("✅ Proxy S26 for %s via %s: used %s (S25) + P(2->3); then P(1->2). F25 FCBud=%d, FCBpre=%d",
                      cohort, prev, key_prev_1_2,
                      as.integer(xt_F25[1,"FCBud"]), as.integer(xt_F25[1,"FCBpre"])))
    }
  }
  
  
  list(fcb_pre_df = fcb_pre_df, fcb_ud_df = fcb_ud_df)
}


## 1) Ensure F25/S26 columns exist and are numeric
for (df in c("fcb_ud_df","fcb_pre_df")) {
  assign(df, within(get(df), {
    if (!"F25" %in% names(get(df))) F25 <- NA_real_
    if (!"S26" %in% names(get(df))) S26 <- NA_real_
  }))
}
fcb_ud_df[,"F25"]  <- as.numeric(fcb_ud_df[,"F25"])
fcb_pre_df[,"F25"] <- as.numeric(fcb_pre_df[,"F25"])

## 2) Set your manual F25 values (example numbers — replace with your latest)
# FCB upper-division (ud) at F25:
fcb_ud_df["FT25","F25"] <- 0      # if you expect 0 UD at entry; else your value
fcb_ud_df["TR25","F25"] <- 850    # put your current best estimate

# FCB pre-major (pre) at F25:
fcb_pre_df["FT25","F25"] <- 1110  # your provisional number
fcb_pre_df["TR25","F25"] <- 200   # your provisional number

## 3) Run the forecast WITH proxy for new cohorts
res <- forecast_with_averaging(
  fcb_pre_df, fcb_ud_df,
  transition_matrices, transition_probabilities,
  term_starts, terms,
  n_reference = 2,
  new_cohort_proxy = "prev_cohort"
)

fcb_pre_df <- res$fcb_pre_df
fcb_ud_df  <- res$fcb_ud_df

recompute_totals <- function(df, sum_from = "F22", key_col = NULL) {
  # 1) Decide how to identify rows (rownames vs a Cohort column)
  has_key_col <- !is.null(key_col) && key_col %in% colnames(df)
  if (has_key_col) {
    cohort_vals <- df[[key_col]]
    total_idx <- which(cohort_vals == "Total")
    data_rows <- which(cohort_vals != "Total" & !is.na(cohort_vals))
  } else {
    rn <- rownames(df)
    if (is.null(rn)) stop("No rownames and no key_col provided.")
    total_idx <- which(rn == "Total")
    data_rows <- which(rn != "Total")
  }
  
  # 2) Ensure sum_from exists and get the target column names
  start <- match(sum_from, colnames(df))
  if (is.na(start)) stop(sprintf("sum_from '%s' not found in df.", sum_from))
  cols <- colnames(df)[start:ncol(df)]
  
  # 3) If Total row doesn’t exist, create it (NA) with same class as df
  if (length(total_idx) == 0) {
    # Create a one-row “Total” with NAs
    tot_row <- df[1, , drop = FALSE]; tot_row[] <- NA
    if (has_key_col) tot_row[[key_col]] <- "Total"
    row_to_bind <- if (has_key_col) tot_row else {
      rownames(tot_row) <- "Total"; tot_row
    }
    df <- rbind(df, row_to_bind)
    # recompute indices
    if (has_key_col) {
      total_idx <- which(df[[key_col]] == "Total")
      data_rows <- which(df[[key_col]] != "Total" & !is.na(df[[key_col]]))
    } else {
      total_idx <- which(rownames(df) == "Total")
      data_rows <- which(rownames(df) != "Total")
    }
  }
  
  # 4) Sum only numeric columns among the target cols
  numeric_cols <- cols[sapply(df[cols], is.numeric)]
  if (length(numeric_cols) == 0) stop("No numeric columns to sum.")
  totals <- colSums(df[data_rows, numeric_cols, drop = FALSE], na.rm = TRUE)
  
  # 5) Write back by column names (prevents positional mismatches)
  df[total_idx, names(totals)] <- totals
  
  df
}


fcb_ud_df  <- recompute_totals(fcb_ud_df,  sum_from = "F22")
fcb_pre_df <- recompute_totals(fcb_pre_df, sum_from = "F22")


# ---- save with forecasts ----
save(fcb_pre_df, fcb_ud_df,
     file = here::here("_outputs/Rdata", "cohorts_with_forecast_dfs.RData"))




