decompose_longstep_matrix <- function(mat_name, mat, 
                                      states = c("FCBud","FCBpre","AL","EDU","ENG","HHS","PSFA","SCI","UNDCL","NE")) {
  # Ensure required states exist (creates 0 where missing)
  missing_states <- setdiff(states, rownames(mat))
  if (length(missing_states)) {
    mat <- rbind(mat, setNames(rep(0, ncol(mat)), colnames(mat)))
    rownames(mat)[nrow(mat)] <- missing_states[1]
    if (length(missing_states) > 1) {
      for (s in missing_states[-1]) {
        mat <- rbind(mat, setNames(rep(0, ncol(mat)), colnames(mat)))
        rownames(mat)[nrow(mat)] <- s
      }
    }
  }
  missing_states <- setdiff(states, colnames(mat))
  if (length(missing_states)) {
    for (s in missing_states) mat <- cbind(mat, setNames(list(rep(0, nrow(mat))), s))
  }
  # Reorder rows/cols
  mat <- mat[states, states, drop = FALSE]
  
  # Parse cohort and k from name like "FT22_1_5" or "TR23_1_1"
  rx <- regexec("^((FT|TR)\\d+)_1_(\\d+)$", mat_name)
  m <- regmatches(mat_name, rx)[[1]]
  if (length(m) == 0) return(NULL)
  cohort <- m[2]; k <- as.integer(m[4])
  
  # Entry (semester 1) composition and term-k composition
  row_tot <- rowSums(mat, na.rm = TRUE)   # entry counts by state
  col_tot <- colSums(mat, na.rm = TRUE)   # term-k counts by state
  
  CohortN_1 <- sum(row_tot, na.rm = TRUE)
  FCB_entry <- row_tot["FCBud"] + row_tot["FCBpre"]
  FCB_termk <- col_tot["FCBud"] + col_tot["FCBpre"]
  
  PctFCB          <- if (CohortN_1 > 0) FCB_entry / CohortN_1 else NA_real_
  PromotionRate   <- if (FCB_termk  > 0) col_tot["FCBud"] / FCB_termk else NA_real_
  
  # Retention over 1→k (for k=1 this yields 1)
  FCB_to_NE <- mat["FCBud","NE"] + mat["FCBpre","NE"]
  Retention <- if (FCB_entry > 0) 1 - (FCB_to_NE / FCB_entry) else NA_real_
  
  # Net internal transfers factor over 1→k:
  # "What fraction of FCB at term k came from *any* entry college (numerator)
  #  relative to what came *only from FCB entry states* (denominator)."
  FCB_from_any <- FCB_termk
  FCB_from_FCBonly <- mat["FCBud","FCBud"] + mat["FCBud","FCBpre"] +
    mat["FCBpre","FCBud"] + mat["FCBpre","FCBpre"]
  TransferFactor <- if (FCB_from_FCBonly > 0) FCB_from_any / FCB_from_FCBonly else NA_real_
  
  # Store in percent to match existing convention
  data.frame(
    Matrix         = mat_name,
    Cohort         = cohort,
    CohortN_1      = CohortN_1,
    PctFCB         = 100 * PctFCB,
    Retention      = 100 * Retention,
    TransferFactor = 100 * TransferFactor,
    PromotionRate  = 100 * PromotionRate,
    FCBud_k        = col_tot["FCBud"],
    stringsAsFactors = FALSE
  )
}

build_fcbud_decomposition_table <- function(matrix_names, transition_matrices,
                                            states = c("FCBud","FCBpre","AL","EDU","ENG","HHS","PSFA","SCI","UNDCL","NE"),
                                            verbose = TRUE) {
  pick <- intersect(matrix_names, names(transition_matrices))
  pick <- grep("^((FT|TR)\\d+)_1_\\d+$", pick, value = TRUE)  # any 1->k
  
  rows <- lapply(pick, function(nm) {
    mat <- transition_matrices[[nm]]
    if (is.null(mat)) return(NULL)
    out <- decompose_longstep_matrix(nm, mat, states = states)
    if (verbose && !is.null(out)) message(sprintf("✓ decomposed %s", nm))
    out
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    warning("No long-step matrices were decomposed.")
    return(data.frame(Matrix=character(), Cohort=character(),
                      CohortN_1=numeric(), PctFCB=numeric(), Retention=numeric(),
                      TransferFactor=numeric(), PromotionRate=numeric(), FCBud_k=numeric(),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Build the decomposition table
all_mats <- names(transition_matrices)              # full set
fcbud_decomp_table <- build_fcbud_decomposition_table(all_mats, transition_matrices)

# View or export
print(fcbud_decomp_table)

save(fcbud_decomp_table,
     file = here::here("_outputs/Rdata", "decomp_table.RData"))
