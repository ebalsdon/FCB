# --- Waterfall decomposition table ---

build_fcbud_waterfall <- function(fcb_ud_df,
                                  fcbud_decomp_table,
                                  term1 = "F23",
                                  term2 = "F24",
                                  term_starts = NULL,
                                  terms = colnames(fcb_ud_df),
                                  include_singletons = TRUE,
                                  verbose = FALSE) {
  stopifnot(term1 %in% terms, term2 %in% terms)
  
  # Helper: entry term for a cohort (or pass map)
  if (is.null(term_starts)) {
    term_start <- function(cohort) paste0("F", sub("^(FT|TR)", "", cohort))
  } else {
    term_start <- function(cohort) unname(term_starts[[cohort]])
  }
  
  parse_id <- function(x) as.integer(sub("^(FT|TR)", "", x))
  newest_first <- function(v) v[order(parse_id(v), decreasing = TRUE)]
  
  all_cohorts <- rownames(fcb_ud_df)
  ft_list <- newest_first(all_cohorts[startsWith(all_cohorts, "FT")])
  tr_list <- newest_first(all_cohorts[startsWith(all_cohorts, "TR")])
  
  make_pairs <- function(v) if (length(v) >= 2) Map(\(a,b) c(prev=b, curr=a), v[-length(v)], v[-1]) else list()
  pairs <- c(make_pairs(ft_list), make_pairs(tr_list))
  
  rows <- list(); r <- 0L
  idx_t1 <- match(term1, terms); idx_t2 <- match(term2, terms)
  
  for (p in pairs) {
    prev <- p["prev"]; curr <- p["curr"]
    if (substr(prev,1,2) != substr(curr,1,2)) next
    
    s_prev <- idx_t1 - match(term_start(prev), terms) + 1
    s_curr <- idx_t2 - match(term_start(curr), terms) + 1
    if (!is.finite(s_prev) || !is.finite(s_curr) || s_prev != s_curr) next
    sem <- as.integer(s_curr)
    
    y1 <- suppressWarnings(as.numeric(fcb_ud_df[prev, term1]))
    y2 <- suppressWarnings(as.numeric(fcb_ud_df[curr, term2]))
    if (is.na(y1) || is.na(y2)) next
    delta <- y2 - y1
    
    A <- B <- C <- D <- E <- NA_real_
    mat_prev <- sprintf("%s_1_%d", prev, sem)
    mat_curr <- sprintf("%s_1_%d", curr, sem)
    d1 <- fcbud_decomp_table[fcbud_decomp_table$Matrix == mat_prev, , drop = FALSE]
    d2 <- fcbud_decomp_table[fcbud_decomp_table$Matrix == mat_curr, , drop = FALSE]
    
    if (nrow(d1) == 1L && nrow(d2) == 1L && y1 > 0 && y2 > 0) {
      if (verbose) message(sprintf("✔ Decomp using %s vs %s", mat_prev, mat_curr))
      ldiff <- log(y2 / y1)
      if (is.finite(ldiff) && abs(ldiff) > .Machine$double.eps) {
        logs <- log(c(
          as.numeric(d2$CohortN_1)      / as.numeric(d1$CohortN_1),
          as.numeric(d2$PctFCB)         / as.numeric(d1$PctFCB),
          as.numeric(d2$Retention)      / as.numeric(d1$Retention),
          as.numeric(d2$TransferFactor) / as.numeric(d1$TransferFactor),
          as.numeric(d2$PromotionRate)  / as.numeric(d1$PromotionRate)
        ))
        contrib <- logs / ldiff * delta
        A <- contrib[1]; B <- contrib[2]; C <- contrib[3]; D <- contrib[4]; E <- contrib[5]
      }
    } else if (verbose) {
      message(sprintf("… No decomp for %s/%s (missing %s or %s)", prev, curr, mat_prev, mat_curr))
    }
    
    r <- r + 1L
    row <- list(
      Cohort = sprintf("%s/%s", prev, curr),
      Semester = as.integer(sem),
      Diff = delta,
      Pct_of_Total = NA_real_,
      A = round(A, 1), B = round(B, 1), C = round(C, 1), D = round(D, 1), E = round(E, 1)
    )
    row[[term1]] <- y1
    row[[term2]] <- y2
    rows[[r]] <- as.data.frame(row[c("Cohort","Semester", term1, term2, "Diff","Pct_of_Total","A","B","C","D","E")],
                               stringsAsFactors = FALSE)
  }
  
  # Bind early rows (may be empty)
  out <- if (r > 0L) do.call(rbind, rows) else {
    data.frame(Cohort=character(), Semester=integer(),
               !!term1 := numeric(), !!term2 := numeric(),
               Diff=numeric(), Pct_of_Total=numeric(),
               A=numeric(), B=numeric(), C=numeric(), D=numeric(), E=numeric(),
               check.names = FALSE)
  }
  
  # --- Add singleton rows for unmatched oldest cohorts (optional) ---
  if (include_singletons) {
    add_singleton <- function(cohort_vec, label_prefix) {
      if (length(cohort_vec) == 0L) return(NULL)
      oldest <- tail(cohort_vec, 1)               # smallest ID
      y2 <- suppressWarnings(as.numeric(fcb_ud_df[oldest, term2]))
      if (is.na(y2) || y2 <= 0) return(NULL)
      
      # semester for the oldest cohort at term2
      sem <- idx_t2 - match(term_start(oldest), terms) + 1
      if (!is.finite(sem)) return(NULL)
      
      srow <- list(
        Cohort = oldest,                          # single label (not a pair)
        Semester = as.integer(sem),
        Diff = y2,                                # treat prior term as 0
        Pct_of_Total = NA_real_,
        A = NA, B = NA, C = NA, D = NA, E = NA
      )
      srow[[term1]] <- NA_real_                   # display as NA per your spec
      srow[[term2]] <- y2
      as.data.frame(srow[c("Cohort","Semester", term1, term2, "Diff","Pct_of_Total","A","B","C","D","E")],
                    stringsAsFactors = FALSE)
    }
    
    # Oldest FT (e.g., FT17) and oldest TR (e.g., TR19)
    ft_single <- add_singleton(ft_list, "FT")
    tr_single <- add_singleton(tr_list, "TR")
    adders <- do.call(rbind, Filter(Negate(is.null), list(ft_single, tr_single)))
    if (!is.null(adders)) out <- rbind(out, adders)
  }
  
  # Sort: FT block (newest→oldest), then TR block (newest→oldest)
  out$Type <- substr(out$Cohort, 1, 2)
  # current cohort id = the "right-hand" cohort id if pair, else the single id
  out$CurrID <- ifelse(grepl("/", out$Cohort),
                       as.integer(sub(".*\\/(FT|TR)(\\d+)$", "\\2", out$Cohort)),
                       parse_id(out$Cohort))
  out <- out[order(out$Type, decreasing = TRUE), ]
  out <- do.call(rbind, by(out, out$Type, function(df) df[order(df$CurrID, decreasing = TRUE), ]))
  out$Type <- NULL; out$CurrID <- NULL; rownames(out) <- NULL
  
  # % of total computed after all rows (including singletons)
  total_diff <- sum(out$Diff, na.rm = TRUE)
  if (is.finite(total_diff) && total_diff != 0) {
    out$Pct_of_Total <- round(100 * out$Diff / total_diff, 2)
    remainder <- 100 - sum(out$Pct_of_Total, na.rm = TRUE)
    if (!is.na(remainder) && abs(remainder) >= 0.01) {
      i <- which.max(abs(out$Pct_of_Total))
      out$Pct_of_Total[i] <- out$Pct_of_Total[i] + remainder
    }
  } else {
    out$Pct_of_Total <- NA_real_
  }
  
  # totals row (Semester = NA, A–E blank)
  total_row <- as.data.frame(
    setNames(list("Total", NA_integer_,
                  sum(out[[term1]], na.rm = TRUE),
                  sum(out[[term2]], na.rm = TRUE),
                  sum(out$Diff, na.rm = TRUE),
                  ifelse(is.finite(total_diff) && total_diff != 0, 100, NA),
                  NA, NA, NA, NA, NA),
             c("Cohort","Semester", term1, term2, "Diff","Pct_of_Total","A","B","C","D","E")),
    stringsAsFactors = FALSE
  )
  out <- rbind(out, total_row); rownames(out) <- NULL
  out
}

fcbud_waterfall <- build_fcbud_waterfall(
  fcb_ud_df = fcb_ud_df,
  fcbud_decomp_table = fcbud_decomp_table,
  term1 = "F23", term2 = "F24",
  term_starts = term_starts, terms = terms,
  include_singletons = TRUE, verbose = TRUE
)

save.image(here::here("_outputs/Rdata/data_image1.RData"))

