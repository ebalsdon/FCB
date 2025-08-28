# --- Libraries ---
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(purrr)

# --- Configuration ---
input_folder <- here::here("_inputs/CARB")
file_pattern <- "\\.csv$"
desired_order <- c("FCBud", "FCBpre", "AL", "EDU", "ENG", "HHS", "PSFA", "SCI", "UNDCL", "NE")

base_map <- c(
  "Arts & Letters" = "AL",
  "Education" = "EDU",
  "Engineering" = "ENG",
  "Health & Human Services" = "HHS",
  "Professional Studies & Fine Arts" = "PSFA",
  "Sciences" = "SCI",
  "DFASS" = "UNDCL",
  "No Data" = "NE"
)

# --- Output containers ---
transition_matrices <- list()
transition_probabilities <- list()

# --- File loop ---
files <- list.files(path = input_folder, pattern = file_pattern, full.names = TRUE)

for (file_path in files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))  # e.g., FT21_4_5
  
  # Step 1: Read the modified Tableau export
  df <- read_tsv(file_path, locale = locale(encoding = "UTF-16"), show_col_types = FALSE)
  
  # Step 2: Extract and map category names and flags
  df_clean <- df %>%
    select(from_college = `Migration Row 1 Dimension`,
           from_flag    = `Migration Row 2 Dimension`,
           to_college   = `Migration Row 3 Dimension`,
           to_flag      = `Migration Row 4 Dimension`,
           count        = `Number of Records`) %>%
    mutate(
      count = as.numeric(count),
      
      from = case_when(
        from_college == "Fowler College of Business" & from_flag == "Pre-Major" ~ "FCBpre",
        from_college == "Fowler College of Business" & from_flag == "No Data"   ~ "FCBud",
        TRUE ~ base_map[from_college]
      ),
      to = case_when(
        to_college == "Fowler College of Business" & to_flag == "Pre-Major" ~ "FCBpre",
        to_college == "Fowler College of Business" & to_flag == "No Data"   ~ "FCBud",
        TRUE ~ base_map[to_college]
      )
    ) %>%
    filter(!is.na(from), !is.na(to), !is.na(count))
  
  # Step 3: Complete and pivot
  df_completed <- df_clean %>%
    group_by(from, to) %>%
    summarize(count = sum(count), .groups = "drop") %>%
    complete(from = desired_order, to = desired_order, fill = list(count = 0))
  
  mat <- df_completed %>%
    pivot_wider(names_from = to, values_from = count, values_fill = 0) %>%
    arrange(factor(from, levels = desired_order)) %>%
    select(from, all_of(desired_order)) %>%
    column_to_rownames("from") %>%
    as.matrix()
  
  # Step 4: Normalize
  prob_mat <- sweep(mat, 1, rowSums(mat), FUN = "/")
  prob_mat[is.na(prob_mat)] <- 0
  prob_mat <- round(prob_mat * 100, 2)  # now formatted as percentage
  
  # Step 5: Store
  transition_matrices[[file_name]] <- mat
  transition_probabilities[[file_name]] <- prob_mat
}

# --- Save results ---
#save(transition_matrices, transition_probabilities, 
#     file = "~/_FCB/Summer25/Markov/Rdata/matrices.RData")

# --- Dashboard tables ---
library(glue)

# --- Define cohort and term structure ---
cohorts <- c(paste0("FT", 17:26), paste0("TR", 19:26))
terms <- c("F17", "S18", "F18", "S19", "F19", "S20", "F20", "S21",
           "F21", "S22", "F22", "S23", "F23", "S24", "F24", "S25", 
           "F25", "S26", "F26")

term_starts <- c(
  FT17 = "F17", FT18 = "F18", FT19 = "F19", FT20 = "F20", FT21 = "F21", 
  FT22 = "F22", FT23 = "F23", FT24 = "F24", FT25 = "F25", FT26 = "F26",
  TR19 = "F19", TR20 = "F20", TR21 = "F21", TR22 = "F22", TR23 = "F23",
  TR24 = "F24", TR25 = "F25", TR26 = "F26"
)

# --- Initialize output matrices ---
fcb_pre_matrix <- matrix(NA, nrow = length(cohorts), ncol = length(terms),
                         dimnames = list(cohorts, terms))
fcb_ud_matrix <- fcb_pre_matrix

# --- Main loop: populate from row sums ---
for (matrix_name in names(transition_matrices)) {
  parts <- strsplit(matrix_name, "_")[[1]]
  cohort <- parts[1]
  from_sem <- as.integer(parts[2])
  to_sem   <- as.integer(parts[3])
  
  if (!cohort %in% names(term_starts)) next
  
  term_start_label <- term_starts[[cohort]]
  start_index <- match(term_start_label, terms)
  term_index  <- start_index + to_sem - 2   # <-- your chosen mapping
  
  # ✅ NEW: never write before entry (e.g., S23 for FT23)
  if (term_index < start_index) {
    message(glue::glue("⏭️ Skipping pre-entry write: {matrix_name} -> {terms[term_index]} for {cohort}"))
    next
  }
  
  if (term_index <= ncol(fcb_pre_matrix)) {
    mat <- transition_matrices[[matrix_name]]
    fcb_pre_matrix[cohort, term_index] <- sum(mat["FCBpre", ])
    fcb_ud_matrix [cohort, term_index] <- sum(mat["FCBud", ])
  }
}


# Step 1: Find the last transition for each cohort
library(dplyr)

# Build a lookup table for (cohort, matrix name, to_sem)
transition_info <- tibble(
  matrix_name = names(transition_matrices),
  parts = strsplit(names(transition_matrices), "_")
) %>%
  mutate(
    cohort = sapply(parts, `[`, 1),
    to_sem = as.integer(sapply(parts, `[`, 3))
  ) %>%
  filter(cohort %in% names(term_starts)) %>%
  group_by(cohort) %>%
  slice_max(to_sem, with_ties = FALSE) %>%
  ungroup()

# Step 2: Loop through the last transitions and fill the table using column sums
for (i in seq_len(nrow(transition_info))) {
  cohort <- transition_info$cohort[i]
  matrix_name <- transition_info$matrix_name[i]
  to_sem <- transition_info$to_sem[i]
  
  start_index <- match(term_starts[[cohort]], terms)
  term_index <- start_index + to_sem - 1  # target term reached
  
  if (term_index >= 1 && term_index <= ncol(fcb_pre_matrix)) {
    mat <- transition_matrices[[matrix_name]]
    
    if ("FCBpre" %in% colnames(mat)) {
      fcb_pre_matrix[cohort, term_index] <- sum(mat[, "FCBpre"], na.rm = TRUE)
    }
    if ("FCBud" %in% colnames(mat)) {
      fcb_ud_matrix[cohort, term_index] <- sum(mat[, "FCBud"], na.rm = TRUE)
    }
    
    message(glue::glue("✅ Filled latest term ({terms[term_index]}) for {cohort} using {matrix_name}"))
  }
}

fcb_pre_df <- as.data.frame(fcb_pre_matrix)
fcb_pre_df <- tibble::rownames_to_column(fcb_pre_df, var = "Cohort")

fcb_ud_df <- as.data.frame(fcb_ud_matrix)
fcb_ud_df <- tibble::rownames_to_column(fcb_ud_df, var = "Cohort")

pre_totals <- colSums(fcb_pre_matrix, na.rm = TRUE)
ud_totals <- colSums(fcb_ud_matrix, na.rm = TRUE)

rownames(fcb_ud_df) <- fcb_ud_df$Cohort
fcb_ud_df$Cohort <- NULL  # Optional: remove the redundant column

rownames(fcb_pre_df) <- fcb_pre_df$Cohort
fcb_pre_df$Cohort <- NULL  # Optional: remove the redundant column

# ---------- helper(s) ----------
`%||%` <- function(x, y) if (is.null(x)) y else x

ensure_col <- function(df, col) {
  if (!(col %in% colnames(df))) {
    df[, col] <- NA_real_
    # put it at the end
    df <- df[, unique(c(setdiff(colnames(df), col), col)), drop = FALSE]
  }
  df
}

add_term_totals <- function(df, from_term = "F22") {
  out <- df
  if (!("Total" %in% rownames(out))) {
    out <- rbind(out, Total = rep(NA_real_, ncol(out)))
  }
  j0 <- match(from_term, colnames(out))
  if (!is.na(j0)) {
    rows <- setdiff(rownames(out), "Total")
    cols <- j0:ncol(out)
    out["Total", cols] <- colSums(out[rows, cols, drop = FALSE], na.rm = TRUE)
  }
  out
}

# ---------- Total row for F22 onward (includes forecasts) ----------
fcb_pre_df <- add_term_totals(fcb_pre_df, from_term = "F22")
fcb_ud_df  <- add_term_totals(fcb_ud_df,  from_term = "F22")

# ---- helpers ----
add_term_totals_range <- function(df, from_term = "F22", to_term = "S25") {
  out <- df
  if (!("Total" %in% rownames(out))) {
    out <- rbind(out, Total = rep(NA_real_, ncol(out)))
  }
  j0 <- match(from_term, colnames(out))
  j1 <- match(to_term,   colnames(out))
  if (!is.na(j0) && !is.na(j1) && j0 <= j1) {
    rows <- setdiff(rownames(out), "Total")
    cols <- j0:j1
    out["Total", cols] <- colSums(out[rows, cols, drop = FALSE], na.rm = TRUE)
  }
  out
}

# ---- add totals F22..S25 to clean (no-forecast) tables ----
fcb_pre_df <- add_term_totals_range(fcb_pre_df, from_term = "F22", to_term = "S25")
fcb_ud_df  <- add_term_totals_range(fcb_ud_df,  from_term = "F22", to_term = "S25")

fcb_tot_df <- fcb_pre_df + fcb_ud_df


# ---- save clean snapshot (no forecasts) ----
save(fcb_pre_df, fcb_ud_df,
     file = here::here("_outputs/Rdata", "cohorts_noforecast_dfs.RData"))




