# load packages needed by ALL chunks
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(scales)
library(reactable)
library(stringr)
library(gt)
library(forcats)   # << this is the missing one

# 1) Prepare (order columns, tidy percent)
tbl <- fcbud_waterfall %>%
  mutate(
    Semester = as.integer(Semester),
    Group = ifelse(str_starts(Cohort, "FT"), "FT cohorts", 
                   ifelse(Cohort == "Total", "Totals", "TR cohorts")),
    Pct_of_Total = Pct_of_Total/100  # convert to proportion if it's in %
  ) %>%
  select(Group, Cohort, Semester, F23, F24, Diff, Pct_of_Total, A, B, C, D, E)

# 2) Build a handsome gt table
gt(tbl, groupname_col = "Group") %>%
  tab_header(
    title = md("**FCBud change: F23 → F24**"),
    subtitle = "Waterfall rows with factor decomposition (A–E)"
  ) %>%
  tab_spanner(label = "Factors (Δ decomposition)", columns = c(A, B, C, D, E)) %>%
  cols_label(
    Cohort = "Cohort (pair)",
    Semester = "S",
    F23 = "F23",
    F24 = "F24",
    Diff = "Δ",
    Pct_of_Total = "% of total",
    A = "Admissions", B = "% FCB", C = "Retention",
    D = "Net transfers", E = "Promotion"
  ) %>%
  fmt_number(columns = c(F23, F24, Diff, A:E), decimals = 0, use_seps = TRUE) %>%
  fmt_percent(columns = Pct_of_Total, decimals = 1) %>%
  data_color(
    columns = Diff,
    colors = scales::col_bin(
      palette = c("#d95f02","#f0f0f0","#1b9e77"),
      domain = c(min(tbl$Diff, na.rm=TRUE), max(tbl$Diff, na.rm=TRUE)),
      bins = c(-Inf, 0, Inf)
    )
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_row_groups(groups = "Totals")
  ) %>%
  tab_options(table.font.size = px(14)) %>%
  opt_row_striping()




tbl <- fcbud_waterfall %>%
  mutate(
    Semester = as.integer(Semester),
    Pct_of_Total = Pct_of_Total/100,
    Label = paste0(Cohort, " (S", Semester, ")")
  ) %>%
  select(Label, F23, F24, Diff, Pct_of_Total, A, B, C, D, E)

reactable(
  tbl,
  striped = TRUE,
  highlight = TRUE,
  filterable = TRUE,
  searchable = TRUE,
  defaultPageSize = 20,
  defaultColDef = colDef(align = "right"),
  columns = list(
    Label = colDef(align = "left", sticky = "left", width = 180),
    F23 = colDef(format = colFormat(separators = TRUE, digits = 0)),
    F24 = colDef(format = colFormat(separators = TRUE, digits = 0)),
    Diff = colDef(format = colFormat(separators = TRUE, digits = 0),
                  style = function(value) ifelse(value >= 0, list(color="#1b9e77"), list(color="#d95f02"))),
    Pct_of_Total = colDef(name = "% of total", format = colFormat(percent = TRUE, digits = 1)),
    A = colDef(name = "Admissions", format = colFormat(separators = TRUE, digits = 0)),
    B = colDef(name = "% FCB", format = colFormat(separators = TRUE, digits = 0)),
    C = colDef(name = "Retention", format = colFormat(separators = TRUE, digits = 0)),
    D = colDef(name = "Net transfers", format = colFormat(separators = TRUE, digits = 0)),
    E = colDef(name = "Promotion", format = colFormat(separators = TRUE, digits = 0))
  )
)



library(gt)
library(scales)
library(dplyr)

make_heat_gt <- function(df, title, forecast_cols = c("F25","S26"),
                         palette = c("#f7fbff", "#6baed6", "#08306b")) {

  # Convert to data frame with Cohort column
  mat <- as.data.frame(df, check.names = FALSE)
  mat[] <- lapply(mat, function(x) suppressWarnings(as.numeric(x)))
  tbl <- tibble::tibble(Cohort = rownames(mat), !!!mat)
  
  # Identify numeric columns for formatting/coloring
  num_cols <- names(tbl)[sapply(tbl, is.numeric)]
  
  # Color domain from full table
  rng <- range(as.matrix(mat), na.rm = TRUE)
  
  gt(tbl) %>%
    tab_header(title = md(paste0("**", title, "**"))) %>%
    fmt_number(columns = all_of(num_cols), decimals = 0, sep_mark = ",") %>%
    fmt_missing(columns = all_of(num_cols), missing_text = "") %>%
    data_color(
      columns = all_of(num_cols),
      colors = col_numeric(palette = palette, domain = rng, na.color = "white")
    ) %>%
    tab_style(
      style = list(cell_text(style = "italic", weight = "bold")),
      locations = cells_body(columns = tidyselect::any_of(forecast_cols))
    ) %>%
    opt_row_striping() %>%
    cols_label(.list = setNames(names(tbl), names(tbl))) %>%
    tab_options(table.font.size = px(14), data_row.padding = px(4))
}



make_heat_gt(fcb_ud_df,  "FCBud by cohort & term (heatmap)", forecast_cols = c("F25","S26"))
make_heat_gt(fcb_pre_df, "FCBpre by cohort & term (heatmap)", forecast_cols = c("F25","S26"))


