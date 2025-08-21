### Native plotly version of waterfall chart
library(dplyr)
library(plotly)
library(scales)

build_plotly_waterfall_topk <- function(fcbud_waterfall,
                                        term1 = "F23",
                                        term2 = "F24",
                                        top_k = 10,
                                        order_by = c("abs_diff","pct","diff"),
                                        y_floor = NULL,     # e.g., 2500 (NULL = baseline from 0)
                                        verbose = TRUE) {
  order_by <- match.arg(order_by)
  
  # Basic checks
  stopifnot(all(c("Cohort","Semester","Diff", term1, term2) %in% names(fcbud_waterfall)))
  
  # Separate out cohort rows (exclude Total)
  cohorts_only <- fcbud_waterfall %>% filter(Cohort != "Total")
  
  # Totals (prefer Total row if present)
  if ("Total" %in% fcbud_waterfall$Cohort) {
    total_t1 <- fcbud_waterfall[[term1]][fcbud_waterfall$Cohort == "Total"][1]
    total_t2 <- fcbud_waterfall[[term2]][fcbud_waterfall$Cohort == "Total"][1]
  } else {
    total_t1 <- sum(cohorts_only[[term1]], na.rm = TRUE)
    total_t2 <- sum(cohorts_only[[term2]], na.rm = TRUE)
  }
  total_delta <- total_t2 - total_t1
  
  # Order key (+ % fallback)
  core <- cohorts_only %>%
    mutate(
      pct_used  = ifelse(is.na(Pct_of_Total), 100 * Diff / total_delta, Pct_of_Total),
      order_key = dplyr::case_when(
        order_by == "abs_diff" ~ abs(Diff),
        order_by == "pct"      ~ abs(pct_used),
        TRUE                   ~ Diff
      )
    ) %>%
    arrange(desc(order_key))
  
  # Top-K + Other
  top_k <- min(top_k, nrow(core))
  top   <- core %>% slice_head(n = top_k)
  
  other_sum <- total_delta - sum(top$Diff, na.rm = TRUE)
  
  make_label <- function(cohort, sem) {
    ifelse(cohort == "Other", "Other", paste0(cohort, " (S", as.integer(sem), ")"))
  }
  
  # Build contributors table then enforce order: |Δ| desc, with "Other" last
  steps <- bind_rows(
    transmute(
      top,
      Cohort, Semester,
      label = make_label(Cohort, Semester),
      value = Diff,
      hover = paste0(
        Cohort, " (S", as.integer(Semester), ")",
        "<br>", term1, ": ", comma(get(term1)),
        " → ", term2, ": ", comma(get(term2)),
        "<br>Δ = ", comma(Diff)
      )
    ),
    tibble(
      Cohort = "Other",
      Semester = NA_integer_,
      label = "Other",
      value = other_sum,
      hover = paste0("Other contributors<br>Δ = ", comma(other_sum))
    )
  ) %>%
    mutate(is_other = Cohort == "Other") %>%
    arrange(desc(abs(value)), is_other) %>%   # desired order
    select(-is_other)
  
  if (verbose) {
    message(sprintf("Top-K: %d of %d. Other = %+d. Sum(steps) = %+d. Full Δ = %+d.",
                    top_k, nrow(core), round(other_sum),
                    round(sum(steps$value, na.rm = TRUE)), round(total_delta)))
  }
  
  # Plotly inputs
  measures <- c("total", rep("relative", nrow(steps)), "total")
  x_vals   <- c(term1, steps$label, term2)
  y_vals   <- c(total_t1, steps$value, total_t2)
  # harden y to numeric; replace NA deltas with 0 to avoid dropped bars
  y_vals <- suppressWarnings(as.numeric(y_vals))
  y_vals[is.na(y_vals)] <- 0
  
  textvals <- c(
    paste0(term1, ": ", comma(round(total_t1))),
    steps$hover,
    paste0(term2, ": ", comma(round(total_t2)))
  )
  
  # Axis range so totals are visible
  y_min <- if (is.null(y_floor)) 0 else y_floor
  y_max <- max(total_t1, total_t2) * 1.05
  
  plot_ly(
    type    = "waterfall",
    x       = x_vals,
    y       = y_vals,
    measure = measures,
    text    = textvals,
    hoverinfo = "text",
    connector = list(line = list(color = "rgba(130,130,130,0.6)", width = 1))
  ) %>%
    style(
      increasing = list(marker = list(color = "#1b9e77", opacity = 0.95)),
      decreasing = list(marker = list(color = "#d95f02", opacity = 0.95)),
      totals     = list(marker = list(color = "#6b6b6b", opacity = 0.95))
    ) %>%
    layout(
      title = paste0("FCBud change — Top-K + Other: ", term1, " → ", term2),
      xaxis = list(
        title = "",
        type = "category",
        categoryorder = "array",
        categoryarray = x_vals   # << force F23 → steps → F24 order
      ),
      yaxis = list(title = "Headcount", range = c(y_min, y_max))
    )
}

wf_plotly <- build_plotly_waterfall_topk(
  fcbud_waterfall,
  term1 = "F23", term2 = "F24",
  top_k = 5,                # try 5 / 8 / 12
  order_by = "abs_diff",     # or "pct" / "diff"
  y_floor = 2500             # NULL for baseline at 0
)
wf_plotly

# Quick debug to see exactly what Plotly receives
print(data.frame(x_vals = c("F23", head(wf_plotly$x$data[[1]]$x, -1), "F24"),
                 measure = c("total", head(wf_plotly$x$data[[1]]$measure, -1), "total"),
                 y_vals = c(wf_plotly$x$data[[1]]$y)))
