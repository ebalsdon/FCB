# Diagnostics: (1) F25->F26 cohort flows, (2) S27 forecast under n_reference = 1,2,3,
# (3) flow decomposition of the S27 forecast, (4) sensitivity to reference-cohort assumptions.
# Read-only with respect to the pipeline: loads data_image1.RData, writes CSVs to _outputs/analysis.

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(tibble) })
options(width = 200, dplyr.width = Inf, pillar.width = Inf)

e <- new.env()
load(here::here("_outputs/Rdata/data_image1.RData"), envir = e)
TM <- e$transition_matrices; TP <- e$transition_probabilities
terms <- e$terms; term_starts <- e$term_starts
pre <- e$fcb_pre_df; ud <- e$fcb_ud_df
states <- colnames(TP[[1]])
colleges <- setdiff(states, c("FCBud", "FCBpre", "NE"))
out_dir <- here::here("_outputs/analysis"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cohorts <- setdiff(rownames(ud), "Total")

sem_at <- function(cohort, term) match(term, terms) - match(term_starts[[cohort]], terms) + 1L
max_sem <- function(cohort) if (startsWith(cohort, "FT")) 16L else 12L
val <- function(df, cohort, term) { v <- df[cohort, term]; if (is.na(v)) 0 else as.numeric(v) }
write_out <- function(df, name) write.csv(df, file.path(out_dir, name), row.names = FALSE)

# ---------------------------------------------------------------- reference matrices
ref_P <- function(cohort, step_from, n_ref) {
  grp <- substr(cohort, 1, 2); yr <- as.integer(substr(cohort, 3, 4)); out <- list()
  for (k in 1:40) {
    cand <- paste0(grp, sprintf("%02d", yr - k))
    key <- paste0(cand, "_", step_from, "_", step_from + 1L)
    if (key %in% names(TP)) {
      P <- TP[[key]]; if (max(P, na.rm = TRUE) > 1.0001) P <- P / 100
      P <- P[states, states, drop = FALSE]; P[is.na(P)] <- 0
      out[[cand]] <- P
      if (length(out) == n_ref) break
    }
  }
  out
}
x_F26 <- function(cohort) {
  s <- sem_at(cohort, "F26"); key <- if (s == 1L) paste0(cohort, "_1_1") else paste0(cohort, "_", s - 1L, "_", s)
  if (!(key %in% names(TM))) return(NULL)
  x <- colSums(TM[[key]], na.rm = TRUE)[states]; x[is.na(x)] <- 0; x
}
eligible <- Filter(function(c) {
  s <- sem_at(c, "F26"); !is.na(s) && s >= 1 && (s + 1L) <= max_sem(c) && !is.null(x_F26(c))
}, cohorts)

# ================================================================ 1. F25 -> F26 actuals
cat("\n==== 1. ACTUAL HEADCOUNT CHANGE F25 -> F26 ====\n")
act <- bind_rows(lapply(cohorts, function(c) {
  tibble(Cohort = c,
         pre_F25 = val(pre, c, "F25"), pre_F26 = val(pre, c, "F26"),
         ud_F25 = val(ud, c, "F25"), ud_F26 = val(ud, c, "F26"))
})) %>% mutate(
  role = ifelse(Cohort %in% c("FT26", "TR26"), "Incoming", ifelse(startsWith(Cohort, "FT"), "Continuing FT", "Continuing TR")),
  d_pre = pre_F26 - pre_F25, d_ud = ud_F26 - ud_F25, tot_F25 = pre_F25 + ud_F25, tot_F26 = pre_F26 + ud_F26, d_tot = tot_F26 - tot_F25)
print(as.data.frame(act), row.names = FALSE)
cat("\nBy role:\n")
by_role <- act %>% group_by(role) %>% summarise(across(c(pre_F25, pre_F26, d_pre, ud_F25, ud_F26, d_ud, tot_F25, tot_F26, d_tot), sum), .groups = "drop")
print(as.data.frame(by_role), row.names = FALSE)
cat("\nGrand total:\n"); print(colSums(act[, c("pre_F25", "pre_F26", "d_pre", "ud_F25", "ud_F26", "d_ud", "tot_F25", "tot_F26", "d_tot")]))
write_out(act, "F25_F26_cohort_change.csv")

# ---------------------------------------------------------------- flow decomposition over a period
period_flows <- function(c, t_from, t_to) {
  s <- sem_at(c, t_from); n_steps <- match(t_to, terms) - match(t_from, terms)
  if (s < 1) return(NULL)
  rows <- list()
  for (j in 0:(n_steps - 1)) {
    key <- paste0(c, "_", s + j, "_", s + j + 1L)
    if (!(key %in% names(TM))) next
    m <- TM[[key]][states, states, drop = FALSE]
    for (tg in c("FCBpre", "FCBud")) {
      oth <- setdiff(c("FCBpre", "FCBud"), tg)
      rows[[length(rows) + 1]] <- tibble(Cohort = c, step = j + 1L, target = tg,
        in_fcb = m[oth, tg], in_college = sum(m[colleges, tg]), in_NE = m["NE", tg],
        out_fcb = m[tg, oth], out_college = sum(m[tg, colleges]), out_NE = m[tg, "NE"])
    }
  }
  if (length(rows)) bind_rows(rows) else NULL
}
flow_table <- function(t_from, t_to) {
  fl <- bind_rows(lapply(cohorts, function(c) period_flows(c, t_from, t_to)))
  fl %>% group_by(Cohort, target) %>% summarise(across(c(in_fcb, in_college, in_NE, out_fcb, out_college, out_NE), sum), .groups = "drop")
}
cat("\n==== 1b. FLOW RECONCILIATION F25 -> F26 (continuing cohorts, two steps: F25->S26->F26) ====\n")
fl <- flow_table("F25", "F26")
recon <- fl %>% left_join(
  bind_rows(act %>% transmute(Cohort, target = "FCBpre", delta = d_pre), act %>% transmute(Cohort, target = "FCBud", delta = d_ud)),
  by = c("Cohort", "target")) %>%
  mutate(net_flows = in_fcb + in_college + in_NE - out_fcb - out_college - out_NE, residual = delta - net_flows)
tot_recon <- recon %>% group_by(target) %>% summarise(across(c(in_fcb, in_college, in_NE, out_fcb, out_college, out_NE, net_flows, delta, residual), sum), .groups = "drop")
print(as.data.frame(tot_recon), row.names = FALSE)
cat("\nCohorts with non-zero residual (flows do not fully reconcile, e.g. beyond tracked window):\n")
print(as.data.frame(recon %>% filter(abs(residual) > 0.5) %>% select(Cohort, target, delta, net_flows, residual)), row.names = FALSE)
write_out(recon, "F25_F26_flow_reconciliation.csv")

cat("\n==== 1c. FCB net flows by cohort, FCBpre+FCBud combined (F25 -> F26) ====\n")
by_c <- recon %>% group_by(Cohort) %>% summarise(
  promo_pre_to_ud_net = 0, in_college = sum(in_college), in_NE = sum(in_NE), out_college = sum(out_college), out_NE = sum(out_NE),
  net = sum(net_flows), delta = sum(delta), residual = sum(residual), .groups = "drop") %>% select(-promo_pre_to_ud_net)
print(as.data.frame(by_c), row.names = FALSE)
write_out(by_c, "F25_F26_flows_by_cohort.csv")

cat("\n==== 1d. CALENDAR-STEP RATES (all cohorts with matrices; mix of stages differs by step) ====\n")
step_rates <- function(t_from, t_to) {
  rows <- list()
  for (c in cohorts) {
    s <- sem_at(c, t_from); key <- paste0(c, "_", s, "_", s + 1L)
    if (s >= 1 && key %in% names(TM)) rows[[length(rows) + 1]] <- tibble(Cohort = c, mat = list(TM[[key]][states, states]))
  }
  tot <- Reduce(`+`, lapply(rows, function(r) r$mat[[1]]))
  rs <- function(from) { r <- tot[from, ]; n <- sum(r)
    list(n = n, stay = unname(r[from] / n),
         promote = if (from == "FCBpre") unname(r["FCBud"] / n) else NA_real_,
         to_college = sum(r[colleges]) / n, to_NE = unname(r["NE"] / n)) }
  p <- rs("FCBpre"); u <- rs("FCBud")
  inflow <- sum(tot[colleges, c("FCBpre", "FCBud")])
  tibble(step = paste0(t_from, "->", t_to),
         pre_n = p$n, pre_stay = p$stay, pre_promote = p$promote, pre_to_college = p$to_college, pre_to_NE = p$to_NE,
         ud_n = u$n, ud_stay = u$stay, ud_to_college = u$to_college, ud_to_NE = u$to_NE,
         in_from_colleges = inflow, in_from_NE = sum(tot["NE", c("FCBpre", "FCBud")]), n_cohorts = length(rows))
}
sr <- bind_rows(step_rates("F24", "S25"), step_rates("S25", "F25"), step_rates("F25", "S26"), step_rates("S26", "F26"))
print(as.data.frame(sr %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
write_out(sr, "calendar_step_rates.csv")

# ================================================================ 2. S27 under n_reference
cat("\n==== 2. S27 FORECAST BY n_reference ====\n")
fc_cohort <- function(c, n_ref) {
  s <- sem_at(c, "F26"); x <- x_F26(c); refs <- ref_P(c, s, n_ref)
  if (!length(refs)) return(NULL)
  P <- Reduce(`+`, refs) / length(refs); x1 <- as.numeric(x %*% P); names(x1) <- states
  tibble(Cohort = c, n_ref_requested = n_ref, n_ref_used = length(refs), refs = paste(names(refs), collapse = ","),
         F26_pre = x["FCBpre"], F26_ud = x["FCBud"], S27_pre = x1["FCBpre"], S27_ud = x1["FCBud"])
}
runs <- bind_rows(lapply(1:3, function(n) bind_rows(lapply(eligible, fc_cohort, n_ref = n))))
runs <- runs %>% mutate(S27_tot = S27_pre + S27_ud, F26_tot = F26_pre + F26_ud)
beyond <- setdiff(cohorts, eligible)
cat("Cohorts forced to 0 / not forecast (beyond program length or no F26 matrix):", paste(beyond, collapse = ", "), "\n")
chk <- runs %>% filter(n_ref_requested == 3) %>% transmute(Cohort, saved_pre = pre[Cohort, "S27"], saved_ud = ud[Cohort, "S27"],
        d_pre = round(S27_pre) - saved_pre, d_ud = round(S27_ud) - saved_ud)
cat("Reproduction check vs saved n_reference=3 forecast (max abs diff pre/ud):", max(abs(chk$d_pre)), max(abs(chk$d_ud)), "\n")
agg <- runs %>% group_by(n_ref_requested) %>% summarise(F26_pre = sum(F26_pre), F26_ud = sum(F26_ud), S27_pre = sum(S27_pre), S27_ud = sum(S27_ud),
        S27_tot = sum(S27_tot), F26_tot = sum(F26_tot), .groups = "drop") %>% mutate(chg_tot = S27_tot - F26_tot)
print(as.data.frame(agg %>% mutate(across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)
cat("\nPer-cohort S27 (pre / ud), by n_reference:\n")
wide <- runs %>% select(Cohort, n_ref_requested, n_ref_used, S27_pre, S27_ud, F26_pre, F26_ud) %>%
  pivot_wider(names_from = n_ref_requested, values_from = c(n_ref_used, S27_pre, S27_ud)) %>%
  mutate(range_ud = pmax(S27_ud_1, S27_ud_2, S27_ud_3) - pmin(S27_ud_1, S27_ud_2, S27_ud_3),
         range_pre = pmax(S27_pre_1, S27_pre_2, S27_pre_3) - pmin(S27_pre_1, S27_pre_2, S27_pre_3))
print(as.data.frame(wide %>% mutate(across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)
write_out(runs, "S27_by_n_reference_long.csv"); write_out(wide, "S27_by_n_reference_wide.csv")

# ---------------------------------------------------------------- reference calendar terms
cat("\n==== 2b. WHICH CALENDAR TRANSITIONS THE REFERENCES REPRESENT (n=3) ====\n")
ref_cal <- bind_rows(lapply(eligible, function(c) {
  s <- sem_at(c, "F26"); refs <- names(ref_P(c, s, 3))
  tibble(Cohort = c, step = paste0(s, "->", s + 1L), ref = refs,
         ref_from_term = vapply(refs, function(r) terms[match(term_starts[[r]], terms) + s - 1L], ""),
         ref_to_term = vapply(refs, function(r) terms[match(term_starts[[r]], terms) + s], ""))
}))
print(as.data.frame(ref_cal %>% group_by(Cohort, step) %>% summarise(ref_transitions = paste0(ref, "(", ref_from_term, ">", ref_to_term, ")", collapse = "  "), .groups = "drop")), row.names = FALSE)
write_out(ref_cal, "S27_reference_calendar_terms.csv")

# ================================================================ 3. flow decomposition of S27 (n=3)
cat("\n==== 3. FLOW DECOMPOSITION OF S27 FORECAST (n_reference = 3), headcounts ====\n")
cat_of <- function(i, j) {
  if (i == "FCBpre" && j == "FCBpre") return("pre stays pre")
  if (i == "FCBpre" && j == "FCBud") return("pre promotes to ud")
  if (i == "FCBud" && j == "FCBud") return("ud stays ud")
  if (i == "FCBud" && j == "FCBpre") return("ud -> pre")
  if (i %in% colleges && j %in% c("FCBpre", "FCBud")) return(paste0("in-transfer from college -> ", j))
  if (i == "NE" && j %in% c("FCBpre", "FCBud")) return(paste0("re-entry NE -> ", j))
  if (i %in% c("FCBpre", "FCBud") && j %in% colleges) return(paste0("out-transfer ", i, " -> college"))
  if (i %in% c("FCBpre", "FCBud") && j == "NE") return(paste0("attrition/grad ", i, " -> NE"))
  NA_character_
}
flow_s27 <- function(n_ref) bind_rows(lapply(eligible, function(c) {
  s <- sem_at(c, "F26"); x <- x_F26(c); refs <- ref_P(c, s, n_ref); P <- Reduce(`+`, refs) / length(refs)
  g <- expand.grid(i = states, j = states, stringsAsFactors = FALSE)
  g$flow <- mapply(function(i, j) x[i] * P[i, j], g$i, g$j)
  g$cohort <- c; g$category <- mapply(cat_of, g$i, g$j); g$x_i <- x[g$i]
  g
}))
fs <- flow_s27(3)
fcb_flows <- fs %>% filter(!is.na(category)) %>% group_by(category) %>% summarise(headcount = sum(flow), .groups = "drop") %>% arrange(desc(headcount))
print(as.data.frame(fcb_flows %>% mutate(headcount = round(headcount, 1))), row.names = FALSE)
cat("\nSum of flows into FCBpre / FCBud (=S27):",
    round(sum(fs$flow[fs$j == "FCBpre"]), 1), round(sum(fs$flow[fs$j == "FCBud"]), 1), "\n")
mass <- bind_rows(lapply(eligible, function(c) {
  s <- sem_at(c, "F26"); x <- x_F26(c); P <- Reduce(`+`, ref_P(c, s, 3)) / length(ref_P(c, s, 3))
  tibble(Cohort = c, state = states, x = x, row_sum = rowSums(P)[states], lost = x * (1 - rowSums(P)[states]))
}))
cat("Headcount whose reference row has no observations (mass silently dropped); FCB states / other states:\n")
print(c(FCB_states = sum(mass$lost[mass$state %in% c("FCBpre", "FCBud")]), other_states = sum(mass$lost[!mass$state %in% c("FCBpre", "FCBud")])))
print(as.data.frame(mass %>% filter(abs(lost) > 0.5) %>% mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
write_out(fcb_flows, "S27_flow_categories.csv")

cat("\nPer-cohort net change F26 -> S27 by component (n=3):\n")
pc <- fs %>% filter(!is.na(category)) %>% mutate(grp = case_when(
    grepl("attrition", category) ~ "to_NE", grepl("out-transfer", category) ~ "out_college",
    grepl("in-transfer", category) ~ "in_college", grepl("re-entry", category) ~ "in_NE", TRUE ~ "within_FCB")) %>%
  group_by(cohort, grp) %>% summarise(h = sum(flow), .groups = "drop")
pc_w <- pc %>% pivot_wider(names_from = grp, values_from = h, values_fill = 0)
x0 <- runs %>% filter(n_ref_requested == 3) %>% select(cohort = Cohort, F26_tot, S27_tot)
pc_w <- left_join(pc_w, x0, by = "cohort") %>% mutate(net = S27_tot - F26_tot)
print(as.data.frame(pc_w %>% mutate(across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)
write_out(pc_w, "S27_net_components_by_cohort.csv")

# ================================================================ 4. sensitivity to reference assumptions
cat("\n==== 4a. RECENCY: aggregate S27 if every cohort used ONLY its k-th most recent reference ====\n")
lag <- bind_rows(lapply(eligible, function(c) {
  s <- sem_at(c, "F26"); x <- x_F26(c); refs <- ref_P(c, s, 5)
  bind_rows(lapply(seq_along(refs), function(k) { x1 <- x %*% refs[[k]]
    tibble(Cohort = c, k = k, ref = names(refs)[k], S27_pre = x1[, "FCBpre"], S27_ud = x1[, "FCBud"]) }))
}))
lag_agg <- lag %>% group_by(k) %>% summarise(cohorts_with_k = n(), S27_pre = sum(S27_pre), S27_ud = sum(S27_ud), .groups = "drop") %>% mutate(S27_tot = S27_pre + S27_ud)
base_sub <- lapply(1:5, function(k) { cs <- lag$Cohort[lag$k == k]; runs %>% filter(n_ref_requested == 1, Cohort %in% cs) %>% summarise(n1_pre = sum(S27_pre), n1_ud = sum(S27_ud)) })
lag_agg <- bind_cols(lag_agg, bind_rows(base_sub)) %>% mutate(ud_minus_n1 = S27_ud - n1_ud, pre_minus_n1 = S27_pre - n1_pre)
print(as.data.frame(lag_agg %>% mutate(across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)
write_out(lag, "S27_by_single_reference_lag.csv")

cat("\n==== 4b. WHICH TRANSITION PROBABILITIES MATTER: flow, headcount at risk (x_i), and cross-reference spread ====\n")
cat("spread_headcount = sum over cohorts of x_i * SD of P[i,t] across up to 5 reference cohorts (headcount-weighted variability);\n")
cat("per_1pp = headcount change in S27 if P[i,t] were 1 percentage point higher in every cohort.\n")
sens <- bind_rows(lapply(eligible, function(c) {
  s <- sem_at(c, "F26"); x <- x_F26(c); refs5 <- ref_P(c, s, 5); refs3 <- refs5[seq_len(min(3, length(refs5)))]
  P3 <- Reduce(`+`, refs3) / length(refs3)
  arr <- simplify2array(refs5)
  sdm <- if (length(refs5) > 1) apply(arr, c(1, 2), sd) else matrix(0, length(states), length(states), dimnames = list(states, states))
  g <- expand.grid(i = states, t = c("FCBpre", "FCBud"), stringsAsFactors = FALSE)
  g$Cohort <- c; g$x_i <- x[g$i]; g$P3 <- mapply(function(i, t) P3[i, t], g$i, g$t)
  g$sd_ref <- mapply(function(i, t) sdm[i, t], g$i, g$t); g$n_refs <- length(refs5)
  g$flow <- g$x_i * g$P3; g$spread <- g$x_i * g$sd_ref; g$per_1pp <- g$x_i * 0.01
  g
}))
sens_cell <- sens %>% group_by(target = t, from = i) %>% summarise(
  flow = sum(flow), spread_linear = sum(spread), spread_rss = sqrt(sum(spread^2)), per_1pp = sum(per_1pp), .groups = "drop")
for (tg in c("FCBud", "FCBpre")) {
  cat("\nTarget", tg, "- top (from-state -> target) cells by headcount-weighted spread:\n")
  print(as.data.frame(sens_cell %>% filter(target == tg) %>% arrange(desc(spread_linear)) %>% head(8) %>%
                        mutate(across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)
}
cat("\nTop (cohort, from-state, target) cells by spread:\n")
print(as.data.frame(sens %>% arrange(desc(spread)) %>% head(15) %>% select(Cohort, from = i, target = t, x_i, P3, sd_ref, n_refs, flow, spread) %>%
                      mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
write_out(sens, "S27_sensitivity_cells.csv"); write_out(sens_cell, "S27_sensitivity_by_cell.csv")

cat("\nSpread by cohort (linear sum across cells, FCBpre and FCBud targets):\n")
sc <- sens %>% group_by(Cohort, t) %>% summarise(spread = sum(spread), flow = sum(flow), .groups = "drop") %>%
  pivot_wider(names_from = t, values_from = c(spread, flow)) %>% arrange(desc(spread_FCBud))
print(as.data.frame(sc %>% mutate(across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)

cat("\nDone. CSVs in", out_dir, "\n")

