# Regime-aware S27 forecast. Policy: from F25>S26 on, switching from other colleges into FCBpre is gated
# (30 units + higher GPA). Old regime = transitions ending at or before F25 (S25>F25 is the last old-regime step).
# Hybrid P for the S27 step: college->FCBpre column from the most recent post-policy same-stage reference
# (the F25>S26 transition); non-switchers remain in their college (row sums preserved); all other cells use the
# n_reference average. Optional pre->UD uplift for stage 3>4 pre students (smaller switcher share => stronger pool).
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(tibble) })
options(width = 220, dplyr.width = Inf, pillar.width = Inf)
e <- new.env(); load(here::here("_outputs/Rdata/data_image1.RData"), envir = e)
TM <- e$transition_matrices; TP <- e$transition_probabilities; terms <- e$terms; term_starts <- e$term_starts
pre <- e$fcb_pre_df; ud <- e$fcb_ud_df
states <- colnames(TP[[1]]); colleges <- setdiff(states, c("FCBud", "FCBpre", "NE"))
out_dir <- here::here("_outputs/analysis"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cohorts <- setdiff(rownames(ud), "Total")
sem_at <- function(c, t) match(t, terms) - match(term_starts[[c]], terms) + 1L
max_sem <- function(c) if (startsWith(c, "FT")) 16L else 12L
from_term <- function(c, s) terms[match(term_starts[[c]], terms) + s - 1L]

getP <- function(key) { P <- TP[[key]]; if (max(P, na.rm = TRUE) > 1.0001) P <- P / 100; P <- P[states, states, drop = FALSE]; P[is.na(P)] <- 0; P }
ref_P <- function(c, s, n_ref) {
  grp <- substr(c, 1, 2); yr <- as.integer(substr(c, 3, 4)); out <- list()
  for (k in 1:40) { cand <- paste0(grp, sprintf("%02d", yr - k)); key <- paste0(cand, "_", s, "_", s + 1L)
    if (key %in% names(TP)) { out[[cand]] <- getP(key); if (length(out) == n_ref) break } }
  out
}
x_F26 <- function(c) { s <- sem_at(c, "F26"); key <- if (s == 1L) paste0(c, "_1_1") else paste0(c, "_", s - 1L, "_", s)
  if (!(key %in% names(TM))) return(NULL); x <- colSums(TM[[key]], na.rm = TRUE)[states]; x[is.na(x)] <- 0; x }
eligible <- Filter(function(c) { s <- sem_at(c, "F26"); !is.na(s) && s >= 1 && s + 1L <= max_sem(c) && !is.null(x_F26(c)) }, cohorts)
POST_FROM <- "F25"   # first term of the first post-policy transition (F25>S26)

# ---------------------------------------------------------------- uplift estimate: promotion vs switcher share at stage 3>4
cat("\n==== U. Estimating the pre->UD uplift at stage 3>4 (spring) from switcher share of the pre pool ====\n")
u_rows <- lapply(setdiff(cohorts, c("FT17")), function(c) {
  if (!startsWith(c, "FT") && !startsWith(c, "TR")) return(NULL)
  k34 <- paste0(c, "_3_4"); k12 <- paste0(c, "_1_2"); k23 <- paste0(c, "_2_3")
  if (!all(c(k34, k12, k23) %in% names(TM))) return(NULL)
  sw <- sum(TM[[k12]][colleges, "FCBpre"]) + sum(TM[[k23]][colleges, "FCBpre"])
  npre <- sum(TM[[k34]]["FCBpre", ]); pr <- TM[[k34]]["FCBpre", "FCBud"] / npre
  tibble(Cohort = c, grp = substr(c, 1, 2), from34 = from_term(c, 3), npre = npre, switchers = sw, share = sw / npre, promote = pr)
})
U <- bind_rows(u_rows) %>% filter(grp == "FT", match(from34, terms) < match(POST_FROM, terms))   # pre-policy FT cohorts only
print(as.data.frame(U %>% mutate(across(c(share, promote), ~ round(.x, 3)))), row.names = FALSE)
fit <- lm(promote ~ share, data = U, weights = npre); sm <- summary(fit)
cat("\nWLS fit: promote = a + b*share\n"); print(round(coef(sm), 4))
FT25_share <- {
  c <- "FT25"; sw <- sum(TM[["FT25_1_2"]][colleges, "FCBpre"]) + sum(TM[["FT25_2_3"]][colleges, "FCBpre"])
  sw / sum(TM[["FT25_3_4"]][ "FCBpre", ], na.rm = TRUE)
}
FT25_pool <- sum(TM[["FT25_2_3"]][, "FCBpre"]); cat("FT25 pre at F26 (to-pre column of FT25_2_3):", FT25_pool, "\n")
FT25_sw <- sum(TM[["FT25_1_2"]][colleges, "FCBpre"]) + sum(TM[["FT25_2_3"]][colleges, "FCBpre"])
FT25_share_now <- FT25_sw / FT25_pool
share_ref <- weighted.mean(U$share, U$npre)
cat(sprintf("FT25 observed switcher share: %.3f (switchers=%d, pool=%d); pre-policy pooled share: %.3f\n", FT25_share_now, FT25_sw, FT25_pool, share_ref))
b <- coef(fit)["share"]; b_se <- sm$coefficients["share", 2]
UPLIFT_POINT <- as.numeric(-b * (share_ref - FT25_share_now)); UPLIFT_SE <- as.numeric(b_se * (share_ref - FT25_share_now))
cat(sprintf("Implied uplift in P(pre->UD, 3>4): %.4f  (SE %.4f; regression on %d cohorts -> treat as a scenario knob)\n", UPLIFT_POINT, UPLIFT_SE, nrow(U)))
uplift_levels <- c(none = 0, point = max(UPLIFT_POINT, 0), high = max(UPLIFT_POINT, 0) + UPLIFT_SE)

# ---------------------------------------------------------------- hybrid forecast
make_P <- function(c, s, n_ref, regime_aware = TRUE, uplift = 0) {
  refs <- ref_P(c, s, n_ref); if (!length(refs)) return(NULL)
  P <- Reduce(`+`, refs) / length(refs); info <- list(post_ref = NA_character_)
  if (regime_aware) {
    post1 <- ref_P(c, s, 1); r1 <- names(post1)[1]
    is_post <- match(from_term(r1, s), terms) >= match(POST_FROM, terms)
    if (is_post) {
      Pp <- post1[[1]]
      for (i in colleges) { d <- P[i, "FCBpre"] - Pp[i, "FCBpre"]; P[i, "FCBpre"] <- Pp[i, "FCBpre"]; P[i, i] <- P[i, i] + d }
      info$post_ref <- r1
    }
  }
  if (uplift != 0 && s == 3L) { P["FCBpre", "FCBud"] <- P["FCBpre", "FCBud"] + uplift; P["FCBpre", "FCBpre"] <- P["FCBpre", "FCBpre"] - uplift }
  attr(P, "info") <- info; P
}
run <- function(n_ref, regime_aware, uplift = 0) bind_rows(lapply(eligible, function(c) {
  s <- sem_at(c, "F26"); x <- x_F26(c); P <- make_P(c, s, n_ref, regime_aware, uplift); if (is.null(P)) return(NULL)
  x1 <- as.numeric(x %*% P); names(x1) <- states
  tibble(Cohort = c, S27_pre = x1["FCBpre"], S27_ud = x1["FCBud"], F26_pre = x["FCBpre"], F26_ud = x["FCBud"],
         in_col_pre = sum(x[colleges] * P[colleges, "FCBpre"]), post_ref = attr(P, "info")$post_ref)
}))

scen <- bind_rows(lapply(1:3, function(n) bind_rows(
  run(n, FALSE, 0) %>% mutate(scenario = "old (as published)"),
  run(n, TRUE, 0) %>% mutate(scenario = "regime-aware"),
  run(n, TRUE, uplift_levels["point"]) %>% mutate(scenario = "regime-aware + uplift (point)"),
  run(n, TRUE, uplift_levels["high"]) %>% mutate(scenario = "regime-aware + uplift (+1SE)")) %>% mutate(n_ref = n)))
if (any(is.na(scen$post_ref[grepl("regime", scen$scenario)]))) cat("NOTE: some cohorts had no post-policy reference\n")
agg <- scen %>% group_by(scenario, n_ref) %>% summarise(S27_pre = sum(S27_pre), S27_ud = sum(S27_ud), in_col_pre = sum(in_col_pre), .groups = "drop") %>%
  mutate(S27_tot = S27_pre + S27_ud)
cat("\n==== R1. S27 totals by scenario and n_reference ====\n")
print(as.data.frame(agg %>% arrange(n_ref, scenario) %>% mutate(across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)

cat("\n==== R2. Spread across n_reference (max-min) by scenario: how much does the n choice still matter? ====\n")
print(as.data.frame(agg %>% group_by(scenario) %>% summarise(pre_range = max(S27_pre) - min(S27_pre), ud_range = max(S27_ud) - min(S27_ud),
  tot_range = max(S27_tot) - min(S27_tot), .groups = "drop") %>% mutate(across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)

cat("\n==== R3. Per cohort, n=3: old vs regime-aware (pre), and the college->pre inflow ====\n")
w <- scen %>% filter(n_ref == 3, scenario %in% c("old (as published)", "regime-aware", "regime-aware + uplift (point)")) %>%
  select(Cohort, scenario, S27_pre, S27_ud, in_col_pre) %>%
  pivot_wider(names_from = scenario, values_from = c(S27_pre, S27_ud, in_col_pre))
print(as.data.frame(w %>% mutate(across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)

cat("\n==== R4. Regime-aware FT26 pre-major S27 across n (the former main swing) ====\n")
print(as.data.frame(scen %>% filter(Cohort == "FT26", scenario %in% c("old (as published)", "regime-aware")) %>%
  select(n_ref, scenario, S27_pre, in_col_pre) %>% mutate(across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)

write.csv(scen, file.path(out_dir, "S27_regime_aware_scenarios.csv"), row.names = FALSE)
write.csv(agg, file.path(out_dir, "S27_regime_aware_totals.csv"), row.names = FALSE)
cat("\nDone\n")
