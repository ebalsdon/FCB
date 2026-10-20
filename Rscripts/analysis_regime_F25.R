# Diagnostics: pre-/post-F25 regime shift in switching into FCBpre (30-unit + GPA gate).
# Calendar step = (from term -> to term); stage = semester step s -> s+1 of the cohort.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(tibble) })
options(width = 220, dplyr.width = Inf, pillar.width = Inf)
e <- new.env(); load(here::here("_outputs/Rdata/data_image1.RData"), envir = e)
TM <- e$transition_matrices; terms <- e$terms; term_starts <- e$term_starts
states <- colnames(TM[[1]]); colleges <- setdiff(states, c("FCBud", "FCBpre", "NE"))
out_dir <- here::here("_outputs/analysis"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

keys <- grep("^(FT|TR)[0-9]+_[0-9]+_[0-9]+$", names(TM), value = TRUE)
meta <- bind_rows(lapply(keys, function(k) {
  p <- strsplit(k, "_")[[1]]; c <- p[1]; a <- as.integer(p[2]); b <- as.integer(p[3])
  if (b != a + 1L) return(NULL)
  i <- match(term_starts[[c]], terms) + a - 1L
  if (is.na(i) || i + 1L > length(terms)) return(NULL)
  tibble(key = k, cohort = c, grp = substr(c, 1, 2), stage = a, from = terms[i], to = terms[i + 1L])
}))
meta <- meta %>% mutate(step = paste0(from, ">", to),
  regime = case_when(to %in% c("F25") ~ "boundary (S25>F25)", match(from, terms) >= match("F25", terms) ~ "post", TRUE ~ "pre"))

row_of <- function(key) {
  m <- TM[[key]][states, states]
  tibble(n_pre = sum(m["FCBpre", ]), n_ud = sum(m["FCBud", ]),
         n_col = sum(m[colleges, ]), n_ne = sum(m["NE", ]),
         col_to_pre = sum(m[colleges, "FCBpre"]), col_to_ud = sum(m[colleges, "FCBud"]),
         ne_to_pre = m["NE", "FCBpre"], ne_to_ud = m["NE", "FCBud"],
         pre_stay = m["FCBpre", "FCBpre"], pre_to_ud = m["FCBpre", "FCBud"],
         pre_to_col = sum(m["FCBpre", colleges]), pre_to_ne = m["FCBpre", "NE"])
}
d <- bind_cols(meta, bind_rows(lapply(meta$key, row_of)))
d <- d %>% mutate(stage_bin = case_when(stage <= 1 ~ "1>2", stage == 2 ~ "2>3", stage == 3 ~ "3>4", stage == 4 ~ "4>5", TRUE ~ "5+"))
write.csv(d, file.path(out_dir, "regime_F25_matrix_level.csv"), row.names = FALSE)

cat("\n==== A. Switching INTO FCBpre by calendar step (FT cohorts, spring vs fall) ====\n")
a <- d %>% filter(grp == "FT") %>% group_by(step) %>% summarise(regime = first(regime), col_n = sum(n_col), col_to_pre = sum(col_to_pre),
  rate_per_1000 = 1000 * sum(col_to_pre) / sum(n_col), ne_to_pre = sum(ne_to_pre), .groups = "drop") %>% arrange(match(substr(step, 1, 3), terms))
print(as.data.frame(a %>% mutate(rate_per_1000 = round(rate_per_1000, 2))), row.names = FALSE)

cat("\n==== B. College->FCBpre by stage, SPRING steps (S: from F -> to S) FT only; count / rate per 1000 college students ====\n")
b <- d %>% filter(grp == "FT", startsWith(from, "F"), startsWith(to, "S")) %>% group_by(step, regime, stage_bin) %>%
  summarise(n_col = sum(n_col), n = sum(col_to_pre), per1000 = round(1000 * n / n_col, 2), .groups = "drop")
print(as.data.frame(b %>% select(-n_col) %>% pivot_wider(names_from = step, values_from = c(n, per1000))), row.names = FALSE)

cat("\n==== C. College->FCBpre by stage, FALL steps (from S -> to F) FT only ====\n")
cc <- d %>% filter(grp == "FT", startsWith(from, "S"), startsWith(to, "F")) %>% group_by(step, regime, stage_bin) %>%
  summarise(n_col = sum(n_col), n = sum(col_to_pre), per1000 = round(1000 * n / n_col, 2), .groups = "drop")
print(as.data.frame(cc %>% select(-n_col) %>% pivot_wider(names_from = step, values_from = c(n, per1000))), row.names = FALSE)

cat("\n==== D. Source college of switchers into FCBpre, FT, by calendar step ====\n")
src <- bind_rows(lapply(seq_len(nrow(d)), function(i) { if (d$grp[i] != "FT") return(NULL); m <- TM[[d$key[i]]][states, states]
  tibble(step = d$step[i], college = colleges, n = m[colleges, "FCBpre"], base = rowSums(m[colleges, , drop = FALSE])) })) %>%
  group_by(step, college) %>% summarise(n = sum(n), base = sum(base), .groups = "drop")
print(as.data.frame(src %>% select(step, college, n) %>% pivot_wider(names_from = step, values_from = n, values_fill = 0)), row.names = FALSE)

cat("\n==== E. NE -> FCBpre/UD by calendar step (re-entry / unclassified inflow, all cohorts) ====\n")
ee <- d %>% group_by(step) %>% summarise(regime = first(regime), ne_to_pre = sum(ne_to_pre), ne_to_ud = sum(ne_to_ud), .groups = "drop") %>%
  arrange(match(substr(step, 1, 3), terms)); print(as.data.frame(ee), row.names = FALSE)

cat("\n==== F. Pre-major progression by stage and calendar step, FT (pre->UD, pre stay, pre->NE rates) ====\n")
f <- d %>% filter(grp == "FT", stage <= 7) %>% group_by(stage, step) %>% summarise(regime = first(regime), n = sum(n_pre),
  promote = sum(pre_to_ud) / sum(n_pre), stay = sum(pre_stay) / sum(n_pre), to_ne = sum(pre_to_ne) / sum(n_pre), .groups = "drop") %>%
  filter(n > 0) %>% arrange(stage, match(substr(step, 1, 3), terms))
print(as.data.frame(f %>% mutate(across(c(promote, stay, to_ne), ~ round(.x, 3)))), row.names = FALSE)
write.csv(f, file.path(out_dir, "regime_F25_pre_progression_by_stage.csv"), row.names = FALSE)
cat("\nDone\n")
