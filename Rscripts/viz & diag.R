

# --- Generate Sankey diagram and alluvial plot for a specific transition matrix ---
# Convert matrix to long format
library(tibble)
library(dplyr)

mat <- transition_matrices[["FT19_4_5"]]
df_long <- as.data.frame(as.table(mat)) %>%
  rename(from = Var1, to = Var2, value = Freq) %>%
  filter(value > 0)

# Create node list and index links
nodes <- tibble(name = unique(c(df_long$from, df_long$to)))
df_long <- df_long %>%
  mutate(
    source = match(from, nodes$name) - 1,
    target = match(to, nodes$name) - 1
  )

# Create Sankey
sankeyNetwork(
  Links = df_long,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  fontSize = 12,
  nodeWidth = 30
)

# Create alluvial plot
mat <- transition_matrices[["FT19_4_5"]]
df_alluvial <- as.data.frame(as.table(mat)) %>%
  rename(from = Var1, to = Var2, value = Freq) %>%
  filter(value > 0)

ggplot(df_alluvial,
       aes(axis1 = from, axis2 = to, y = value)) +
  geom_alluvium(aes(fill = from), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Term 4", "Term 5"), expand = c(.05, .05)) +
  theme_minimal() +
  labs(title = "FT19: Transition from Term 4 to 5")

# --- Same viz tools, collapsing "Other Colleges" ---
collapse_colleges <- function(label) {
  if (label %in% c("AL", "EDU", "ENG", "HHS", "PSFA", "SCI")) {
    return("Other Colleges")
  } else {
    return(label)
  }
}

df_long <- as.data.frame(as.table(transition_matrices[["FT19_4_5"]])) %>%
  rename(from = Var1, to = Var2, value = Freq) %>%
  filter(value > 0) %>%
  mutate(
    from = sapply(from, collapse_colleges),
    to = sapply(to, collapse_colleges)
  ) %>%
  group_by(from, to) %>%
  summarize(value = sum(value), .groups = "drop")

nodes <- tibble(name = unique(c(df_long$from, df_long$to)))
df_long <- df_long %>%
  mutate(
    source = match(from, nodes$name) - 1,
    target = match(to, nodes$name) - 1
  )

sankeyNetwork(
  Links = df_long,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  fontSize = 12,
  nodeWidth = 30
)

df_alluvial <- as.data.frame(as.table(transition_matrices[["FT19_4_5"]])) %>%
  rename(from = Var1, to = Var2, value = Freq) %>%
  filter(value > 0) %>%
  mutate(
    from = sapply(from, collapse_colleges),
    to = sapply(to, collapse_colleges)
  ) %>%
  group_by(from, to) %>%
  summarize(value = sum(value), .groups = "drop")

ggplot(df_alluvial,
       aes(axis1 = from, axis2 = to, y = value)) +
  geom_alluvium(aes(fill = from), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Term 4", "Term 5"), expand = c(.05, .05)) +
  theme_minimal() +
  labs(title = "FT19: Transition from Term 4 to 5 (Simplified)")

display_matrix <- addmargins(transition_matrices[["FT21_4_5"]])
View(display_matrix)

check_transition_matrices <- function(matrices, expected_dim = 10, ne_label = "NE") {
  tibble::tibble(
    name = names(matrices),
    has_na     = sapply(matrices, function(m) any(is.na(m))),
    all_zero   = sapply(matrices, function(m) all(m == 0)),
    wrong_dim  = sapply(matrices, function(m) !all(dim(m) == c(expected_dim, expected_dim))),
    total_sum  = sapply(matrices, sum),
    sum_excl_NE = sapply(matrices, function(m) {
      if (ne_label %in% colnames(m)) {
        sum(m[, setdiff(colnames(m), ne_label)])
      } else {
        NA_real_
      }
    })
  ) %>%
    arrange(desc(all_zero), has_na, wrong_dim)
}

x <- check_transition_matrices(transition_matrices)


which_target_s25 <- sapply(names(transition_matrices), function(name) {
  parts <- strsplit(name, "_")[[1]]
  cohort <- parts[1]
  to_sem <- as.integer(parts[3])
  start_index <- match(term_starts[[cohort]], terms)
  term_index <- start_index + to_sem - 2
  terms[term_index] == "S25"
})

names(transition_matrices)[which_target_s25]
