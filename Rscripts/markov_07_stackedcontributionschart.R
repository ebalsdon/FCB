## Stacked contributions chart (A - E)

library(tidyr)

contrib_long <- fcbud_waterfall %>%
  filter(Cohort != "Total") %>%
  filter(!is.na(A) & !is.na(B) & !is.na(C) & !is.na(D) & !is.na(E)) %>%
  mutate(CohortLabel = sprintf("%s (S%d)", Cohort, as.integer(Semester))) %>%
  pivot_longer(A:E, names_to = "Factors", values_to = "Value") %>%
  group_by(CohortLabel) %>%
  mutate(AbsDiff = first(abs(Diff))) %>%
  ungroup()

topK <- 8
keep_rows <- contrib_long %>%
  distinct(CohortLabel, AbsDiff) %>%
  slice_head(n = topK) %>%
  pull(CohortLabel)

plot_df <- contrib_long %>%
  filter(CohortLabel %in% keep_rows) %>%
  mutate(
    Factors = factor(Factors, levels = c("A","B","C","D","E"),
                       labels = c("Admissions","% FCB","Retention","Net transfers","Promotion")),
    CohortLabel = fct_reorder(CohortLabel, AbsDiff, .desc = TRUE)
  )

ggplot(plot_df, aes(x = CohortLabel, y = Value, fill = Factors)) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
  scale_y_continuous(labels = comma) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = NULL, y = "Headcount contribution",
    title = "Decomposition of cohort contributions to change F23 → F24",
    subtitle = "Positive (above zero) and negative (below zero)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "top"
  )



library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(scales)

# --- 1) Build a tidy long table from your waterfall (A–E columns) ---
# Assumes fcbud_waterfall has: Cohort, Semester, Diff, A,B,C,D,E, Pct_of_Total
# Pick the same rows you plotted earlier (e.g., Top-K by |Diff|)
top_k <- 8
wf_long <- fcbud_waterfall %>%
  filter(Cohort != "Total") %>%
  arrange(desc(abs(Diff))) %>%
  slice_head(n = top_k) %>%
  mutate(Label = paste0(Cohort, " (S", Semester, ")")) %>%
  select(Label, Cohort, Semester, Diff, `Admissions`=A, `% FCB`=B, `Retention`=C, `Net transfers`=D, `Promotion`=E) %>%
  pivot_longer(cols = c(`Admissions`,`% FCB`,`Retention`,`Net transfers`,`Promotion`),
               names_to = "Factor", values_to = "Contribution") %>%
  # Nice hover text
  mutate(hover = paste0(
    Label, "<br>",
    Factor, ": ", comma(round(Contribution)), "<br>",
    "Total contribution (Δ): ", comma(round(Diff))
  ))

# order x by total absolute contribution, newest on left if you prefer
wf_long <- wf_long %>%
  group_by(Label) %>% mutate(abs_tot = first(abs(Diff))) %>% ungroup() %>%
  mutate(Label = factor(Label, levels = unique(Label[order(-abs_tot)])))

# --- 2) ggplot (your existing stacked chart), just add aes(text=hover) ---
p_stack <- ggplot(wf_long, aes(x = Label, y = Contribution, fill = Factor, text = hover)) +
  geom_col(position = "stack", width = 0.75) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Decomposition of cohort contributions to change F23 → F24",
    subtitle = "Positive (above zero) and negative (below zero)",
    x = NULL, y = "Headcount contribution", fill = "Factors"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 28, hjust = 1)
  )

# --- 3) Make it interactive ---
ggplotly(p_stack, tooltip = "text")



