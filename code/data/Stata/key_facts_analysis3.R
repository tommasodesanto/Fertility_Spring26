################################################################################
# ANALYSIS 3: Within-metro vs across-metro moves for new parents
# Uses extract27 which has met2013 + migmet131 (same 2013 CBSA delineation)
################################################################################

library(haven)
library(dplyr)
library(ggplot2)
library(fixest)

outdir <- "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/KeyFacts"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Loading extract27 (selective columns)...\n")

cols_needed <- c("year", "age", "sex", "race", "perwt", "gq",
                 "metro", "met2013", "migmet131",
                 "migrate1",
                 "ownershp", "nchild", "nchlt5", "eldch")

df <- read_dta("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta",
               col_select = all_of(cols_needed))

cat(sprintf("Loaded %d rows\n", nrow(df)))

# Basic cleaning
df <- df %>%
  filter(age >= 22, age <= 45,
         gq %in% c(1, 2))

# Key variables
df <- df %>%
  mutate(
    newparent = as.integer(eldch < 4 & eldch != 99 & nchild > 0),
    moved1y = as.integer(migrate1 >= 2 & !is.na(migrate1)),
    in_principal = case_when(
      metro == 2 ~ 1L,
      metro %in% c(1, 3) ~ 0L,
      TRUE ~ NA_integer_
    )
  )

cat(sprintf("After cleaning: %d rows\n", nrow(df)))

################################################################################
# Within-metro vs across-metro moves
################################################################################

cat("\n=== ANALYSIS 3: Within vs across MSA moves ===\n")

# Focus on movers with identifiable current AND previous MSA
movers <- df %>%
  filter(moved1y == 1,
         met2013 > 0, migmet131 > 0)

movers <- movers %>%
  mutate(
    same_msa = as.integer(met2013 == migmet131),
    diff_msa = as.integer(met2013 != migmet131)
  )

cat(sprintf("Total movers with identifiable current & previous MSA: %d\n", nrow(movers)))

# Quick sanity check on codes
cat(sprintf("Unique current MSAs: %d\n", n_distinct(movers$met2013)))
cat(sprintf("Unique previous MSAs: %d\n", n_distinct(movers$migmet131)))
cat(sprintf("Sample of same-MSA matches: %d\n", sum(movers$same_msa)))

# Overall shares
cat("\n--- All movers: within vs across MSA ---\n")
overall <- movers %>%
  summarise(
    within_share = weighted.mean(same_msa, perwt),
    across_share = weighted.mean(diff_msa, perwt),
    n = n()
  )
cat(sprintf("Within MSA: %.1f%%\n", overall$within_share * 100))
cat(sprintf("Across MSA: %.1f%%\n", overall$across_share * 100))

# By parent status
cat("\n--- By new parent status ---\n")
by_parent <- movers %>%
  group_by(newparent) %>%
  summarise(
    within_share = weighted.mean(same_msa, perwt),
    across_share = weighted.mean(diff_msa, perwt),
    n = n(),
    .groups = "drop"
  )

cat("Non-parents:\n")
cat(sprintf("  Within MSA: %.1f%% | Across MSA: %.1f%% | N = %d\n",
            by_parent$within_share[by_parent$newparent == 0] * 100,
            by_parent$across_share[by_parent$newparent == 0] * 100,
            by_parent$n[by_parent$newparent == 0]))
cat("New parents:\n")
cat(sprintf("  Within MSA: %.1f%% | Across MSA: %.1f%% | N = %d\n",
            by_parent$within_share[by_parent$newparent == 1] * 100,
            by_parent$across_share[by_parent$newparent == 1] * 100,
            by_parent$n[by_parent$newparent == 1]))

# Regression: P(across MSA) ~ newparent
cat("\n--- Regression: P(across-MSA move) on new parent ---\n")
reg_across <- feols(diff_msa ~ newparent + i(sex) | age + year,
                    data = movers, weights = ~perwt)
print(summary(reg_across))

# Among within-MSA movers: principal city status
cat("\n--- Within-MSA movers: where do they end up? ---\n")
within_movers <- movers %>%
  filter(same_msa == 1, !is.na(in_principal))

cat(sprintf("Within-MSA movers with identifiable metro status: %d\n", nrow(within_movers)))

if (nrow(within_movers) > 50) {
  by_parent_pc <- within_movers %>%
    group_by(newparent) %>%
    summarise(
      share_in_principal = weighted.mean(in_principal, perwt),
      n = n(),
      .groups = "drop"
    )
  cat("Non-parents in principal city: ")
  cat(sprintf("%.1f%% (N=%d)\n",
              by_parent_pc$share_in_principal[by_parent_pc$newparent == 0] * 100,
              by_parent_pc$n[by_parent_pc$newparent == 0]))
  cat("New parents in principal city: ")
  cat(sprintf("%.1f%% (N=%d)\n",
              by_parent_pc$share_in_principal[by_parent_pc$newparent == 1] * 100,
              by_parent_pc$n[by_parent_pc$newparent == 1]))

  # Regression
  cat("\n--- Regression: P(in principal city) among within-MSA movers ---\n")
  tryCatch({
    reg_pc <- feols(in_principal ~ newparent + i(sex) | age + year,
                    data = within_movers, weights = ~perwt)
    print(summary(reg_pc))
  }, error = function(e) {
    cat("feols failed, trying OLS:\n")
    reg_pc <- lm(in_principal ~ newparent + factor(sex) + factor(age) + factor(year),
                 data = within_movers, weights = perwt)
    cat(sprintf("newparent coefficient: %.4f (se: %.4f)\n",
                coef(reg_pc)["newparent"],
                summary(reg_pc)$coefficients["newparent", "Std. Error"]))
  })
}

# Also look at ALL movers (not just within-MSA) and principal city status
cat("\n--- ALL movers: principal city status by parent status ---\n")
all_with_pc <- movers %>% filter(!is.na(in_principal))
by_parent_all <- all_with_pc %>%
  group_by(newparent) %>%
  summarise(
    share_in_principal = weighted.mean(in_principal, perwt),
    n = n(),
    .groups = "drop"
  )
cat("Non-parents in principal city: ")
cat(sprintf("%.1f%% (N=%d)\n",
            by_parent_all$share_in_principal[by_parent_all$newparent == 0] * 100,
            by_parent_all$n[by_parent_all$newparent == 0]))
cat("New parents in principal city: ")
cat(sprintf("%.1f%% (N=%d)\n",
            by_parent_all$share_in_principal[by_parent_all$newparent == 1] * 100,
            by_parent_all$n[by_parent_all$newparent == 1]))

# --- Bar chart ---
bar_data <- by_parent %>%
  mutate(parent_label = ifelse(newparent == 1, "New Parents", "Non-Parents")) %>%
  tidyr::pivot_longer(cols = c(within_share, across_share),
                      names_to = "move_type", values_to = "share") %>%
  mutate(move_type = ifelse(move_type == "within_share", "Within MSA", "Across MSA"))

p_bar <- ggplot(bar_data, aes(x = parent_label, y = share, fill = move_type)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", share * 100)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Within MSA" = "navy", "Across MSA" = "firebrick")) +
  labs(x = "", y = "Share of movers", fill = "",
       title = "Move Type by Parent Status",
       subtitle = "ACS movers ages 22-45, matched 2013 CBSA codes") +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")
ggsave(file.path(outdir, "bar_move_type_by_parent.pdf"), p_bar, width = 7, height = 5)
ggsave(file.path(outdir, "bar_move_type_by_parent.png"), p_bar, width = 7, height = 5, dpi = 150)
cat("\nSaved bar_move_type_by_parent\n")

# --- Stacked bar: move type decomposition ---
# Classify moves more finely
movers_classified <- movers %>%
  filter(!is.na(in_principal)) %>%
  mutate(
    move_detail = case_when(
      diff_msa == 1 ~ "Across MSA",
      same_msa == 1 & in_principal == 1 ~ "Within MSA, in principal city",
      same_msa == 1 & in_principal == 0 ~ "Within MSA, not in principal city",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(move_detail))

detail_data <- movers_classified %>%
  group_by(newparent, move_detail) %>%
  summarise(wt = sum(perwt), .groups = "drop") %>%
  group_by(newparent) %>%
  mutate(share = wt / sum(wt),
         parent_label = ifelse(newparent == 1, "New Parents", "Non-Parents"))

cat("\n--- Detailed move classification ---\n")
print(detail_data %>% select(parent_label, move_detail, share))

p_stack <- ggplot(detail_data, aes(x = parent_label, y = share, fill = move_detail)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", share * 100)),
            position = position_stack(vjust = 0.5), size = 3.5, color = "white") +
  scale_fill_manual(values = c("Across MSA" = "firebrick",
                                "Within MSA, in principal city" = "steelblue",
                                "Within MSA, not in principal city" = "navy")) +
  labs(x = "", y = "Share of movers", fill = "",
       title = "Detailed Move Classification by Parent Status",
       subtitle = "ACS movers ages 22-45") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
ggsave(file.path(outdir, "stacked_move_detail_by_parent.pdf"), p_stack, width = 8, height = 6)
ggsave(file.path(outdir, "stacked_move_detail_by_parent.png"), p_stack, width = 8, height = 6, dpi = 150)
cat("Saved stacked_move_detail_by_parent\n")

cat("\n=== DONE ===\n")
