################################################################################
# MOVE DECOMPOSITION: Cross the MSA margin with the principal city margin
# For movers: 4 categories × parent status
################################################################################

library(haven)
library(dplyr)
library(ggplot2)
library(fixest)
library(tidyr)

outdir <- "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/KeyFacts"

cat("Loading extract27...\n")
df <- read_dta(
  "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta",
  col_select = c("year","age","sex","race","perwt","gq",
                 "metro","met2013","migmet131",
                 "migrate1","migrate1d","nchild","nchlt5","eldch")
)

df <- df %>%
  filter(age >= 22, age <= 45, gq %in% c(1,2)) %>%
  mutate(
    newparent = as.integer(eldch < 4 & eldch != 99 & nchild > 0),
    moved1y = as.integer(migrate1 >= 2 & !is.na(migrate1)),
    stayer = as.integer(migrate1 == 1),
    in_principal = case_when(
      metro == 2 ~ 1L, metro %in% c(1,3) ~ 0L, TRUE ~ NA_integer_
    ),
    in_metro = as.integer(metro %in% c(2,3))
  )

cat(sprintf("Sample: %d\n", nrow(df)))


################################################################################
# 1. Full population: analogous to the principal city table but for MSA
################################################################################

cat("\n=== PART 1: Full population outcomes ===\n")

# For the full population, define outcomes parallel to the principal city table:
# (a) Currently in principal city
# (b) Moved AND now in principal city ("moved to principal city")
# (c) Moved AND now NOT in principal city ("moved to suburb/non-metro")
# (d) Changed MSA in last year
# (e) Moved within same MSA

full <- df %>%
  filter(!is.na(in_principal)) %>%
  mutate(
    has_msa_info = as.integer(met2013 > 0 & (stayer == 1 | migmet131 > 0)),
    changed_msa = case_when(
      stayer == 1 ~ 0L,
      moved1y == 1 & met2013 > 0 & migmet131 > 0 & met2013 != migmet131 ~ 1L,
      moved1y == 1 & met2013 > 0 & migmet131 > 0 & met2013 == migmet131 ~ 0L,
      TRUE ~ NA_integer_
    ),
    moved_within_msa = case_when(
      stayer == 1 ~ 0L,
      moved1y == 1 & met2013 > 0 & migmet131 > 0 & met2013 == migmet131 ~ 1L,
      moved1y == 1 & met2013 > 0 & migmet131 > 0 & met2013 != migmet131 ~ 0L,
      TRUE ~ NA_integer_
    )
  )

# Table: means by parent status (like the principal city table in slides)
cat("\n--- Weighted means ---\n")
means <- full %>%
  filter(!is.na(changed_msa)) %>%
  group_by(newparent) %>%
  summarise(
    `Lives in principal city` = weighted.mean(in_principal, perwt),
    `Moved (any)` = weighted.mean(moved1y, perwt),
    `Changed MSA` = weighted.mean(changed_msa, perwt),
    `Moved within MSA` = weighted.mean(moved_within_msa, perwt),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = ifelse(newparent == 1, "New Parents", "Non-Parents"))

cat("\n")
print(means %>% select(label, everything(), -newparent, -n))

# Regressions
cat("\n--- Regressions: outcome ~ newparent | age + year ---\n")
reg_data <- full %>% filter(!is.na(changed_msa))

outcomes <- c("in_principal", "moved1y", "changed_msa", "moved_within_msa")
labels <- c("Lives in principal city", "Moved (any)", "Changed MSA", "Moved within MSA")

for (i in seq_along(outcomes)) {
  f <- as.formula(paste(outcomes[i], "~ newparent + i(sex) + i(race) | age + year"))
  r <- feols(f, data = reg_data, weights = ~perwt)
  cat(sprintf("  %s: newparent = %+.4f (se %.4f, p = %.4f)\n",
              labels[i], coef(r)["newparent"], se(r)["newparent"],
              pvalue(r)["newparent"]))
}


################################################################################
# 2. Among MOVERS only: 4-way decomposition
################################################################################

cat("\n\n=== PART 2: Movers — 4-way decomposition ===\n")

movers <- df %>%
  filter(moved1y == 1, met2013 > 0, migmet131 > 0, !is.na(in_principal)) %>%
  mutate(
    same_msa = as.integer(met2013 == migmet131),
    diff_msa = as.integer(met2013 != migmet131),
    # 4 categories
    cat4 = case_when(
      same_msa == 1 & in_principal == 1 ~ "Within MSA, principal city",
      same_msa == 1 & in_principal == 0 ~ "Within MSA, suburb",
      diff_msa == 1 & in_principal == 1 ~ "Across MSA, principal city",
      diff_msa == 1 & in_principal == 0 ~ "Across MSA, suburb",
      TRUE ~ NA_character_
    )
  )

# Shares by parent status
cat("\n--- Share of movers in each category ---\n")
shares <- movers %>%
  group_by(newparent, cat4) %>%
  summarise(wt = sum(perwt), .groups = "drop") %>%
  group_by(newparent) %>%
  mutate(share = wt / sum(wt)) %>%
  ungroup() %>%
  mutate(label = ifelse(newparent == 1, "New Parents", "Non-Parents"))

shares_wide <- shares %>%
  select(label, cat4, share) %>%
  pivot_wider(names_from = label, values_from = share)

print(shares_wide)

# Differences
cat("\n--- Difference (New Parents - Non-Parents) ---\n")
diff_table <- shares %>%
  select(newparent, cat4, share) %>%
  pivot_wider(names_from = newparent, values_from = share) %>%
  mutate(diff = `1` - `0`)
print(diff_table)

# Regressions for each category indicator
cat("\n--- Regressions: P(category) ~ newparent | age + year ---\n")
movers <- movers %>%
  mutate(
    within_principal = as.integer(cat4 == "Within MSA, principal city"),
    within_suburb    = as.integer(cat4 == "Within MSA, suburb"),
    across_principal = as.integer(cat4 == "Across MSA, principal city"),
    across_suburb    = as.integer(cat4 == "Across MSA, suburb")
  )

cat_outcomes <- c("within_principal", "within_suburb", "across_principal", "across_suburb")
cat_labels <- c("Within MSA → principal", "Within MSA → suburb",
                "Across MSA → principal", "Across MSA → suburb")

regs <- list()
for (i in seq_along(cat_outcomes)) {
  f <- as.formula(paste(cat_outcomes[i], "~ newparent + i(sex) + i(race) | age + year"))
  r <- feols(f, data = movers, weights = ~perwt)
  regs[[i]] <- r
  cat(sprintf("  %s: newparent = %+.4f (se %.4f, p = %.4f)\n",
              cat_labels[i], coef(r)["newparent"], se(r)["newparent"],
              pvalue(r)["newparent"]))
}

# Plot: stacked bar with 4 categories
p1 <- ggplot(shares, aes(x = label, y = share, fill = cat4)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", share * 100)),
            position = position_stack(vjust = 0.5), size = 3.2, color = "white") +
  scale_fill_manual(values = c(
    "Within MSA, principal city" = "steelblue",
    "Within MSA, suburb" = "navy",
    "Across MSA, principal city" = "lightsalmon",
    "Across MSA, suburb" = "firebrick"
  )) +
  labs(x = "", y = "Share of movers", fill = "",
       title = "Where Do Movers Go? By Parent Status",
       subtitle = "ACS movers ages 22-45, 4-way decomposition") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9))
ggsave(file.path(outdir, "stacked_4way_decomposition.pdf"), p1, width = 8, height = 6)
ggsave(file.path(outdir, "stacked_4way_decomposition.png"), p1, width = 8, height = 6, dpi = 150)
cat("\nSaved stacked_4way_decomposition\n")

# Grouped bar (easier to read differences)
p2 <- ggplot(shares, aes(x = cat4, y = share, fill = label)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", share * 100)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("New Parents" = "firebrick", "Non-Parents" = "navy")) +
  labs(x = "", y = "Share of movers", fill = "",
       title = "Move Destination: New Parents vs Non-Parents",
       subtitle = "ACS movers ages 22-45") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 9, angle = 15, hjust = 1))
ggsave(file.path(outdir, "grouped_4way_decomposition.pdf"), p2, width = 9, height = 5.5)
ggsave(file.path(outdir, "grouped_4way_decomposition.png"), p2, width = 9, height = 5.5, dpi = 150)
cat("Saved grouped_4way_decomposition\n")


################################################################################
# 3. Coefficient plot: newparent effect on each destination type
################################################################################

coef_data <- data.frame(
  outcome = cat_labels,
  estimate = sapply(regs, function(r) coef(r)["newparent"]),
  se = sapply(regs, function(r) se(r)["newparent"])
) %>%
  mutate(
    ci_lo = estimate - 1.96 * se,
    ci_hi = estimate + 1.96 * se,
    outcome = factor(outcome, levels = rev(cat_labels))
  )

p3 <- ggplot(coef_data, aes(x = estimate, y = outcome)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2, linewidth = 0.8) +
  geom_point(size = 3, color = "navy") +
  labs(x = "Effect of new parenthood (pp)", y = "",
       title = "Where Do New Parents Move?",
       subtitle = "Coefficients from P(destination) ~ new parent | age + year + sex + race FE") +
  theme_minimal(base_size = 13)
ggsave(file.path(outdir, "coefplot_4way.pdf"), p3, width = 8, height = 4)
ggsave(file.path(outdir, "coefplot_4way.png"), p3, width = 8, height = 4, dpi = 150)
cat("Saved coefplot_4way\n")


################################################################################
# Summary
################################################################################

cat("\n\n")
cat("################################################################\n")
cat("#           MOVE DECOMPOSITION: SUMMARY                        #\n")
cat("################################################################\n\n")

cat("Full population (analogous to your principal city table):\n")
print(means %>% select(label, `Lives in principal city`, `Moved (any)`,
                         `Changed MSA`, `Moved within MSA`))

cat("\n4-way decomposition of movers:\n")
print(shares_wide)

cat("\nNew parent coefficients (from regression):\n")
for (i in seq_along(cat_labels)) {
  cat(sprintf("  %s: %+.3f pp\n", cat_labels[i], coef(regs[[i]])["newparent"]*100))
}

cat("\n=== DONE ===\n")
