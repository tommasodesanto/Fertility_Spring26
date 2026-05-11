################################################################################
# KEY FACTS: Two analyses to guide model design
#
# Analysis 1: MSA-level fertility vs ownership/price-to-income
# Analysis 3: Within-metro vs across-metro moves for new parents
################################################################################

library(haven)
library(dplyr)
library(ggplot2)
library(fixest)

outdir <- "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/KeyFacts"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Loading data (selective columns)...\n")

# Only read the columns we need to save memory
cols_needed <- c("year", "age", "sex", "race", "perwt", "gq",
                 "metro", "metarea", "met2013", "migmet1",
                 "migrate1", "migrate1d",
                 "ownershp", "nchild", "nchlt5", "eldch", "fertyr",
                 "valueh", "rent", "hhincome", "rooms")

df <- read_dta("/Users/tommasodesanto/Desktop/Projects/Datasets/acs_with_hybrid_codes_andHPI.dta",
               col_select = all_of(cols_needed))

cat(sprintf("Loaded %d rows\n", nrow(df)))

# Basic cleaning
df <- df %>%
  filter(age >= 22, age <= 45,
         gq %in% c(1, 2)) # Households only

# Key variables
df <- df %>%
  mutate(
    owner = as.integer(ownershp == 1),
    has_youngchild = as.integer(nchlt5 > 0 & !is.na(nchlt5)),
    newparent = as.integer(eldch < 4 & eldch != 99 & nchild > 0),
    moved1y = as.integer(migrate1 >= 2 & !is.na(migrate1)),
    # Metro status
    in_principal = case_when(
      metro == 2 ~ 1L,
      metro %in% c(1, 3) ~ 0L,
      TRUE ~ NA_integer_
    ),
    in_metro = as.integer(metro %in% c(2, 3)),
    # Price to income (owners only)
    pti = ifelse(owner == 1 & valueh > 0 & valueh < 9999999 &
                   hhincome > 0 & hhincome < 9999999,
                 valueh / hhincome, NA_real_),
    # Rent to income
    rti = ifelse(rent > 0 & !is.na(rent) & hhincome > 0 & hhincome < 9999999,
                 rent * 12 / hhincome, NA_real_)
  )

cat(sprintf("After cleaning: %d rows\n", nrow(df)))

################################################################################
# ANALYSIS 1: MSA-level fertility vs ownership vs housing costs
################################################################################

cat("\n=== ANALYSIS 1: MSA-level fertility vs affordability ===\n")

# Collapse to MSA-year
msa_yr <- df %>%
  filter(met2013 > 0) %>%
  group_by(met2013, year) %>%
  summarise(
    fert_rate = weighted.mean(has_youngchild, perwt, na.rm = TRUE),
    own_rate = weighted.mean(owner, perwt, na.rm = TRUE),
    med_pti = median(pti, na.rm = TRUE),
    med_rti = median(rti, na.rm = TRUE),
    mean_rent = weighted.mean(rent, perwt, na.rm = TRUE),
    mean_hhinc = weighted.mean(hhincome, perwt, na.rm = TRUE),
    pop = sum(perwt),
    .groups = "drop"
  ) %>%
  filter(pop >= 50000)

# Average across years for stable MSA estimates
msa <- msa_yr %>%
  group_by(met2013) %>%
  summarise(
    fert_rate = weighted.mean(fert_rate, pop),
    own_rate = weighted.mean(own_rate, pop),
    med_pti = median(med_pti, na.rm = TRUE),
    med_rti = median(med_rti, na.rm = TRUE),
    pop = mean(pop),
    .groups = "drop"
  )

cat(sprintf("Number of MSAs: %d\n", nrow(msa)))
cat("\nSummary statistics:\n")
print(summary(msa[, c("fert_rate", "own_rate", "med_pti", "med_rti")]))

# Correlations
cat("\nWeighted correlations:\n")
cat(sprintf("  Fertility vs Ownership:  r = %.3f\n",
            cov.wt(msa[, c("fert_rate", "own_rate")], wt = msa$pop, cor = TRUE)$cor[1,2]))
cor_pti <- msa %>% filter(!is.na(med_pti))
cat(sprintf("  Fertility vs PTI:        r = %.3f\n",
            cov.wt(cor_pti[, c("fert_rate", "med_pti")], wt = cor_pti$pop, cor = TRUE)$cor[1,2]))

# --- Scatter 1: Fertility vs Ownership ---
p1 <- ggplot(msa, aes(x = own_rate, y = fert_rate, size = pop)) +
  geom_point(alpha = 0.4, color = "navy") +
  geom_smooth(method = "lm", aes(weight = pop), color = "firebrick",
              se = TRUE, linewidth = 1) +
  labs(x = "Homeownership rate",
       y = "Fertility rate (share with child < 5)",
       title = "Fertility vs Ownership across MSAs",
       subtitle = "Ages 22-45, ACS pooled") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
ggsave(file.path(outdir, "scatter_fertility_ownership.pdf"), p1, width = 8, height = 6)
ggsave(file.path(outdir, "scatter_fertility_ownership.png"), p1, width = 8, height = 6, dpi = 150)
cat("Saved scatter_fertility_ownership\n")

# --- Scatter 2: Fertility vs Price-to-income ---
p2 <- ggplot(msa %>% filter(!is.na(med_pti)), aes(x = med_pti, y = fert_rate, size = pop)) +
  geom_point(alpha = 0.4, color = "navy") +
  geom_smooth(method = "lm", aes(weight = pop), color = "firebrick",
              se = TRUE, linewidth = 1) +
  labs(x = "Median price-to-income ratio",
       y = "Fertility rate (share with child < 5)",
       title = "Fertility vs Housing Costs across MSAs",
       subtitle = "Ages 22-45, ACS pooled") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
ggsave(file.path(outdir, "scatter_fertility_pricetoincome.pdf"), p2, width = 8, height = 6)
ggsave(file.path(outdir, "scatter_fertility_pricetoincome.png"), p2, width = 8, height = 6, dpi = 150)
cat("Saved scatter_fertility_pricetoincome\n")

# --- Scatter 3: Fertility vs Rent-to-income ---
p3 <- ggplot(msa %>% filter(!is.na(med_rti)), aes(x = med_rti, y = fert_rate, size = pop)) +
  geom_point(alpha = 0.4, color = "navy") +
  geom_smooth(method = "lm", aes(weight = pop), color = "firebrick",
              se = TRUE, linewidth = 1) +
  labs(x = "Median rent-to-income ratio",
       y = "Fertility rate (share with child < 5)",
       title = "Fertility vs Rental Costs across MSAs",
       subtitle = "Ages 22-45, ACS pooled") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
ggsave(file.path(outdir, "scatter_fertility_renttoincome.pdf"), p3, width = 8, height = 6)
ggsave(file.path(outdir, "scatter_fertility_renttoincome.png"), p3, width = 8, height = 6, dpi = 150)
cat("Saved scatter_fertility_renttoincome\n")

# --- Scatter 4: Ownership vs Price-to-income ---
p4 <- ggplot(msa %>% filter(!is.na(med_pti)), aes(x = med_pti, y = own_rate, size = pop)) +
  geom_point(alpha = 0.4, color = "navy") +
  geom_smooth(method = "lm", aes(weight = pop), color = "firebrick",
              se = TRUE, linewidth = 1) +
  labs(x = "Median price-to-income ratio",
       y = "Homeownership rate",
       title = "Ownership vs Housing Costs across MSAs",
       subtitle = "Ages 22-45, ACS pooled") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
ggsave(file.path(outdir, "scatter_ownership_pricetoincome.pdf"), p4, width = 8, height = 6)
ggsave(file.path(outdir, "scatter_ownership_pricetoincome.png"), p4, width = 8, height = 6, dpi = 150)
cat("Saved scatter_ownership_pricetoincome\n")

# Weighted regressions
cat("\n--- MSA-level regressions (weighted by pop) ---\n")
cat("\nFertility ~ Ownership:\n")
print(summary(lm(fert_rate ~ own_rate, data = msa, weights = pop)))
cat("\nFertility ~ PTI:\n")
print(summary(lm(fert_rate ~ med_pti, data = msa %>% filter(!is.na(med_pti)), weights = pop)))
cat("\nFertility ~ Ownership + PTI:\n")
print(summary(lm(fert_rate ~ own_rate + med_pti, data = msa %>% filter(!is.na(med_pti)), weights = pop)))


################################################################################
# ANALYSIS 3: Within-metro vs across-metro moves
################################################################################

cat("\n\n=== ANALYSIS 3: Within vs across MSA moves ===\n")

movers <- df %>%
  filter(moved1y == 1,
         metarea > 0, migmet1 > 0,
         !is.na(metarea), !is.na(migmet1))

movers <- movers %>%
  mutate(
    same_msa = as.integer(metarea == migmet1),
    diff_msa = as.integer(metarea != migmet1)
  )

cat(sprintf("Total movers with identifiable MSAs: %d\n", nrow(movers)))

# Overall shares
cat("\n--- All movers: within vs across MSA ---\n")
overall <- movers %>%
  summarise(
    within_share = weighted.mean(same_msa, perwt),
    across_share = weighted.mean(diff_msa, perwt),
    n = n()
  )
print(overall)

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
print(by_parent)

cat(sprintf("\nDifference in across-MSA share: new parents %.1f%% vs non-parents %.1f%%\n",
            by_parent$across_share[by_parent$newparent == 1] * 100,
            by_parent$across_share[by_parent$newparent == 0] * 100))

# Regression: P(across MSA) ~ newparent
cat("\n--- Regression: P(across MSA move) on new parent ---\n")
reg_across <- feols(diff_msa ~ newparent + i(sex) + i(race) | age + year,
                    data = movers, weights = ~perwt)
print(summary(reg_across))

# Among within-MSA movers: principal city status
cat("\n--- Within-MSA movers: principal city status by parent status ---\n")
within_movers <- movers %>%
  filter(same_msa == 1, !is.na(in_principal))

by_parent_pc <- within_movers %>%
  group_by(newparent) %>%
  summarise(
    share_in_principal = weighted.mean(in_principal, perwt),
    n = n(),
    .groups = "drop"
  )
print(by_parent_pc)

cat(sprintf("\nNew parents in principal city: %.1f%% vs non-parents: %.1f%%\n",
            by_parent_pc$share_in_principal[by_parent_pc$newparent == 1] * 100,
            by_parent_pc$share_in_principal[by_parent_pc$newparent == 0] * 100))

# Regression: P(in principal city) ~ newparent, within-MSA movers
cat("\n--- Regression: P(in principal city) among within-MSA movers ---\n")
reg_pc <- feols(in_principal ~ newparent + i(sex) + i(race) | age + year,
                data = within_movers, weights = ~perwt)
print(summary(reg_pc))

# --- Bar chart ---
bar_data <- movers %>%
  group_by(newparent) %>%
  summarise(
    `Within MSA` = weighted.mean(same_msa, perwt),
    `Across MSA` = weighted.mean(diff_msa, perwt),
    .groups = "drop"
  ) %>%
  mutate(parent_label = ifelse(newparent == 1, "New Parents", "Non-Parents")) %>%
  tidyr::pivot_longer(cols = c(`Within MSA`, `Across MSA`),
                      names_to = "move_type", values_to = "share")

p_bar <- ggplot(bar_data, aes(x = parent_label, y = share, fill = move_type)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("Within MSA" = "navy", "Across MSA" = "firebrick")) +
  labs(x = "", y = "Share of movers", fill = "",
       title = "Move Type by Parent Status",
       subtitle = "ACS movers ages 22-45, identifiable MSAs") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")
ggsave(file.path(outdir, "bar_move_type_by_parent.pdf"), p_bar, width = 7, height = 5)
ggsave(file.path(outdir, "bar_move_type_by_parent.png"), p_bar, width = 7, height = 5, dpi = 150)
cat("Saved bar_move_type_by_parent\n")


################################################################################
# BONUS: Age split for Analysis 1
################################################################################

cat("\n\n=== BONUS: Age-split fertility-ownership gradient ===\n")

msa_age <- df %>%
  filter(met2013 > 0) %>%
  mutate(age_group = ifelse(age <= 32, "22-32", "33-45")) %>%
  group_by(met2013, age_group) %>%
  summarise(
    fert_rate = weighted.mean(has_youngchild, perwt, na.rm = TRUE),
    own_rate = weighted.mean(owner, perwt, na.rm = TRUE),
    pop = sum(perwt),
    .groups = "drop"
  ) %>%
  filter(pop >= 20000)

p_age <- ggplot(msa_age, aes(x = own_rate, y = fert_rate, size = pop, color = age_group)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", aes(weight = pop), se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = c("22-32" = "navy", "33-45" = "firebrick")) +
  labs(x = "Homeownership rate",
       y = "Fertility rate (share with child < 5)",
       title = "Fertility vs Ownership by Age Group",
       color = "Age group") +
  theme_minimal(base_size = 14) +
  guides(size = "none")
ggsave(file.path(outdir, "scatter_fertility_ownership_byage.pdf"), p_age, width = 9, height = 6)
ggsave(file.path(outdir, "scatter_fertility_ownership_byage.png"), p_age, width = 9, height = 6, dpi = 150)

cat("\nYoung (22-32):\n")
print(summary(lm(fert_rate ~ own_rate, data = msa_age %>% filter(age_group == "22-32"), weights = pop)))
cat("\nOlder (33-45):\n")
print(summary(lm(fert_rate ~ own_rate, data = msa_age %>% filter(age_group == "33-45"), weights = pop)))

cat("\n=== DONE ===\n")
cat(sprintf("Outputs saved to: %s\n", outdir))
