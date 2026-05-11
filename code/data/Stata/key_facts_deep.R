################################################################################
# DEEP DIVE: Tenure × Move × Parenthood interactions
# All output to Outputs/KeyFacts/ — nothing else touched
################################################################################

library(haven)
library(dplyr)
library(ggplot2)
library(fixest)
library(tidyr)

outdir <- "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/KeyFacts"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# PART A: Load extract27 for move × tenure analysis
# ============================================================================

cat("=== Loading extract27 ===\n")
df <- read_dta(
  "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta",
  col_select = c("year","age","sex","race","perwt","gq",
                 "metro","met2013","migmet131",
                 "migrate1","ownershp","nchild","nchlt5","eldch",
                 "valueh","rent","hhincome")
)
cat(sprintf("Loaded %d rows\n", nrow(df)))

df <- df %>%
  filter(age >= 22, age <= 45, gq %in% c(1,2)) %>%
  mutate(
    owner    = as.integer(ownershp == 1),
    renter   = as.integer(ownershp == 2),
    newparent = as.integer(eldch < 4 & eldch != 99 & nchild > 0),
    has_youngchild = as.integer(nchlt5 > 0 & !is.na(nchlt5)),
    moved1y  = as.integer(migrate1 >= 2 & !is.na(migrate1)),
    in_principal = case_when(
      metro == 2 ~ 1L, metro %in% c(1,3) ~ 0L, TRUE ~ NA_integer_
    ),
    pti = ifelse(owner==1 & valueh>0 & valueh<9999999 &
                   hhincome>0 & hhincome<9999999, valueh/hhincome, NA_real_)
  )

cat(sprintf("Working sample: %d rows\n", nrow(df)))

# ============================================================================
# PART 1: Tenure × move type × parenthood cross-tab
# ============================================================================

cat("\n========================================\n")
cat("PART 1: Tenure × Move × Parenthood\n")
cat("========================================\n")

movers <- df %>%
  filter(moved1y == 1, met2013 > 0, migmet131 > 0, !is.na(owner)) %>%
  mutate(
    same_msa = as.integer(met2013 == migmet131),
    diff_msa = as.integer(met2013 != migmet131)
  )

# Cross-tab: ownership rate by move type × parent status
crosstab <- movers %>%
  filter(!is.na(in_principal)) %>%
  mutate(
    move_detail = case_when(
      diff_msa == 1 ~ "Across MSA",
      same_msa == 1 & in_principal == 1 ~ "Within MSA → principal city",
      same_msa == 1 & in_principal == 0 ~ "Within MSA → suburb",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(move_detail)) %>%
  group_by(newparent, move_detail) %>%
  summarise(
    own_rate = weighted.mean(owner, perwt, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(parent_label = ifelse(newparent==1, "New Parents", "Non-Parents"))

cat("\nOwnership rate at destination, by move type and parent status:\n")
crosstab_wide <- crosstab %>%
  select(parent_label, move_detail, own_rate) %>%
  pivot_wider(names_from = parent_label, values_from = own_rate)
print(crosstab_wide)

# Plot
p1 <- ggplot(crosstab, aes(x = move_detail, y = own_rate, fill = parent_label)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = sprintf("%.0f%%", own_rate*100)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("New Parents" = "firebrick", "Non-Parents" = "navy")) +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult=c(0,0.15))) +
  labs(x = "", y = "Homeownership rate at destination",
       fill = "", title = "Ownership Rate by Move Type and Parenthood") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top", axis.text.x = element_text(size = 10))
ggsave(file.path(outdir, "ownership_by_movetype_parent.pdf"), p1, width = 9, height = 5.5)
ggsave(file.path(outdir, "ownership_by_movetype_parent.png"), p1, width = 9, height = 5.5, dpi = 150)
cat("Saved ownership_by_movetype_parent\n")


# ============================================================================
# PART 2: Regression — P(owner) among movers, by move type × parent
# ============================================================================

cat("\n========================================\n")
cat("PART 2: Regressions — P(owner) among movers\n")
cat("========================================\n")

movers_reg <- movers %>%
  filter(!is.na(in_principal)) %>%
  mutate(
    to_suburb = as.integer(same_msa == 1 & in_principal == 0),
    to_principal = as.integer(same_msa == 1 & in_principal == 1)
  )

cat("\n--- P(owner) ~ newparent, all movers ---\n")
r1 <- feols(owner ~ newparent + i(sex) | age + year, data = movers_reg, weights = ~perwt)
print(summary(r1))

cat("\n--- P(owner) ~ newparent, within-MSA suburb movers ---\n")
r2 <- feols(owner ~ newparent + i(sex) | age + year,
            data = movers_reg %>% filter(to_suburb == 1), weights = ~perwt)
print(summary(r2))

cat("\n--- P(owner) ~ newparent × to_suburb ---\n")
r3 <- feols(owner ~ newparent * to_suburb + i(sex) | age + year,
            data = movers_reg %>% filter(same_msa == 1), weights = ~perwt)
print(summary(r3))


# ============================================================================
# PART 3: Binscatter — MSA-level fertility vs ownership
# ============================================================================

cat("\n========================================\n")
cat("PART 3: Binscatter fertility vs ownership\n")
cat("========================================\n")

# Rebuild MSA-level data
msa <- df %>%
  filter(met2013 > 0) %>%
  group_by(met2013) %>%
  summarise(
    fert_rate = weighted.mean(has_youngchild, perwt, na.rm = TRUE),
    own_rate  = weighted.mean(owner, perwt, na.rm = TRUE),
    med_pti   = median(pti, na.rm = TRUE),
    pop       = sum(perwt),
    .groups   = "drop"
  ) %>%
  filter(pop >= 100000) # Larger MSAs only for cleaner binscatter

# Residualize (not strictly needed for bivariate, but clean)
# Create 20 equal-sized bins
msa <- msa %>%
  mutate(own_bin = ntile(own_rate, 20))

binned <- msa %>%
  group_by(own_bin) %>%
  summarise(
    own_rate  = weighted.mean(own_rate, pop),
    fert_rate = weighted.mean(fert_rate, pop),
    pop       = sum(pop),
    .groups   = "drop"
  )

p3 <- ggplot(binned, aes(x = own_rate, y = fert_rate)) +
  geom_point(size = 3.5, color = "navy") +
  geom_smooth(method = "lm", color = "firebrick", se = FALSE, linewidth = 1) +
  labs(x = "Homeownership rate (vingtile mean)",
       y = "Fertility rate (share with child < 5)",
       title = "Fertility vs Homeownership across MSAs",
       subtitle = "Binscatter (20 bins), ages 22-45, MSAs with pop > 100k") +
  theme_minimal(base_size = 14)
ggsave(file.path(outdir, "binscatter_fertility_ownership.pdf"), p3, width = 7, height = 5)
ggsave(file.path(outdir, "binscatter_fertility_ownership.png"), p3, width = 7, height = 5, dpi = 150)
cat("Saved binscatter_fertility_ownership\n")

# Same for PTI
msa_pti <- msa %>% filter(!is.na(med_pti)) %>% mutate(pti_bin = ntile(med_pti, 20))
binned_pti <- msa_pti %>%
  group_by(pti_bin) %>%
  summarise(
    med_pti   = weighted.mean(med_pti, pop),
    fert_rate = weighted.mean(fert_rate, pop),
    .groups   = "drop"
  )

p3b <- ggplot(binned_pti, aes(x = med_pti, y = fert_rate)) +
  geom_point(size = 3.5, color = "navy") +
  geom_smooth(method = "lm", color = "firebrick", se = FALSE, linewidth = 1) +
  labs(x = "Median price-to-income (vingtile mean)",
       y = "Fertility rate (share with child < 5)",
       title = "Fertility vs Housing Costs across MSAs",
       subtitle = "Binscatter (20 bins), ages 22-45, MSAs with pop > 100k") +
  theme_minimal(base_size = 14)
ggsave(file.path(outdir, "binscatter_fertility_pti.pdf"), p3b, width = 7, height = 5)
ggsave(file.path(outdir, "binscatter_fertility_pti.png"), p3b, width = 7, height = 5, dpi = 150)
cat("Saved binscatter_fertility_pti\n")


# ============================================================================
# PART 4: Fertility-ownership gradient separately for owners vs renters
# ============================================================================

cat("\n========================================\n")
cat("PART 4: MSA fertility gradient by tenure\n")
cat("========================================\n")

msa_tenure <- df %>%
  filter(met2013 > 0, !is.na(owner)) %>%
  group_by(met2013, owner) %>%
  summarise(
    fert_rate = weighted.mean(has_youngchild, perwt, na.rm = TRUE),
    pop       = sum(perwt),
    .groups   = "drop"
  ) %>%
  filter(pop >= 30000)

# Get MSA-level ownership rate (computed from everyone)
msa_own <- df %>%
  filter(met2013 > 0, !is.na(owner)) %>%
  group_by(met2013) %>%
  summarise(msa_own_rate = weighted.mean(owner, perwt), .groups = "drop")

msa_tenure <- msa_tenure %>%
  left_join(msa_own, by = "met2013") %>%
  mutate(tenure_label = ifelse(owner == 1, "Owners", "Renters"))

cat("\nRegression: fert_rate ~ msa_own_rate, OWNERS only:\n")
reg_own <- lm(fert_rate ~ msa_own_rate,
              data = msa_tenure %>% filter(owner == 1), weights = pop)
cat(sprintf("  Coefficient: %.4f (se: %.4f, p: %.4f)\n",
            coef(reg_own)[2],
            summary(reg_own)$coefficients[2,2],
            summary(reg_own)$coefficients[2,4]))

cat("\nRegression: fert_rate ~ msa_own_rate, RENTERS only:\n")
reg_ren <- lm(fert_rate ~ msa_own_rate,
              data = msa_tenure %>% filter(owner == 0), weights = pop)
cat(sprintf("  Coefficient: %.4f (se: %.4f, p: %.4f)\n",
            coef(reg_ren)[2],
            summary(reg_ren)$coefficients[2,2],
            summary(reg_ren)$coefficients[2,4]))

p4 <- ggplot(msa_tenure, aes(x = msa_own_rate, y = fert_rate, color = tenure_label)) +
  geom_point(aes(size = pop), alpha = 0.3) +
  geom_smooth(method = "lm", aes(weight = pop), se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = c("Owners" = "navy", "Renters" = "firebrick")) +
  labs(x = "MSA-level homeownership rate",
       y = "Fertility rate (share with child < 5)",
       color = "",
       title = "Fertility-Ownership Gradient: Owners vs Renters",
       subtitle = "Each dot = MSA × tenure group, ages 22-45") +
  theme_minimal(base_size = 14) +
  guides(size = "none")
ggsave(file.path(outdir, "fertility_ownership_by_tenure.pdf"), p4, width = 8, height = 6)
ggsave(file.path(outdir, "fertility_ownership_by_tenure.png"), p4, width = 8, height = 6, dpi = 150)
cat("Saved fertility_ownership_by_tenure\n")


# ============================================================================
# PART 5: Individual-level — P(young child) ~ owner, with MSA FE
# ============================================================================

cat("\n========================================\n")
cat("PART 5: Individual-level fertility ~ tenure\n")
cat("========================================\n")

indiv <- df %>% filter(met2013 > 0, !is.na(owner))

cat("\n--- P(child<5) ~ owner, with MSA + year FE ---\n")
r_indiv <- feols(has_youngchild ~ owner + i(sex) + i(race) | met2013 + year + age,
                 data = indiv, weights = ~perwt)
print(summary(r_indiv))

cat("\n--- Same regression, by age group ---\n")
indiv <- indiv %>% mutate(young = as.integer(age <= 32))

r_young <- feols(has_youngchild ~ owner + i(sex) + i(race) | met2013 + year + age,
                 data = indiv %>% filter(young == 1), weights = ~perwt)
r_old <- feols(has_youngchild ~ owner + i(sex) + i(race) | met2013 + year + age,
               data = indiv %>% filter(young == 0), weights = ~perwt)
cat(sprintf("  Young (22-32): owner coeff = %.4f (se = %.4f)\n",
            coef(r_young)["owner"], se(r_young)["owner"]))
cat(sprintf("  Older (33-45): owner coeff = %.4f (se = %.4f)\n",
            coef(r_old)["owner"], se(r_old)["owner"]))


# ============================================================================
# SUMMARY
# ============================================================================

cat("\n\n")
cat("################################################################\n")
cat("#                        SUMMARY                               #\n")
cat("################################################################\n\n")

cat("PART 1 — Tenure at destination by move type:\n")
print(crosstab_wide)

cat("\nPART 3 — Binscatter gradients saved (cleaner for slides)\n")

cat("\nPART 4 — Fertility-ownership MSA gradient by tenure:\n")
cat(sprintf("  Owners:  coeff = %.4f (p = %.4f)\n",
            coef(reg_own)[2], summary(reg_own)$coefficients[2,4]))
cat(sprintf("  Renters: coeff = %.4f (p = %.4f)\n",
            coef(reg_ren)[2], summary(reg_ren)$coefficients[2,4]))

cat(sprintf("\nPART 5 — Individual P(child<5) ~ owner | MSA+year+age FE:\n"))
cat(sprintf("  Full sample: %.4f (se %.4f)\n", coef(r_indiv)["owner"], se(r_indiv)["owner"]))
cat(sprintf("  Young 22-32: %.4f (se %.4f)\n", coef(r_young)["owner"], se(r_young)["owner"]))
cat(sprintf("  Older 33-45: %.4f (se %.4f)\n", coef(r_old)["owner"], se(r_old)["owner"]))

cat("\n=== ALL DONE ===\n")
cat(sprintf("Outputs in: %s\n", outdir))
