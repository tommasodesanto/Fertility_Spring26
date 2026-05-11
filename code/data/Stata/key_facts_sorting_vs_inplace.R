################################################################################
# SORTING vs IN-PLACE: Is the MSA fertility gradient driven by who moves
# where, or by what people do where they live?
################################################################################

library(haven)
library(dplyr)
library(ggplot2)
library(tidyr)

outdir <- "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/KeyFacts"

cat("Loading data...\n")
df <- read_dta(
  "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta",
  col_select = c("year","age","sex","race","perwt","gq",
                 "metro","met2013","migmet131",
                 "migrate1","ownershp","nchild","nchlt5","eldch",
                 "valueh","hhincome")
)
cat(sprintf("Loaded %d rows\n", nrow(df)))

df <- df %>%
  filter(age >= 22, age <= 45, gq %in% c(1,2)) %>%
  mutate(
    owner = as.integer(ownershp == 1),
    has_youngchild = as.integer(nchlt5 > 0 & !is.na(nchlt5)),
    moved1y = as.integer(migrate1 >= 2 & !is.na(migrate1)),
    stayer = as.integer(migrate1 == 1),
    pti = ifelse(owner==1 & valueh>0 & valueh<9999999 &
                   hhincome>0 & hhincome<9999999, valueh/hhincome, NA_real_)
  )

################################################################################
# 1. MSA-level fertility gradient: stayers vs movers
################################################################################

cat("\n=== Fertility gradient: stayers vs movers ===\n")

# Collapse to MSA × mover status
msa_by_move <- df %>%
  filter(met2013 > 0, !is.na(moved1y)) %>%
  group_by(met2013, moved1y) %>%
  summarise(
    fert_rate = weighted.mean(has_youngchild, perwt, na.rm = TRUE),
    own_rate  = weighted.mean(owner, perwt, na.rm = TRUE),
    med_pti   = median(pti, na.rm = TRUE),
    pop       = sum(perwt),
    .groups   = "drop"
  ) %>%
  filter(pop >= 30000)

# Get MSA-level ownership (from everyone) to use as x-axis
msa_own <- df %>%
  filter(met2013 > 0) %>%
  group_by(met2013) %>%
  summarise(msa_own_rate = weighted.mean(owner, perwt),
            msa_pti = median(pti, na.rm = TRUE),
            .groups = "drop")

msa_by_move <- msa_by_move %>%
  left_join(msa_own, by = "met2013") %>%
  mutate(move_label = ifelse(moved1y == 0, "Stayers", "Movers (arrived last year)"))

# Regressions
cat("\nFertility ~ MSA ownership, STAYERS only:\n")
reg_stay <- lm(fert_rate ~ msa_own_rate,
               data = msa_by_move %>% filter(moved1y == 0), weights = pop)
cat(sprintf("  coeff = %.4f, se = %.4f, p = %.6f\n",
            coef(reg_stay)[2], summary(reg_stay)$coefficients[2,2],
            summary(reg_stay)$coefficients[2,4]))

cat("\nFertility ~ MSA ownership, MOVERS only:\n")
reg_move <- lm(fert_rate ~ msa_own_rate,
               data = msa_by_move %>% filter(moved1y == 1), weights = pop)
cat(sprintf("  coeff = %.4f, se = %.4f, p = %.6f\n",
            coef(reg_move)[2], summary(reg_move)$coefficients[2,2],
            summary(reg_move)$coefficients[2,4]))

# Plot
p1 <- ggplot(msa_by_move, aes(x = msa_own_rate, y = fert_rate,
                                color = move_label)) +
  geom_point(aes(size = pop), alpha = 0.25) +
  geom_smooth(method = "lm", aes(weight = pop), se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = c("Stayers" = "navy",
                                 "Movers (arrived last year)" = "firebrick")) +
  labs(x = "MSA-level homeownership rate",
       y = "Fertility rate (share with child < 5)",
       color = "",
       title = "MSA Fertility Gradient: Stayers vs Recent Movers",
       subtitle = "Ages 22-45, ACS") +
  theme_minimal(base_size = 14) +
  guides(size = "none")
ggsave(file.path(outdir, "fertility_gradient_stayers_vs_movers.pdf"), p1, width = 8, height = 6)
ggsave(file.path(outdir, "fertility_gradient_stayers_vs_movers.png"), p1, width = 8, height = 6, dpi = 150)
cat("Saved fertility_gradient_stayers_vs_movers\n")

# Same with PTI on x-axis
cat("\nFertility ~ MSA PTI, STAYERS only:\n")
reg_stay_pti <- lm(fert_rate ~ msa_pti,
                   data = msa_by_move %>% filter(moved1y == 0, !is.na(msa_pti)), weights = pop)
cat(sprintf("  coeff = %.4f, se = %.4f, p = %.6f\n",
            coef(reg_stay_pti)[2], summary(reg_stay_pti)$coefficients[2,2],
            summary(reg_stay_pti)$coefficients[2,4]))

cat("\nFertility ~ MSA PTI, MOVERS only:\n")
reg_move_pti <- lm(fert_rate ~ msa_pti,
                   data = msa_by_move %>% filter(moved1y == 1, !is.na(msa_pti)), weights = pop)
cat(sprintf("  coeff = %.4f, se = %.4f, p = %.6f\n",
            coef(reg_move_pti)[2], summary(reg_move_pti)$coefficients[2,2],
            summary(reg_move_pti)$coefficients[2,4]))

p1b <- ggplot(msa_by_move %>% filter(!is.na(msa_pti)),
              aes(x = msa_pti, y = fert_rate, color = move_label)) +
  geom_point(aes(size = pop), alpha = 0.25) +
  geom_smooth(method = "lm", aes(weight = pop), se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = c("Stayers" = "navy",
                                 "Movers (arrived last year)" = "firebrick")) +
  labs(x = "MSA-level median price-to-income",
       y = "Fertility rate (share with child < 5)",
       color = "",
       title = "MSA Fertility Gradient (PTI): Stayers vs Recent Movers",
       subtitle = "Ages 22-45, ACS") +
  theme_minimal(base_size = 14) +
  guides(size = "none")
ggsave(file.path(outdir, "fertility_gradient_pti_stayers_vs_movers.pdf"), p1b, width = 8, height = 6)
ggsave(file.path(outdir, "fertility_gradient_pti_stayers_vs_movers.png"), p1b, width = 8, height = 6, dpi = 150)
cat("Saved fertility_gradient_pti_stayers_vs_movers\n")


################################################################################
# 2. Binscatter versions (cleaner)
################################################################################

cat("\n=== Binscatter versions ===\n")

# Stayers binscatter
stayers_msa <- msa_by_move %>% filter(moved1y == 0)
stayers_msa <- stayers_msa %>% mutate(own_bin = ntile(msa_own_rate, 15))
bin_stay <- stayers_msa %>%
  group_by(own_bin) %>%
  summarise(msa_own_rate = weighted.mean(msa_own_rate, pop),
            fert_rate = weighted.mean(fert_rate, pop), .groups = "drop") %>%
  mutate(group = "Stayers")

movers_msa <- msa_by_move %>% filter(moved1y == 1)
movers_msa <- movers_msa %>% mutate(own_bin = ntile(msa_own_rate, 15))
bin_move <- movers_msa %>%
  group_by(own_bin) %>%
  summarise(msa_own_rate = weighted.mean(msa_own_rate, pop),
            fert_rate = weighted.mean(fert_rate, pop), .groups = "drop") %>%
  mutate(group = "Movers")

bin_both <- bind_rows(bin_stay, bin_move)

p2 <- ggplot(bin_both, aes(x = msa_own_rate, y = fert_rate, color = group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_color_manual(values = c("Stayers" = "navy", "Movers" = "firebrick")) +
  labs(x = "MSA homeownership rate (bin mean)",
       y = "Fertility rate (share with child < 5)",
       color = "",
       title = "Fertility-Ownership Gradient: Stayers vs Movers",
       subtitle = "Binscatter (15 bins), ages 22-45") +
  theme_minimal(base_size = 14)
ggsave(file.path(outdir, "binscatter_stayers_vs_movers.pdf"), p2, width = 8, height = 5.5)
ggsave(file.path(outdir, "binscatter_stayers_vs_movers.png"), p2, width = 8, height = 5.5, dpi = 150)
cat("Saved binscatter_stayers_vs_movers\n")


################################################################################
# 3. Decomposition: how much of MSA fertility variance is within vs between?
################################################################################

cat("\n=== Variance decomposition ===\n")

# Individual-level: how much does MSA explain vs mover status vs their interaction?
indiv <- df %>%
  filter(met2013 > 0, !is.na(moved1y)) %>%
  left_join(msa_own, by = "met2013")

# MSA-level fertility: overall vs stayers-only
msa_all <- df %>%
  filter(met2013 > 0) %>%
  group_by(met2013) %>%
  summarise(
    fert_all = weighted.mean(has_youngchild, perwt),
    fert_stayers = weighted.mean(has_youngchild[migrate1 == 1],
                                  perwt[migrate1 == 1], na.rm = TRUE),
    fert_movers = weighted.mean(has_youngchild[migrate1 >= 2 & !is.na(migrate1)],
                                 perwt[migrate1 >= 2 & !is.na(migrate1)], na.rm = TRUE),
    share_movers = weighted.mean(migrate1 >= 2 & !is.na(migrate1), perwt),
    pop = sum(perwt),
    .groups = "drop"
  ) %>%
  filter(pop >= 100000) %>%
  left_join(msa_own, by = "met2013")

cat("\nCorrelation between stayer fertility and overall fertility across MSAs:\n")
cat(sprintf("  r = %.4f\n", cor(msa_all$fert_stayers, msa_all$fert_all,
                                  use = "complete.obs")))

cat("\nCorrelation between mover fertility and overall fertility:\n")
cat(sprintf("  r = %.4f\n", cor(msa_all$fert_movers, msa_all$fert_all,
                                  use = "complete.obs")))

cat("\nDo MSAs with more movers have different fertility?\n")
reg_share <- lm(fert_all ~ share_movers + msa_own_rate, data = msa_all, weights = pop)
cat(sprintf("  Share movers coeff: %.4f (p = %.4f)\n",
            coef(reg_share)["share_movers"],
            summary(reg_share)$coefficients["share_movers", 4]))

# Scatter: stayer fertility vs mover fertility
p3 <- ggplot(msa_all, aes(x = fert_stayers, y = fert_movers, size = pop)) +
  geom_point(alpha = 0.4, color = "navy") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", aes(weight = pop), color = "firebrick",
              se = FALSE, linewidth = 1) +
  labs(x = "Fertility rate among stayers",
       y = "Fertility rate among recent movers",
       title = "MSA Fertility: Stayers vs Movers",
       subtitle = "Each dot = MSA. Dashed = 45-degree line.") +
  theme_minimal(base_size = 14) +
  guides(size = "none")
ggsave(file.path(outdir, "scatter_stayer_vs_mover_fertility.pdf"), p3, width = 7, height = 6)
ggsave(file.path(outdir, "scatter_stayer_vs_mover_fertility.png"), p3, width = 7, height = 6, dpi = 150)
cat("Saved scatter_stayer_vs_mover_fertility\n")

# Key question: are movers above or below the 45-degree line?
msa_all <- msa_all %>% mutate(mover_premium = fert_movers - fert_stayers)

cat("\nMover fertility minus stayer fertility across MSAs:\n")
cat(sprintf("  Mean: %.4f\n", weighted.mean(msa_all$mover_premium, msa_all$pop, na.rm = TRUE)))
cat(sprintf("  Median: %.4f\n", median(msa_all$mover_premium, na.rm = TRUE)))

cat("\nDoes the mover premium vary with MSA affordability?\n")
reg_prem <- lm(mover_premium ~ msa_own_rate, data = msa_all, weights = pop)
cat(sprintf("  coeff = %.4f, p = %.4f\n",
            coef(reg_prem)[2], summary(reg_prem)$coefficients[2,4]))

reg_prem2 <- lm(mover_premium ~ msa_pti, data = msa_all %>% filter(!is.na(msa_pti)), weights = pop)
cat(sprintf("  (PTI) coeff = %.4f, p = %.4f\n",
            coef(reg_prem2)[2], summary(reg_prem2)$coefficients[2,4]))


################################################################################
# 4. Summary
################################################################################

cat("\n\n")
cat("################################################################\n")
cat("#              SORTING vs IN-PLACE: SUMMARY                    #\n")
cat("################################################################\n\n")

cat("GRADIENT FOR STAYERS (ownership):\n")
cat(sprintf("  coeff = %.4f (p = %.6f)\n",
            coef(reg_stay)[2], summary(reg_stay)$coefficients[2,4]))

cat("GRADIENT FOR MOVERS (ownership):\n")
cat(sprintf("  coeff = %.4f (p = %.6f)\n",
            coef(reg_move)[2], summary(reg_move)$coefficients[2,4]))

cat("GRADIENT FOR STAYERS (PTI):\n")
cat(sprintf("  coeff = %.4f (p = %.6f)\n",
            coef(reg_stay_pti)[2], summary(reg_stay_pti)$coefficients[2,4]))

cat("GRADIENT FOR MOVERS (PTI):\n")
cat(sprintf("  coeff = %.4f (p = %.6f)\n",
            coef(reg_move_pti)[2], summary(reg_move_pti)$coefficients[2,4]))

cat(sprintf("\nMover-stayer fertility correlation: r = %.4f\n",
            cor(msa_all$fert_stayers, msa_all$fert_all, use = "complete.obs")))

cat(sprintf("Mean mover premium: %.4f\n",
            weighted.mean(msa_all$mover_premium, msa_all$pop, na.rm = TRUE)))

cat("\n=== DONE ===\n")
