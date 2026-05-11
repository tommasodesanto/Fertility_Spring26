################################################################################
# COMMUTE TIME: Do within-MSA movers / new parents have longer commutes?
# Uses acs_with_hybrid_codes_andHPI.dta (has trantime, pwmetro)
################################################################################

library(haven)
library(dplyr)
library(fixest)

outdir <- "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/KeyFacts"

cat("Loading data...\n")
df <- read_dta(
  "/Users/tommasodesanto/Desktop/Projects/Datasets/acs_with_hybrid_codes_andHPI.dta",
  col_select = c("year","age","sex","race","perwt","gq",
                 "metro","met2013","metarea","migmet1",
                 "migrate1","nchild","nchlt5","eldch",
                 "trantime","pwmetro","empstatd")
)
cat(sprintf("Loaded %d rows\n", nrow(df)))

df <- df %>%
  filter(age >= 22, age <= 45, gq %in% c(1,2)) %>%
  mutate(
    newparent = as.integer(eldch < 4 & eldch != 99 & nchild > 0),
    moved1y = as.integer(migrate1 >= 2 & !is.na(migrate1)),
    stayer = as.integer(migrate1 == 1),
    in_principal = case_when(
      metro == 2 ~ 1L, metro %in% c(1,3) ~ 0L, TRUE ~ NA_integer_
    ),
    employed = as.integer(empstatd %in% c(10,12)),
    # Clean commute time (0 = works at home or NA)
    commute = ifelse(trantime > 0 & trantime < 200 & employed == 1, trantime, NA_real_)
  )

cat(sprintf("Sample: %d, with valid commute: %d\n", nrow(df), sum(!is.na(df$commute))))

################################################################################
# 1. Simple means: commute by mover status × parent status
################################################################################

cat("\n=== Commute time means ===\n")

means <- df %>%
  filter(!is.na(commute)) %>%
  group_by(moved1y, newparent) %>%
  summarise(
    mean_commute = weighted.mean(commute, perwt, na.rm = TRUE),
    med_commute = median(commute, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    mover_label = ifelse(moved1y == 1, "Mover", "Stayer"),
    parent_label = ifelse(newparent == 1, "New Parent", "Non-Parent")
  )

cat("\n")
print(means %>% select(mover_label, parent_label, mean_commute, med_commute, n))

################################################################################
# 2. Regressions: commute ~ moved + newparent + interaction
################################################################################

cat("\n=== Regressions ===\n")

cat("\n--- (a) Commute ~ moved | age + year ---\n")
r1 <- feols(commute ~ moved1y + i(sex) + i(race) | age + year,
            data = df, weights = ~perwt)
cat(sprintf("  moved1y: %.2f min (se %.2f, p = %.4f)\n",
            coef(r1)["moved1y"], se(r1)["moved1y"], pvalue(r1)["moved1y"]))

cat("\n--- (b) Commute ~ newparent | age + year ---\n")
r2 <- feols(commute ~ newparent + i(sex) + i(race) | age + year,
            data = df, weights = ~perwt)
cat(sprintf("  newparent: %.2f min (se %.2f, p = %.4f)\n",
            coef(r2)["newparent"], se(r2)["newparent"], pvalue(r2)["newparent"]))

cat("\n--- (c) Commute ~ moved × newparent | age + year ---\n")
r3 <- feols(commute ~ moved1y * newparent + i(sex) + i(race) | age + year,
            data = df, weights = ~perwt)
print(summary(r3, keep = c("moved1y", "newparent", "moved1y:newparent")))

cat("\n--- (d) Same with MSA FE ---\n")
r4 <- feols(commute ~ moved1y * newparent + i(sex) + i(race) | age + year + met2013,
            data = df %>% filter(met2013 > 0), weights = ~perwt)
print(summary(r4, keep = c("moved1y", "newparent", "moved1y:newparent")))

################################################################################
# 3. Among movers: commute by destination type
################################################################################

cat("\n=== Among movers: commute by destination ===\n")

# Use metarea/migmet1 for within/across (these match in this dataset)
movers <- df %>%
  filter(moved1y == 1, !is.na(commute), !is.na(in_principal),
         metarea > 0, migmet1 > 0) %>%
  mutate(
    same_msa = as.integer(metarea == migmet1),
    diff_msa = as.integer(metarea != migmet1)
  )

# Hmm, check if metarea/migmet1 matching works in this dataset
cat(sprintf("Movers with MSA info: %d\n", nrow(movers)))
cat(sprintf("Same MSA: %d (%.1f%%)\n", sum(movers$same_msa),
            mean(movers$same_msa) * 100))

# If matching doesn't work (like before), fall back to principal city only
if (mean(movers$same_msa) < 0.05) {
  cat("WARNING: metarea/migmet1 matching broken, using met2013 proxy\n")
  # Can't do within/across with this dataset's metro codes
  # Just do commute by principal city status × parent
  cat("\n--- Commute by destination (principal city) × parent ---\n")
  means_dest <- movers %>%
    group_by(in_principal, newparent) %>%
    summarise(mean_commute = weighted.mean(commute, perwt),
              n = n(), .groups = "drop") %>%
    mutate(dest = ifelse(in_principal == 1, "Principal city", "Suburb"),
           parent = ifelse(newparent == 1, "New Parent", "Non-Parent"))
  print(means_dest %>% select(dest, parent, mean_commute, n))

  cat("\n--- Regression: commute ~ newparent × in_principal | age + year ---\n")
  r5 <- feols(commute ~ newparent * in_principal + i(sex) | age + year,
              data = movers, weights = ~perwt)
  print(summary(r5, keep = c("newparent", "in_principal", "newparent:in_principal")))

} else {
  # Full 4-way
  movers <- movers %>%
    mutate(
      cat4 = case_when(
        same_msa == 1 & in_principal == 1 ~ "Within MSA, principal",
        same_msa == 1 & in_principal == 0 ~ "Within MSA, suburb",
        diff_msa == 1 & in_principal == 1 ~ "Across MSA, principal",
        diff_msa == 1 & in_principal == 0 ~ "Across MSA, suburb"
      )
    )

  cat("\n--- Commute by 4-way destination × parent ---\n")
  means_4way <- movers %>%
    group_by(cat4, newparent) %>%
    summarise(mean_commute = weighted.mean(commute, perwt),
              n = n(), .groups = "drop") %>%
    mutate(parent = ifelse(newparent == 1, "New Parent", "Non-Parent"))
  print(means_4way %>% select(cat4, parent, mean_commute, n))

  cat("\n--- Regression: commute ~ newparent, within-MSA suburb movers ---\n")
  r5 <- feols(commute ~ newparent + i(sex) | age + year,
              data = movers %>% filter(same_msa == 1 & in_principal == 0),
              weights = ~perwt)
  cat(sprintf("  newparent: %.2f min (se %.2f, p = %.4f)\n",
              coef(r5)["newparent"], se(r5)["newparent"], pvalue(r5)["newparent"]))

  cat("\n--- Regression: commute ~ newparent × same_msa × in_principal ---\n")
  r6 <- feols(commute ~ newparent * same_msa * in_principal + i(sex) | age + year,
              data = movers, weights = ~perwt)
  print(summary(r6))
}

cat("\n=== DONE ===\n")
