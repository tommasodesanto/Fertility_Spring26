#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(ggplot2)
  library(haven)
})

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  weighted.mean(x[ok], w[ok])
}

weighted_quantile_safe <- function(x, w, prob = 0.5) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) / sum(w) >= prob)[1]]
}

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[startsWith(args, "--file=")]
  if (length(file_arg) == 0L) stop("Run this file with Rscript.")
  dirname(normalizePath(sub("--file=", "", file_arg[1L], fixed = TRUE)))
}

make_sample <- function(dt, sample_name) {
  if (sample_name == "legacy_person_pooled_rep_6575") {
    out <- dt[age >= 65 & age <= 75 & is.finite(n_children_rep)]
    out[, n_children := n_children_rep]
    return(out)
  }
  if (sample_name == "reference_pooled_num_6575") {
    out <- dt[reference_person & age >= 65 & age <= 75 & is.finite(n_children_num)]
    out[, n_children := n_children_num]
    return(out)
  }
  if (sample_name == "reference_closest70_num_6575") {
    out <- dt[reference_person & age >= 65 & age <= 75 & is.finite(n_children_num)]
    out[, age_distance := abs(age - 70)]
    setorder(out, id, age_distance, -year)
    out <- out[, .SD[1L], by = id]
    out[, age_distance := NULL]
    out[, n_children := n_children_num]
    return(out)
  }
  if (sample_name == "reference_pooled_num_7684") {
    out <- dt[reference_person & age >= 76 & age <= 84 & is.finite(n_children_num)]
    out[, n_children := n_children_num]
    return(out)
  }
  if (sample_name == "reference_closest80_num_7684") {
    out <- dt[reference_person & age >= 76 & age <= 84 & is.finite(n_children_num)]
    out[, age_distance := abs(age - 80)]
    setorder(out, id, age_distance, -year)
    out <- out[, .SD[1L], by = id]
    out[, age_distance := NULL]
    out[, n_children := n_children_num]
    return(out)
  }
  stop("Unknown sample: ", sample_name)
}

summarize_block <- function(block, sample_name, asset_name, asset_col, group_name, group_value) {
  x <- block[[asset_col]]
  ratio <- x / block$income
  ok_asset <- is.finite(x)
  ok_ratio <- is.finite(ratio) & block$income > 1000
  data.table(
    sample = sample_name,
    asset = asset_name,
    fertility_group = group_name,
    fertility_value = group_value,
    n_rows_asset = sum(ok_asset),
    n_people_asset = uniqueN(block$id[ok_asset]),
    n_rows_ratio = sum(ok_ratio),
    n_people_ratio = uniqueN(block$id[ok_ratio]),
    weight_ratio = sum(block$weight[ok_ratio], na.rm = TRUE),
    wealth_mean = weighted_mean_safe(x, block$weight),
    wealth_median = weighted_quantile_safe(x, block$weight),
    income_mean = weighted_mean_safe(block$income, block$weight),
    ratio_mean_of_ratios = weighted_mean_safe(ratio[ok_ratio], block$weight[ok_ratio]),
    ratio_median = weighted_quantile_safe(ratio[ok_ratio], block$weight[ok_ratio]),
    ratio_of_means = weighted_mean_safe(x[ok_ratio], block$weight[ok_ratio]) /
      weighted_mean_safe(block$income[ok_ratio], block$weight[ok_ratio])
  )
}

root <- normalizePath(file.path(script_dir(), "..", "..", ".."))
psid_path <- file.path(dirname(root), "PSID", "PSIDSHELF_MOBILITY.dta")
outdir <- file.path(
  root,
  "code", "data", "psid_followup_mar2026", "output",
  "intergen_bequest_family_size_audit"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(psid_path)) stop("Missing PSID shelf: ", psid_path)

message("Reading selected PSID columns...")
raw <- as.data.table(read_dta(
  psid_path,
  col_select = c(
    "ID", "year", "AGEREP", "DEATHYEAR", "RELTOHEAD_",
    "RELCHIREP", "RELCHINUM", "INCFAMR", "NETWORTHR",
    "NETWORTH2R", "HOMEEQUITYR", "HOMEOWN", "IW",
    "FAMMARRIED", "EDUYEAR", "SEX", "RACE"
  )
))

dt <- raw[, .(
  id = as.numeric(ID),
  year = as.integer(year),
  age = as.numeric(AGEREP),
  death_year = as.numeric(DEATHYEAR),
  relation_to_head = as.numeric(RELTOHEAD_),
  n_children_rep = as.numeric(RELCHIREP),
  n_children_num = as.numeric(RELCHINUM),
  income = as.numeric(INCFAMR),
  total_nw = as.numeric(NETWORTHR),
  nonhousing_nw = as.numeric(NETWORTH2R),
  home_equity = as.numeric(HOMEEQUITYR),
  own = fifelse(as.numeric(HOMEOWN) == 1, 1, fifelse(as.numeric(HOMEOWN) == 2, 0, NA_real_)),
  weight = as.numeric(IW),
  married = as.numeric(FAMMARRIED),
  education_years = as.numeric(EDUYEAR),
  sex = as.numeric(SEX),
  race = as.numeric(RACE)
)]
rm(raw)

dt <- dt[
  year >= 1984 & year <= 2019 &
    (is.na(death_year) | year <= death_year) &
    age >= 65 & age <= 84 &
    is.finite(weight) & weight > 0
]
dt[, reference_person := relation_to_head == 10]
dt[, child_count_agrees := n_children_rep == n_children_num]

consistency <- dt[, .(
  n_rows = .N,
  n_people = uniqueN(id),
  n_reference_rows = sum(reference_person, na.rm = TRUE),
  n_both_child_counts = sum(is.finite(n_children_rep) & is.finite(n_children_num)),
  n_child_counts_disagree = sum(
    is.finite(n_children_rep) & is.finite(n_children_num) & !child_count_agrees
  ),
  n_total_identity_observed = sum(
    is.finite(total_nw) & is.finite(nonhousing_nw) & is.finite(home_equity)
  ),
  max_abs_total_identity_gap = max(
    abs(total_nw - nonhousing_nw - home_equity),
    na.rm = TRUE
  )
)]
fwrite(consistency, file.path(outdir, "definition_consistency.csv"))

sample_names <- c(
  "legacy_person_pooled_rep_6575",
  "reference_pooled_num_6575",
  "reference_closest70_num_6575",
  "reference_pooled_num_7684",
  "reference_closest80_num_7684"
)
assets <- c(
  nonhousing = "nonhousing_nw",
  home_equity = "home_equity",
  total = "total_nw"
)

summary_rows <- list()
sample_count_rows <- list()
for (sample_name in sample_names) {
  sample <- copy(make_sample(dt, sample_name))
  sample[, fertility_bin := fcase(
    n_children == 0, "0",
    n_children == 1, "1",
    n_children == 2, "2",
    n_children >= 3, "3+",
    default = NA_character_
  )]
  sample_count_rows[[sample_name]] <- sample[, .(
    n_rows = .N,
    n_people = uniqueN(id),
    mean_age = weighted_mean_safe(age, weight),
    childless_share = weighted_mean_safe(n_children == 0, weight),
    mean_children = weighted_mean_safe(n_children, weight),
    owner_rate = weighted_mean_safe(own, weight)
  )][, sample := sample_name]

  groups <- c("all", "childless", "parent", "0", "1", "2", "3+")
  for (asset_name in names(assets)) {
    asset_col <- unname(assets[[asset_name]])
    for (group_name in groups) {
      block <- if (group_name == "all") {
        sample
      } else if (group_name == "childless") {
        sample[n_children == 0]
      } else if (group_name == "parent") {
        sample[n_children > 0]
      } else {
        sample[fertility_bin == group_name]
      }
      summary_rows[[length(summary_rows) + 1L]] <- summarize_block(
        block, sample_name, asset_name, asset_col, group_name,
        if (group_name %in% c("0", "1", "2", "3+")) group_name else NA_character_
      )
    }
  }
}

summaries <- rbindlist(summary_rows, fill = TRUE)
sample_counts <- rbindlist(sample_count_rows, fill = TRUE, use.names = TRUE)
setcolorder(sample_counts, c("sample", setdiff(names(sample_counts), "sample")))
fwrite(summaries, file.path(outdir, "wealth_income_by_fertility.csv"))
fwrite(sample_counts, file.path(outdir, "sample_counts.csv"))

gap_rows <- list()
for (sample_name in sample_names) {
  for (asset_name in names(assets)) {
    tab <- summaries[sample == sample_name & asset == asset_name]
    val <- function(group, stat) tab[fertility_group == group][[stat]][1L]
    gap_rows[[length(gap_rows) + 1L]] <- data.table(
      sample = sample_name,
      asset = asset_name,
      parent_minus_childless_mean_ratio =
        val("parent", "ratio_mean_of_ratios") - val("childless", "ratio_mean_of_ratios"),
      parent_minus_childless_median_ratio =
        val("parent", "ratio_median") - val("childless", "ratio_median"),
      two_minus_one_mean_ratio =
        val("2", "ratio_mean_of_ratios") - val("1", "ratio_mean_of_ratios"),
      threeplus_minus_one_mean_ratio =
        val("3+", "ratio_mean_of_ratios") - val("1", "ratio_mean_of_ratios"),
      two_minus_one_median_ratio =
        val("2", "ratio_median") - val("1", "ratio_median"),
      threeplus_minus_one_median_ratio =
        val("3+", "ratio_median") - val("1", "ratio_median")
    )
  }
}
gaps <- rbindlist(gap_rows)
fwrite(gaps, file.path(outdir, "fertility_wealth_gaps.csv"))

# Diagnostic regressions only: these show whether the raw completed-fertility
# pattern is explained by observable composition. They are not candidate SMM
# targets because the model does not contain every control used here.
reg_sample <- copy(make_sample(dt, "reference_pooled_num_6575"))
reg_sample[, fertility_bin := factor(fcase(
  n_children == 0, "0",
  n_children == 1, "1",
  n_children == 2, "2",
  n_children >= 3, "3+",
  default = NA_character_
), levels = c("0", "1", "2", "3+"))]
reg_sample[, `:=`(age_sq = age^2, log_income = NA_real_)]
reg_sample[is.finite(income) & income > 0, log_income := log(income)]
regression_rows <- list()
for (asset_name in names(assets)) {
  asset_col <- unname(assets[[asset_name]])
  reg_sample[, outcome := asinh(get(asset_col) / income)]
  for (spec_name in c("age_year_income", "full_observables")) {
    rhs <- if (spec_name == "age_year_income") {
      "i(fertility_bin, ref = '0') + age + age_sq + log_income"
    } else {
      paste(
        "i(fertility_bin, ref = '0') + age + age_sq + log_income +",
        "married + education_years + i(sex) + i(race)"
      )
    }
    fml <- as.formula(paste0("outcome ~ ", rhs, " | year"))
    fit <- feols(
      fml,
      data = reg_sample[is.finite(outcome) & income > 1000],
      weights = ~weight,
      vcov = ~id,
      notes = FALSE
    )
    ct <- as.data.table(coeftable(fit), keep.rownames = "term")
    setnames(ct, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
             c("estimate", "std_error", "t_value", "p_value"))
    ct <- ct[grepl("^fertility_bin::", term)]
    ct[, `:=`(
      asset = asset_name,
      specification = spec_name,
      outcome_definition = "asinh(wealth / annual family income)",
      reference_group = "0 children",
      n_obs = nobs(fit)
    )]
    regression_rows[[length(regression_rows) + 1L]] <- ct
  }
}
regressions <- rbindlist(regression_rows, fill = TRUE)
setcolorder(
  regressions,
  c(
    "asset", "specification", "outcome_definition", "reference_group",
    "term", "estimate", "std_error", "t_value", "p_value", "n_obs"
  )
)
fwrite(regressions, file.path(outdir, "controlled_fertility_gradients.csv"))

published <- data.table(
  moment = c(
    "old_nonhousing_wealth_to_income_median_6575",
    "old_parent_childless_nonhousing_wealth_to_income_gap_6575",
    "old_total_wealth_to_income_median_6575",
    "old_parent_childless_total_wealth_to_income_gap_6575"
  ),
  stored_value = c(2.23046078, 1.00744952, 5.26375360, 1.28413467)
)
comparison <- rbindlist(lapply(sample_names[1:3], function(sample_name) {
  nh <- summaries[sample == sample_name & asset == "nonhousing"]
  tw <- summaries[sample == sample_name & asset == "total"]
  gg <- gaps[sample == sample_name]
  data.table(
    sample = sample_name,
    moment = published$moment,
    stored_value = published$stored_value,
    recomputed_value = c(
      nh[fertility_group == "all", ratio_median][1L],
      gg[asset == "nonhousing", parent_minus_childless_mean_ratio][1L],
      tw[fertility_group == "all", ratio_median][1L],
      gg[asset == "total", parent_minus_childless_mean_ratio][1L]
    )
  )
}))
comparison[, difference := recomputed_value - stored_value]
fwrite(comparison, file.path(outdir, "stored_target_comparison.csv"))

plot_data <- summaries[
  sample %in% c("reference_pooled_num_6575", "reference_closest70_num_6575") &
    fertility_group %in% c("0", "1", "2", "3+")
]
plot_data[, fertility_group := factor(fertility_group, levels = c("0", "1", "2", "3+"))]
plot_data[, sample_label := fifelse(
  sample == "reference_pooled_num_6575",
  "Pooled family-years",
  "One observation per reference person"
)]
p <- ggplot(plot_data, aes(fertility_group, ratio_median, color = asset, group = asset)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~sample_label) +
  labs(
    x = "Completed children",
    y = "Weighted median wealth / annual family income",
    color = "Wealth object"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(
  file.path(outdir, "wealth_income_by_completed_fertility.png"),
  p, width = 9, height = 4.8, dpi = 180
)

primary_gap <- gaps[sample == "reference_pooled_num_6575"]
primary_all <- summaries[
  sample == "reference_pooled_num_6575" & fertility_group == "all"
]
legacy_gap <- gaps[sample == "legacy_person_pooled_rep_6575"]
full_total <- regressions[asset == "total" & specification == "full_observables"]
readme <- c(
  "# Intergenerational bequest family-size target audit",
  "",
  "This packet recomputes the old-age wealth targets before any bequest recalibration.",
  "The primary sample uses one PSID reference person per family-year, completed-child",
  "records (`RELCHINUM`), ages 65--75, annual family income above $1,000 for ratios,",
  "and the longitudinal individual weight `IW`.",
  "",
  "## Mechanical findings",
  "",
  sprintf(
    "- `RELCHIREP` and `RELCHINUM` disagreements when both are observed: %d.",
    consistency$n_child_counts_disagree
  ),
  sprintf(
    "- Maximum absolute gap in `NETWORTHR = NETWORTH2R + HOMEEQUITYR`: %.6g.",
    consistency$max_abs_total_identity_gap
  ),
  sprintf(
    "- Legacy nonhousing parent-minus-childless mean-ratio gap: %.6f.",
    legacy_gap[asset == "nonhousing", parent_minus_childless_mean_ratio]
  ),
  sprintf(
    "- Reference-person nonhousing parent-minus-childless mean-ratio gap: %.6f.",
    primary_gap[asset == "nonhousing", parent_minus_childless_mean_ratio]
  ),
  sprintf(
    "- Reference-person total-wealth parent-minus-childless mean-ratio gap: %.6f.",
    primary_gap[asset == "total", parent_minus_childless_mean_ratio]
  ),
  sprintf(
    "- Reference-person nonhousing median wealth/income: %.6f.",
    primary_all[asset == "nonhousing", ratio_median]
  ),
  sprintf(
    "- Reference-person total median wealth/income: %.6f.",
    primary_all[asset == "total", ratio_median]
  ),
  sprintf(
    paste0(
      "- Fully controlled total-wealth coefficients relative to childless: ",
      "one child %.3f (p=%.3f), two children %.3f (p=%.3f), ",
      "three-plus %.3f (p=%.3f)."
    ),
    full_total[term == "fertility_bin::1", estimate],
    full_total[term == "fertility_bin::1", p_value],
    full_total[term == "fertility_bin::2", estimate],
    full_total[term == "fertility_bin::2", p_value],
    full_total[term == "fertility_bin::3+", estimate],
    full_total[term == "fertility_bin::3+", p_value]
  ),
  "",
  "## Interpretation boundary",
  "",
  "These are descriptive candidate moments, not approved SMM targets. The fertility-bin",
  "table is the decision object for whether a free child bequest shifter has empirical",
  "support. The controlled regressions are diagnostic only because the model lacks some",
  "of their covariates. Sampling uncertainty, marital-unit alignment, and a model-side",
  "Jacobian remain necessary before promotion.",
  "",
  "## Files",
  "",
  "- `stored_target_comparison.csv`: old constants versus legacy and corrected samples.",
  "- `wealth_income_by_fertility.csv`: levels by asset and completed-fertility bin.",
  "- `fertility_wealth_gaps.csv`: parent-childless and within-parent fertility gaps.",
  "- `controlled_fertility_gradients.csv`: composition-adjusted descriptive gradients.",
  "- `sample_counts.csv`: sample sizes and basic composition.",
  "- `definition_consistency.csv`: fertility-variable and wealth-identity checks.",
  "- `wealth_income_by_completed_fertility.png`: compact visual comparison."
)
writeLines(readme, file.path(outdir, "README.md"))

message("Wrote audit packet to ", outdir)
