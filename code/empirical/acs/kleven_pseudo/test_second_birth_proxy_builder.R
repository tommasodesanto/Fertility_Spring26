# Tiny deterministic tests for build_second_birth_proxy.R.  No ACS file is
# loaded; all rows below are hand-built.

source("build_second_birth_proxy.R")
library(data.table)

expect_error <- function(expr, pattern) {
  got <- tryCatch({force(expr); NULL}, error = function(e) conditionMessage(e))
  stopifnot(!is.null(got), grepl(pattern, got, fixed = TRUE))
}

row <- function(year, serial, pernum, momloc, age, sex = 2L, nchild = 0L,
                fertyr = 0L, perwt = 100, educ = 3L, marst = 1L,
                race = 1L, statefip = 50L) {
  data.table(YEAR = year, SAMPLE = 202301L, SERIAL = serial,
             PERNUM = pernum, MOMLOC = momloc, AGE = age, SEX = sex,
             NCHILD = nchild, FERTYR = fertyr, PERWT = perwt,
             EDUC = educ, MARST = marst, RACE = race, STATEFIP = statefip)
}

# 2022 household: strict event-0 anchor, gap five, FERTYR unknown (0=N/A).
d <- rbindlist(list(
  row(2022, 100, 1, 0, 34, nchild = 2),
  row(2022, 100, 2, 1, 5, sex = 1L),
  row(2022, 100, 3, 1, 0, sex = 1L),
  # 2024 cross-section: post row with the same structural IDs but t=2.
  row(2024, 100, 1, 0, 36, nchild = 2, fertyr = 0L),
  row(2024, 100, 2, 1, 7, sex = 1L),
  row(2024, 100, 3, 1, 2, sex = 1L),
  # Gap-one event-0 anchor: valid at -1, impossible at reference -2.
  row(2021, 101, 1, 0, 30, nchild = 2, fertyr = 2L),
  row(2021, 101, 2, 1, 1, sex = 1L),
  row(2021, 101, 3, 1, 0, sex = 1L),
  # Third-child/twin ambiguity: wide sensitivity can see it, strict cannot.
  row(2020, 102, 1, 0, 31, nchild = 3, fertyr = 2L),
  row(2020, 102, 2, 1, 4, sex = 1L),
  row(2020, 102, 3, 1, 0, sex = 1L),
  row(2020, 102, 4, 1, 0, sex = 1L),
  # NCHILD/link mismatch: only one valid linked child.
  row(2020, 103, 1, 0, 32, nchild = 2, fertyr = 2L),
  row(2020, 103, 2, 1, 3, sex = 1L),
  # Child older than mother: link must be retained in links but invalidated.
  row(2020, 104, 1, 0, 20, nchild = 2, fertyr = 2L),
  row(2020, 104, 2, 1, 25, sex = 1L),
  row(2020, 104, 3, 1, 0, sex = 1L),
  # Observed FERTYR No at t=0: excluded.
  row(2020, 105, 1, 0, 30, nchild = 2, fertyr = 1L),
  row(2020, 105, 2, 1, 4, sex = 1L),
  row(2020, 105, 3, 1, 0, sex = 1L),
  # Age sentinel 999 is reported as missing, not treated as a real age.
  row(2020, 106, 1, 0, 999, nchild = 2, fertyr = 2L),
  row(2020, 106, 2, 1, 4, sex = 1L),
  row(2020, 106, 3, 1, 0, sex = 1L),
  # One-child donor with unchanged PERWT and author covariates.
  row(2020, 107, 1, 0, 30, nchild = 1L),
  row(2020, 107, 2, 1, 4, sex = 1L),
  # Valid links but mother age at proxy birth is below the strict band.
  row(2020, 108, 1, 0, 24, nchild = 2, fertyr = 2L),
  row(2020, 108, 2, 1, 4, sex = 1L),
  row(2020, 108, 3, 1, 0, sex = 1L),
  # Self MOMLOC link is retained and invalidated.
  row(2020, 109, 1, 1, 30, nchild = 2, fertyr = 2L),
  row(2020, 109, 2, 1, 4, sex = 1L),
  row(2020, 109, 3, 1, 0, sex = 1L),
  # t=11 is outside the configured nonnegative post window (+10).
  row(2020, 110, 1, 0, 40, nchild = 2, fertyr = 0L),
  row(2020, 110, 2, 1, 16, sex = 1L),
  row(2020, 110, 3, 1, 11, sex = 1L),
  # Event-0 FERTYR missing, unknown 8, and unrecognized 7.
  row(2020, 111, 1, 0, 32, nchild = 2, fertyr = NA_integer_),
  row(2020, 111, 2, 1, 3, sex = 1L),
  row(2020, 111, 3, 1, 0, sex = 1L),
  row(2020, 112, 1, 0, 32, nchild = 2, fertyr = 8L),
  row(2020, 112, 2, 1, 3, sex = 1L),
  row(2020, 112, 3, 1, 0, sex = 1L),
  row(2020, 113, 1, 0, 32, nchild = 2, fertyr = 7L),
  row(2020, 113, 2, 1, 3, sex = 1L),
  row(2020, 113, 3, 1, 0, sex = 1L),
  # Invalid raw person weight is flagged and excluded from strict output.
  row(2020, 114, 1, 0, 32, nchild = 2, fertyr = 2L, perwt = 0),
  row(2020, 114, 2, 1, 3, sex = 1L),
  row(2020, 114, 3, 1, 0, sex = 1L),
  # Missing child AGE is retained in links with an explicit reason.
  row(2020, 115, 1, 0, 32, nchild = 2, fertyr = 2L),
  row(2020, 115, 2, 1, NA_integer_, sex = 1L),
  row(2020, 115, 3, 1, 0, sex = 1L),
  # Female household row with no MOMLOC-linked child: placeholder is not a child link.
  row(2020, 116, 1, 0, 30, nchild = 0L),
  # A child in a different household has the same PERNUM/MOMLOC and must not link.
  row(2022, 999, 2, 1, 9, sex = 1L)
))

out <- build_second_birth_proxy(
  d,
  fertyr_codes = list(yes = 2L, no = 1L, unknown = c(0L, 8L)),
  match_covariates = c("SEX", "EDUC", "MARST", "RACE", "STATEFIP")
)

# Exact household-scoped linkage and unchanged person weights.
anchor_key <- "2022/202301/100/1"
anchor <- out$anchors[person_key == anchor_key]
stopifnot(nrow(anchor) == 1L, anchor$event_time == 0,
          anchor$event_year == 2022, anchor$age_at_event == 34,
          anchor$birth_gap == 5, anchor$FERTYR_status == "unknown",
          anchor$PERWT == 100)
stopifnot(out$mother_rows[person_key == anchor_key, linked_child_count] == 2L)
stopifnot(!any(!is.na(out$links$mother_person_key) &
               out$links$mother_person_key == anchor_key &
               out$links$SERIAL == 999L))
stopifnot(nrow(out$one_child_donors[person_key == "2020/202301/107/1"]) == 1L,
          out$one_child_donors[person_key == "2020/202301/107/1", sole_child_age] == 4,
          out$one_child_donors[person_key == "2020/202301/107/1", PERWT] == 100)

# Full pre-window support is available at gap five; targets are anchor-only.
stopifnot(anchor$gap_full_pre, anchor$gap_reference)
stopifnot(nrow(out$donor_targets[anchor_person_key == anchor_key]) == 5L)
stopifnot(all(out$donor_targets[anchor_person_key == anchor_key,
                                  target_child_age] >= 0))
stopifnot(!any(out$donor_targets$anchor_person_key == "2024/202301/100/1"))

# Gap one contributes at -1 but not at reference -2; no negative target exists.
gap1_key <- "2021/202301/101/1"
gap1 <- out$anchors[person_key == gap1_key]
stopifnot(nrow(gap1) == 1L, gap1$birth_gap == 1,
          !gap1$gap_full_pre, !gap1$gap_reference)
stopifnot(nrow(out$donor_targets[anchor_person_key == gap1_key]) == 1L,
          out$donor_targets[anchor_person_key == gap1_key, target_event_time] == -1)

# Third child/twin ambiguity, mismatch, impossible age link, FERTYR No, and age
# sentinel are all explicit exclusions rather than silent rewrites.
stopifnot(!any(out$anchors$person_key == "2020/202301/102/1"))
stopifnot(out$mother_rows[person_key == "2020/202301/102/1", exclusion_flag] ==
            "ambiguous_age_tie")
stopifnot(out$mother_rows[person_key == "2020/202301/103/1", exclusion_flag] ==
            "nchild_link_mismatch")
stopifnot(any(out$links$SERIAL == 104 & out$links$link_invalid_reason ==
              "child age not younger than mother"))
stopifnot(out$mother_rows[person_key == "2020/202301/105/1", exclusion_flag] ==
            "fertyr_observed_no_at_event0")
stopifnot(is.na(out$mother_rows[person_key == "2020/202301/106/1", AGE_norm]))
stopifnot(out$mother_rows[person_key == "2020/202301/108/1", exclusion_flag] ==
            "age_at_event_outside_band")
stopifnot(any(out$links$SERIAL == 109 & out$links$momloc_self),
          any(out$links$SERIAL == 109 & out$links$link_invalid_reason ==
                "MOMLOC self-link"))
stopifnot(out$mother_rows[person_key == "2020/202301/110/1", event_time] == 11,
          out$mother_rows[person_key == "2020/202301/110/1", strict_eligible],
          !any(out$post_rows$person_key == "2020/202301/110/1"),
          out$audit[category == "strict_rows_outside_event_window", count] >= 1L)
stopifnot(out$mother_rows[person_key == "2020/202301/111/1", strict_eligible],
          !out$mother_rows[person_key == "2020/202301/111/1", fertyr_event0_complete],
          out$mother_rows[person_key == "2020/202301/112/1", strict_eligible],
          out$mother_rows[person_key == "2020/202301/112/1", FERTYR_status] == "unknown",
          !out$mother_rows[person_key == "2020/202301/113/1", strict_eligible],
          out$mother_rows[person_key == "2020/202301/113/1", exclusion_flag] ==
            "fertyr_unrecognized_raw_code")
stopifnot(!out$mother_rows[person_key == "2020/202301/114/1", strict_eligible],
          out$mother_rows[person_key == "2020/202301/114/1", exclusion_flag] ==
            "missing_or_nonpositive_perwt")
stopifnot(any(out$links$SERIAL == 115 & out$links$link_invalid_reason ==
              "missing child age"),
          out$mother_rows[person_key == "2020/202301/115/1", exclusion_flag] ==
            "nchild_link_mismatch")
stopifnot(any(out$links$SERIAL == 116 & out$links$link_invalid_reason ==
              "mother_without_link"),
          out$audit[category == "mother_without_link_rows", count] >= 1L,
          out$audit[category == "invalid_missing_child_age_links", count] == 1L)
stopifnot(all(is.na(out$mother_rows[strict_eligible == TRUE, exclusion_flag])))

# Required metadata, explicit FERTYR codes, and duplicate person keys fail loudly.
missing_weight <- copy(d)
missing_weight[, PERWT := NULL]
expect_error(build_second_birth_proxy(missing_weight,
                                      fertyr_codes = list(yes = 2L, no = 1L, unknown = 0L),
                                      match_covariates = "SEX"),
             "missing required metadata fields")
expect_error(build_second_birth_proxy(d,
                                      fertyr_codes = list(yes = 2L, no = 1L),
                                      match_covariates = "SEX"),
             "fertyr_codes must explicitly contain")
dup <- rbind(d, d[1])
expect_error(build_second_birth_proxy(dup,
                                      fertyr_codes = list(yes = 2L, no = 1L, unknown = 0L),
                                      match_covariates = "SEX"),
             "person key")
missing_key <- copy(d)
missing_key[1, YEAR := NA_integer_]
expect_error(build_second_birth_proxy(missing_key,
                                      fertyr_codes = list(yes = 2L, no = 1L, unknown = 0L),
                                      match_covariates = "SEX"),
             "missing/nonfinite")
nonfinite_key <- copy(d)
nonfinite_key[1, SERIAL := Inf]
expect_error(build_second_birth_proxy(nonfinite_key,
                                      fertyr_codes = list(yes = 2L, no = 1L, unknown = 0L),
                                      match_covariates = "SEX"),
             "missing/nonfinite")

cat("second_birth_proxy_builder tests: OK (all deterministic checks passed)\n")
