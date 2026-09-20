source("second_birth_matching.R")

expect <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
dt <- data.table::data.table(
  anchor_person_key = rep("a", 5), anchor_household_key = rep("ha", 5), target_event_time = -5:-1,
  target_year = 2010L + (-5:-1), target_mother_age = 30L + (-5:-1),
  target_child_age = 0:4,
  target_gender.num = 2L, target_edlevel.num = 4L, target_marst.num = 1L,
  target_race.num = 1L, target_statefip.num = 50L,
  birth_gap = 5L, gap_full_pre = TRUE, gap_reference = TRUE
)
donors <- data.table::data.table(
  person_key = paste0("d", 1:7), household_key = paste0("h", 1:7),
  PERWT = c(10, 20, NA, 30, 40, 25, 12), YEAR = c(2005:2009, 2007L, 2005L),
  AGE_norm = c(25:29, 27L, 25L), sole_child_age = c(0:4, 2L, 0L),
  gender.num = 2L, edlevel.num = 4L, marst.num = 1L,
  race.num = 1L, statefip.num = 50L
)
donors[person_key == "d3", `:=`(YEAR = 2020L, AGE_norm = 40L)]
out <- second_birth_match_exact(dt, donors,
  match_covariates = c("gender.num", "edlevel.num", "marst.num", "race.num", "statefip.num"))
expect(nrow(out$links) == 6L, "exact-cell matching did not retain all tied controls")
expect(nrow(out$donor_weights) == 6L, "invalid-weight donor was not excluded")
expect(out$links[donor_person_key == "d1" & target_event_time == -5,
                 wgt_match] == 0.5, "tied controls were not normalized")
expect(out$donor_weights[person_key == "d1" & target_event_time == -5,
                         wgt] == 5, "vendor donor weight formula changed")
expect(out$anchor_support$full_pre_supported && out$anchor_support$reference_supported,
       "complete full-pre/reference support was not recognized")

# Two target rows in one cell create all ties and multiply donor matching weight.
dt2 <- rbind(dt, dt[1, ])
dt2$anchor_person_key[6] <- "b"
out2 <- second_birth_match_exact(dt2, donors,
  match_covariates = c("gender.num", "edlevel.num", "marst.num", "race.num", "statefip.num"))
expect(out2$donor_weights[person_key == "d1" & target_event_time == -5,
                          wgt_match] == 1, "donor-by-event reuse was not aggregated")
expect(out2$donor_weights[person_key == "d1" & target_event_time == -5,
                          wgt] == 10, "PERWT times normalized matching weight failed")

cells <- second_birth_author_cells(data.table::data.table(
  gender = c("Men", "Women"), edlevel = c("Below HS", "College +"),
  marst = c("Married, spouse present", "Divorced"),
  race = c("White, non-hispanic", "Hispanic"), statefip = c(50, 9)))
expect(identical(cells$gender.num, c(1L, 2L)) &&
       identical(cells$edlevel.num, c(1L, 2L)) &&
       identical(cells$marst.num, c(4L, 1L)) &&
       identical(cells$race.num, c(4L, 2L)) &&
       identical(cells$statefip.num, c("50", "09")),
       "author transformed-cell order changed")
unknown_cells <- second_birth_author_cells(data.table::data.table(
  gender = "Men", edlevel = "Below HS", marst = "Married, spouse present",
  race = "Other", statefip = 50))
expect(unknown_cells$author_cell_missing, "unrecognized race label was silently coarsened")
raw_cells <- second_birth_author_cells_raw(data.table::data.table(
  sex = c(1, 2), educd = c(62, 101), marst = c(1, 6), race = c(1, 2),
  hispan = c(0, 0), statefip = c(50, 9)))
expect(identical(raw_cells$gender.num, c(1L, 2L)) &&
       identical(raw_cells$edlevel.num, c(3L, 2L)) &&
       identical(raw_cells$marst.num, c(4L, 5L)) &&
       identical(raw_cells$race.num, c(4L, 1L)),
       "raw author recode does not match vendor cell definitions")

bad <- dt; bad$target_event_time[1] <- 0L
tryCatch({ second_birth_match_exact(bad, donors, c("gender.num", "edlevel.num", "marst.num", "race.num", "statefip.num"));
           stop("nonnegative target accepted", call. = FALSE) }, error = function(e) invisible(e))
tryCatch({ second_birth_match_exact(dt, donors, c("SEX", "EDUC", "MARST", "RACE", "STATEFIP"));
           stop("raw covariates accepted without author recodes", call. = FALSE) }, error = function(e) invisible(e))
cat("PASS: second-birth exact-cell ties, support gates, exclusions, and vendor weights\n")
