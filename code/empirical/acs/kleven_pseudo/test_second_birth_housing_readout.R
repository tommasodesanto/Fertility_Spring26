#!/usr/bin/env Rscript
# Tiny fixture for the readout's common-cohort support gate.
suppressPackageStartupMessages(library(data.table))
path <- Sys.getenv("READOUT_SCRIPT", file.path("code", "empirical", "acs", "kleven_pseudo", "run_second_birth_housing_readout.R"))
exprs <- parse(path)
helper <- exprs[vapply(exprs, function(e) is.call(e) && identical(as.character(e[[1L]]), "<-"), logical(1L))]
helper <- helper[vapply(helper, function(e) identical(as.character(e[[2L]]), "common_cohort_gate"), logical(1L))]
stopifnot(length(helper) == 1L)
eval(helper[[1L]], envir=.GlobalEnv)

panel <- CJ(implied_event_year=2007:2016, event_time=c(-2L,-1L,0L,1L,2L,3L))
panel[,YEAR:=implied_event_year+event_time][,implied_event_year:=NULL]
panel[,`:=`(rooms9=1, bedrooms5=2, ownership_lw=.75, weight=1)]
panel <- rbind(panel, data.table(YEAR=2006L, event_time=0L, rooms9=99, bedrooms5=99, ownership_lw=99, weight=1))
g <- common_cohort_gate(panel, c("rooms9","bedrooms5","ownership_lw"))
stopifnot(nrow(g$rows) == 180L, nrow(g$panel) == 60L,
          all(g$summary$requested_event_cells == 60L),
          all(g$summary$observed_event_cells == 60L), all(g$summary$supported_event_cells == 60L),
          all(g$summary$all_six_supported),
          all(g$panel$implied_event_year >= 2007L & g$panel$implied_event_year <= 2016L),
          !any(g$panel$rooms9 == 99))

panel_missing <- panel[event_time != 3L]
gm <- common_cohort_gate(panel_missing, c("rooms9","bedrooms5","ownership_lw"))
stopifnot(all(gm$summary$requested_event_cells == 60L),
          all(gm$summary$observed_event_cells == 50L),
          all(gm$summary$supported_event_cells == 50L),
          all(!gm$summary$all_six_supported),
          all(gm$rows[event_time == 3L, n_rows] == 0L),
          all(gm$rows[event_time == 3L, n_observed] == 0L))
cat("READOUT_GATE_FIXTURE_PASS\n")
