#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("cannot resolve test path")
driver <- file.path(dirname(normalizePath(sub("^--file=", "", file_arg))),
                    "run_first_birth_short_window.R")
status <- system2(Sys.which("Rscript"), c(driver), env = "SHORT_WINDOW_TEST=1")
if (!identical(status, 0L)) stop("short-window fixture failed with status ", status)
cat("test_first_birth_short_window PASS\n")
