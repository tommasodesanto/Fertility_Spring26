options(repos=c(CRAN="https://cloud.r-project.org"), timeout=180, Ncpus=1)
stopifnot(grepl("aarch64|arm64", R.version$arch))
user_lib <- path.expand(Sys.getenv("R_LIBS_USER"))
dir.create(user_lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(user_lib, .libPaths()))
packages <- c("data.table", "haven", "fixest", "ggplot2", "digest", "dplyr", "readr",
              "tidyr", "scales", "stringr", "broom", "ggrepel", "here", "ipumsr",
              "jsonlite", "knitr", "matrixStats", "modelsummary", "patchwork", "plm",
              "purrr", "readxl", "sandwich", "sf", "tidycensus", "tidyselect",
              "tigris", "viridis", "GGally", "Hmisc", "kableExtra", "tidyverse", "usethis")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) install.packages(missing, lib=user_lib, type="binary")
ok <- vapply(packages, requireNamespace, logical(1), quietly=TRUE)
receipt <- data.frame(package=packages, load_ok=ok,
 version=vapply(packages, function(p) if(requireNamespace(p,quietly=TRUE)) as.character(packageVersion(p)) else NA_character_, character(1)))
write.csv(receipt,"output/setup/laptop_20260919/r/packages.csv",row.names=FALSE)
writeLines(capture.output(sessionInfo()),"output/setup/laptop_20260919/r/session.txt")
stopifnot(all(ok))
cat("FERTILITY_R_PACKAGES_PASS",length(packages),"packages\n")
