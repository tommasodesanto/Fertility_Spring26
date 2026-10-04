#!/usr/bin/env Rscript
# Same assessment fields and sample as the historical cache, selecting 2023.
suppressPackageStartupMessages({library(haven);library(data.table);library(jsonlite)})
out <- commandArgs(trailingOnly=TRUE)[1]
source <- "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
fields <- c("year","RELTOHEAD_","AGEREP","IW","NETWORTHR","EARNINDRRC","HOMEOWN","HOMEVALUER","ACTUALROOMS_")
x <- as.data.table(read_dta(source,col_select=tidyselect::all_of(fields)))
survey_year <- max(x[RELTOHEAD_==10 & AGEREP>=18 & AGEREP<=85 & is.finite(IW) & IW>0 & is.finite(NETWORTHR),year],na.rm=TRUE)
cat("Latest PSID wave:",survey_year,"\n")
x <- x[year==survey_year & RELTOHEAD_==10 & AGEREP>=18 & AGEREP<=85 & is.finite(IW) & IW>0]
if (!nrow(x)) stop("PSID shelf lacks eligible recent observations")
y <- x[,.(year=as.integer(year),age=as.numeric(AGEREP),weight=as.numeric(IW),total_net_wealth=as.numeric(NETWORTHR),
 annual_gross_labor_earnings=as.numeric(EARNINDRRC),owner=fifelse(HOMEOWN==1,1,fifelse(HOMEOWN==2,0,NA_real_)),
 gross_home_value=fifelse(HOMEOWN==1,as.numeric(HOMEVALUER),fifelse(HOMEOWN==2,0,NA_real_)),
 rooms=fifelse(is.finite(ACTUALROOMS_) & ACTUALROOMS_>0 & ACTUALROOMS_!=99,as.numeric(ACTUALROOMS_),NA_real_))]
p <- file.path(out,"psid_recent.csv");fwrite(y,p,na="")
sha <- sub(" .*","",system2("shasum",c("-a","256",p),stdout=TRUE)[1])
write_json(list(selected_cache_sha256=sha,source_path=source,source_size=file.info(source)$size,
 source_mtime=as.character(file.info(source)$mtime),year=survey_year,rows=nrow(y),
 filters="RELTOHEAD_=10; age18..85; finite IW>0; same fields and outcome masks as historical assessment",
 currency="2022 USD as labeled in source; plots normalize within sample"),file.path(out,"metadata.json"),auto_unbox=TRUE,pretty=TRUE)
cat("PSID recent wave cached:",nrow(y),"rows\n")
