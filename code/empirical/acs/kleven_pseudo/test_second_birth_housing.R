#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(fixest)})
source("second_birth_housing.R", local=TRUE)
expect <- function(x,msg) if (!isTRUE(x)) stop(msg,call.=FALSE)
expect_error <- function(expr, msg) { ok <- FALSE; tryCatch({ force(expr) }, error=function(e) ok <<- TRUE); expect(ok, msg) }
source_rows <- data.table(person_key=c("a0","p1","p2","p3","d2","d1","unused"), YEAR=c(2010L,2010L,2010L,2010L,2008L,2008L,2010L), SAMPLE=1L, SERIAL=c(1:6,99L), PERNUM=1L, AGE_norm=c(30,30,30,30,28,28,40), STATEFIP=31L, SEX=2L, PERWT=c(10,10,10,10,5,5,NA_real_), rooms9=c(4,4,5,5,3,3,NA_real_), bedrooms5=c(2,2,3,3,1,1,NA_real_), ownership_lw=c(1,1,1,1,0,0,NA_real_))
mk <- function(k,t) data.table(anchor_person_key="a0",anchor_household_key="h0",target_event_time=t,target_year=2010+t,target_mother_age=30+t,target_child_age=5+t,donor_person_key=k,donor_household_key=paste0("h",k),donor_PERWT=5,wgt_match=1)
links <- rbindlist(list(mk("d2",-2),mk("d1",-1))); target_support <- data.table(anchor_person_key="a0",target_event_time=c(-2L,-1L),has_donor=TRUE)
proxy <- list(input=source_rows,anchors=data.table(person_key="a0",birth_gap=5,FERTYR_status="yes"),post_rows=data.table(person_key=c("a0","p1","p2","p3"),event_time=0:3,birth_gap=5,PERWT=10,FERTYR_status="yes"),config=list())
matched <- list(links=links,target_support=target_support)
allp <- prepare_second_birth_housing(proxy,matched,source_rows)
expect(nrow(allp$post)==4L,"post rows were dropped or anchors appended"); expect(sum(allp$panel$pseudo_role=="post" & allp$panel$event_time==0L)==1L,"event-0 duplicated"); expect(all(allp$post[order(event_time),weight]==10),"post PERWT changed")
joint <- prepare_second_birth_housing(proxy,matched,source_rows,support_spec="joint_negative")
expect(identical(joint$support$n_joint_negative,1L),"joint negative support gate failed"); expect(nrow(joint$negative_donor)==2L,"joint negative donor cells missing"); expect(all(joint$negative_donor$weight==c(5,5)),"donor weights changed")
no_proxy <- proxy; no_proxy$anchors <- copy(proxy$anchors); no_proxy$anchors[, FERTYR_status := "no"]; no_proxy$post_rows <- copy(proxy$post_rows); no_proxy$post_rows[event_time == 0, FERTYR_status := "no"]
no <- prepare_second_birth_housing(no_proxy, matched, source_rows, fertyr_spec="event0_yes")
expect(nrow(no$post) == 3L && nrow(no$negative_donor) == 0L && all(no$post$anchor_population == "event0_yes_anchor_validation"), "FERTYR event-0 sensitivity did not restrict anchors and t0 only")
sup <- second_birth_housing_support(allp); expect(all(sup[n_observed > 0, weight_ess] >= 1),"weight ESS missing")
expected <- c(-2L,-1L,0L,1L,2L,3L)
fit_times <- rep(expected, each=2L)
fit_panel <- rbindlist(lapply(seq_along(fit_times), function(i) { t <- fit_times[i]; data.table(
  person_key = paste0("f", i), YEAR = 2010L, SAMPLE = 1L, SERIAL = i,
  AGE_norm = 30, STATEFIP = 31L, event_time = t, weight = 1, rooms9 = t + 3,
  pseudo_role = "post", source_household_cluster = paste0("f", i)) }))
fit <- fit_second_birth_housing(list(panel=fit_panel),outcomes="rooms9")
expect("rooms9" %in% names(fit$fits),"pooled fit missing")
expect(abs(fit$contrasts$rooms9$estimate - 4) < 1e-8,"known +3 minus -1 contrast failed")
expect(all(fit$curves[outcome == "rooms9", event_time] == expected),"exported six-cell curves missing")
expect(abs(fit$curves[outcome == "rooms9" & event_time == -2L, baseline_raw_mean] - 1) < 1e-8 &&
       fit$curves[outcome == "rooms9" & event_time == -2L, reference_event] &&
       fit$curves[outcome == "rooms9" & event_time == -2L, std_error] == 0 &&
       abs(fit$curves[outcome == "rooms9" & event_time == -1L, pre_minus1_coefficient] - 1) < 1e-8 &&
       all(c("weight_ess", "raw_event_mean", "contrast_plus3_minus_minus1") %in% names(fit$curves)),
       "baseline, ESS, pre-minus-one, or contrast diagnostics missing")
expect(all(dim(fit$fits$rooms9$contrast$vcov) == c(2L,2L)), "full contrast covariance missing")
expect(fit$fits$rooms9$nobs == stats::nobs(fit$fits$rooms9$fit) && fit$fits$rooms9$nobs_prefit >= fit$fits$rooms9$nobs, "fit observation counts missing")
expect(grepl("source-household",fit$specification,fixed=TRUE),"cluster contract missing")
bad_panel <- fit_panel[event_time != 3L]
expect_error(fit_second_birth_housing(list(panel=bad_panel),outcomes="rooms9"))
bad_link <- copy(matched); bad_link$links <- copy(matched$links); bad_link$links[1, donor_person_key := "missing_source"]
expect_error(prepare_second_birth_housing(proxy, bad_link, source_rows), "missing matched source key did not fail")
bad_proxy <- proxy; bad_proxy$post_rows <- copy(proxy$post_rows); bad_proxy$post_rows[person_key == "p1", PERWT := NA_real_]
expect_error(prepare_second_birth_housing(bad_proxy, matched, source_rows), "invalid used post weight did not fail")
cat("PASS: second-birth pooled housing population, support sensitivity, weights, and fit\n")
