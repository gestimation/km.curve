test_that("ci.curve yields the same outputs as cif of mets", {
  testdata <- createTestData(20, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  library(mets)
  library(Rcpp)
  Surv <- km.curve:::Surv
  e <- cif(Event(t,epsilon)~1, data=testdata, cause=1)
  e_surv <- 1-e$mu
  e_time <- e$times
  expected <- as.numeric(c(e_surv[c(1,2,4,6,8,10,12,14)], e_time[c(1,2,4,6,8,10,12,14)]))
  t <- ci.curve(Surv(t, epsilon)~1, data=testdata, code.event=c(1,2), report.ggsurvfit = FALSE, conf.type = "none")
  tested <- as.numeric(c(t$surv[c(2:9)], t$time[c(2:9)]))
  expect_equal(tested, expected)
})

test_that("cuminc of cmprsk yields the same outputs as cif of mets", {
  testdata <- createTestData(20, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  library(mets)
  library(cmprsk)
  e <- cif(Event(t,epsilon)~1, data=testdata, cause=1)
  e_surv <- e$mu
  e_time <- e$times
  expected <- as.numeric(c(e_surv[c(1,2,4,6,8,10,12,14)], e_time[c(1,2,4,6,8,10,12,14)]))
#  print(expected)
  t <- cuminc(testdata$t, testdata$epsilon, cencode=0)
  tested <- t[1]
  expect_equal(tested, expected)
})



test_that("ci.curve yields the same outputs as cif of mets with strata", {
  testdata <- createTestData(20, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  library(mets)
  library(Rcpp)
  Surv <- km.curve:::Surv
  e <- cif(Event(t,epsilon)~strata(strata), data=testdata, cause=1)
  e_surv <- 1-e$mu
  e_time <- e$times
  e_strata <- e$strata
#  print(e_surv)
#  print(e_time)
#  print(e_strata)
#  expected <- as.numeric(c(e_surv[c(2,5,6,9,10,13,14)], e_time[c(2,5,6,9,10,13,14)]))
  expected <- as.numeric(c(e_surv[c(2,6,10,14)], e_time[c(2,6,10,14)]))
  t <- ci.curve(Surv(t, epsilon)~strata, data=testdata, code.event=c(1,2), report.ggsurvfit = FALSE, conf.type = "none")
  tested <- as.numeric(c(t$surv[c(2:5)], t$time[c(2:5)]))
  expect_equal(tested, expected)
})

test_that("ci.curve by strata yields the same outputs as subsetting", {
  testdata <- createTestData(20, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  testdata_f <- testdata[testdata$strata==FALSE,]
  library(Rcpp)
  Surv <- km.curve:::Surv
  e <- cif(Event(t,epsilon)~1, data=testdata, cause=1)
  e <- ci.curve(Surv(t, epsilon)~1, data=testdata_f, code.event=c(1,2), report.ggsurvfit = FALSE, conf.type = "none")
  expected <- c(e$surv, e$time)
  t <- ci.curve(Surv(t, epsilon)~strata, data=testdata, code.event=c(1,2), report.ggsurvfit = FALSE, conf.type = "none")
  tested <- c(t$surv[c(1:5)], t$time[c(1:5)])
  expect_equal(tested, expected)
})


test_that("coxphyields the same outputs as cifreg of mets for Fine-Gray model", {
  library(survival)
  library(mets)
  Surv <- survival:::Surv
  etime <- with(mgus2, ifelse(pstat==0, futime, ptime))
  event <- with(mgus2, ifelse(pstat==0, 2*death, 1))
  e <- cifreg(Event(etime,event)~age + sex,data=mgus2,cause=1,propodds=NULL)
  expected <- round(c(e$coef, e$se.coef),digit=3)
  etime <- with(mgus2, ifelse(pstat==0, futime, ptime))
  event <- with(mgus2, ifelse(pstat==0, 2*death, 1))
  event <- factor(event, 0:2, labels=c("censor", "pcm", "death"))
  pdata <- finegray(Surv(etime, event) ~ ., data=mgus2)
  t <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex, weight=fgwt, data=pdata)
  #  t <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex, weight=fgwt, cluster=id, data=pdata)
  tested <- round(c(t$coefficients, sqrt(t$var[1,1]), sqrt(t$var[2,2])),digit=3)
  expect_equal(tested, expected)
})

test_that("syntax in Examples of finegray of survival", {
  library(survival)
  Surv <- survival:::Surv
  etime <- with(mgus2, ifelse(pstat==0, futime, ptime))
  event <- with(mgus2, ifelse(pstat==0, 2*death, 1))
  event <- factor(event, 0:2, labels=c("censor", "pcm", "death"))
  pdata <- finegray(Surv(etime, event) ~ ., data=mgus2)
  tested <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex, weight=fgwt, data=pdata)
  expected <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex, weight=fgwt, cluster=id, data=pdata)
  expect_equal(round(tested$var,digit=3), round(expected$var,digit=3))
})



