test_that("Surv yields the same outputs as Surv of survfit", {
  library(survival)
  df_test <- createTestData(100, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  expected <- survival::Surv(df_test$t, df_test$d)
  tested <- Surv(df_test$t, df_test$d)
  expect_equal(expected, tested)
})

test_that("Surv yields the same outputs as Surv of survfit with NA", {
  library(survival)
  df_test <- createTestData(100, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  df_test$t[1] <- NA
  df_test$d[1] <- NA
  expected <- survival::Surv(df_test$t, df_test$d)
  tested <- Surv(df_test$t, df_test$d)
  expect_equal(expected, tested)
})

test_that("Surv yields the same outputs as Surv of survfit with a factor", {
  library(survival)
  df_test <- createTestData(100, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  f <- as.factor(df_test$d)
  expected <- survival::Surv(df_test$t, f)
  tested <- Surv(df_test$t, f)
  expect_equal(expected[,2], tested[,2])
})

test_that("km.curve yields the same outputs as survfit", {
  library(survival)
  library(ggsurvfit)
  library(Rcpp)
  df_test <- createTestData(20, 2, first_zero=TRUE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  e <- survfit(Surv(t, d)~1, df_test, weight=w, na.action=na.omit, conf.type = "none")
  t <- km.curve(Surv(t, d)~1, df_test, weight="w", na.action=na.omit, conf.type = "none", use.ggsurvfit = FALSE)
  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("km.curve yields the same outputs as survfit when strata is logical", {
  library(survival)
  library(ggsurvfit)
  library(Rcpp)
  df_test <- createTestData(20, 1, first_zero=FALSE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  e <- survfit(Surv(t, d)~strata, df_test, weight=w, conf.type = "none")
  t <- km.curve(Surv(t, d)~strata, df_test, weight="w", conf.type = "none", use.ggsurvfit = FALSE)
  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("Jack Knife standard error of km.curve yields the same outputs as separate jack knifing when strata is present", {
  library(ggsurvfit)
  library(Rcpp)
  df_test <- createTestData(20, 1, first_zero=TRUE, last_zero=FALSE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  df_test1 <- subset(df_test, strata==FALSE)
  df_test2 <- subset(df_test, strata==TRUE)
  t <- km.curve(Surv(t, d)~strata, df_test, weight="w", error = "jackknife", use.ggsurvfit = FALSE)
  e1 <- km.curve(Surv(t, d)~1, df_test1, weight="w", error = "jackknife", use.ggsurvfit = FALSE)
  e2 <- km.curve(Surv(t, d)~1, df_test2, weight="w", error = "jackknife", use.ggsurvfit = FALSE)
  expected <- c(e1$std.err, e2$std.err)
  tested <- t$std.err
  expect_equal(expected, tested)
})

