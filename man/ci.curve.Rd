% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ci.curve.R
\name{ci.curve}
\alias{ci.curve}
\title{#' Title Estimate and plot cumulative incidence curves}
\usage{
ci.curve(
  formula,
  data,
  weights = NULL,
  subset = NULL,
  code.event = c(1, 2),
  code.censoring = 0,
  na.action = na.pass,
  conf.int = 0.95,
  error = "delta",
  conf.type = "arcsine-square root",
  report.survfit.std.err = FALSE,
  report.ggsurvfit = TRUE,
  label.x = "Time",
  label.y = "Cumulative incidence probability",
  label.strata = NULL,
  lims.x = NULL,
  lims.y = c(0, 1),
  font.family = "sans",
  font.size = 14,
  legend.position = "top"
)
}
\arguments{
\item{formula}{formula Model formula representing outcome and strata}

\item{data}{data.frame Input dataset containing survival data.}

\item{weights}{character Column name representing the weights. The weights must be nonnegative and it is strongly recommended that they be strictly positive, since zero weights are ambiguous, compared to use of the subset argument.}

\item{subset}{character Specifies a condition for subsetting the data. Defaults to NULL.}

\item{code.event}{integer Specifies the code of event. Defaults to 1.}

\item{code.censoring}{integer Specifies the code of censoring. Defaults to 0.}

\item{na.action}{character Specifies a missing-data filter function, applied to the model frame, after any subset argument has been used. Defaults to na.pass.}

\item{conf.int}{numeric The level for a two-sided confidence interval on the survival probabilities. Defaults to 0.95.}

\item{error}{character Specifies standard error calculation. "aalen" for the Aalen formula, "delta" for the delta method or "jackknife" for the jack knife method. Defaults to "delta".}

\item{conf.type}{character Specifies transformation used to construct the confidence interval on the probabilities. Defaults to "arcsine-square root".}

\item{report.survfit.std.err}{logical Report standard error of log of survival probabilities. If this is not specified, the SE of survival probabilities is stored in std.err, unlike the original survfitThis is  according to the original survfit. Defaults to FALSE.}

\item{report.ggsurvfit}{logical Draw a survival plot using ggsurvfit. Defaults to TRUE.}

\item{label.x}{character Labels of x axis. Defaults to "Cumulative incidence probability".}

\item{label.y}{character Labels of y axis. Defaults to "Time".}

\item{label.strata}{character Labels of strata. Defaults to NULL.}

\item{lims.x}{vector Range of x axis. Defaults to NULL.}

\item{lims.y}{vector Range of y axis. Defaults to c(0, 1).}

\item{font.family}{character Specifies font family of the plot. Defaults to "sans".}

\item{font.size}{numeric Specifies font size of the plot. Defaults to 14.}

\item{legend.position}{character Specifies position of the legend of curves. Defaults to "top".}
}
\value{
An object consists of Aalen-Johansen estimator and related statistics. This object is formatted to conform to the survfit class, so note that surv contains estimates corresponds to 1-cumulative incidence. Some methods for the class (e.g. residuals.survfit) are not supported.
}
\description{
#' Title Estimate and plot cumulative incidence curves
}
\examples{
library(km.curve)
library(gtsummary)
library(dplyr)
library(labelled)
data(prostate)
prostate <- prostate \%>\% mutate(epsilon=ifelse(status=="alive",0,
                                       ifelse(status=="dead - prostatic ca",1,
                                       ifelse(status=="dead - other ca",1,
                                       ifelse(status=="dead - heart or vascular",2,
                                       ifelse(status=="dead - cerebrovascular",2,2))))))
prostate <- prostate \%>\% mutate(a=ifelse(rx=="placebo","Placebo","Experimental"))
prostate$t <- prostate$dtime/12
attr(prostate$a, "label") <- "Treatment"
survfit_by_group <- ci.curve(Event(t, epsilon)~a, data=prostate, code.event = c(1,2), label.y = "Cumulative incidence of cancer-specific death", label.x = "Years from randomization")
ci.curve(Event(t, epsilon)~a, data=prostate, code.event = c(2,1), label.y = "Cumulative incidence of non-cancer-specific death", label.x = "Years from randomization", conf.type = "n")

survfit_overall <- ci.curve(Event(t, epsilon)~1, data=prostate, code.event = c(1,2), report.ggsurvfit=FALSE)
survfit_list <- list(survfit_overall, survfit_by_group)
table_from_survfit <- tbl_survfit(survfit_list, times = c(2, 4, 6), label_header = "**{time} years**", type="risk") |>
  modify_spanning_header(all_stat_cols() ~ "**Cumulative incidence of cancer-specific death**")
table_from_survfit
}
