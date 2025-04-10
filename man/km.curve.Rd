% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/km.curve.R
\name{km.curve}
\alias{km.curve}
\title{#' Title Kaplan-Meier curve}
\usage{
km.curve(
  formula,
  data,
  weights = NULL,
  subset = NULL,
  code.event = c(1, 2),
  code.censoring = 0,
  na.action = na.pass,
  conf.int = 0.95,
  error = "greenwood",
  conf.type = "arcsine-square root",
  use.ggsurvfit = TRUE,
  label.x = "Time",
  label.y = "Survival probability",
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

\item{error}{character Specifies standard error calculation. "greenwood" for the Greenwood formula, "tsiatis" for the Tsiatis formula or "jackknife" for the jack knife method. Defaults to "greenwood".}

\item{conf.type}{character Specifies transformation used to construct the confidence interval on the probabilities. Defaults to "arcsine-square root".}

\item{use.ggsurvfit}{logical Draw a survival plot using ggsurvfit. Defaults to TRUE.}

\item{label.x}{character Labels of x axis. Defaults to "Survival probability".}

\item{label.y}{character Labels of y axis. Defaults to "Time".}

\item{label.strata}{character Labels of strata. Defaults to NULL.}

\item{lims.x}{vector Range of x axis. Defaults to NULL.}

\item{lims.y}{vector Range of y axis. Defaults to c(0, 1).}

\item{font.family}{character Specifies font family of the plot. Defaults to "sans".}

\item{font.size}{numeric Specifies font size of the plot. Defaults to 14.}

\item{legend.position}{character Specifies position of the legend of curves. Defaults to "top".}
}
\value{
An object consists of Kaplan-Meier estimator and related statistics. This object is formatted to conform to the suvrfit class.
}
\description{
#' Title Kaplan-Meier curve
}
\examples{
library(ggsurvfit)
library(Rcpp)
testdata <- createTestData(200, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
km.curve(Surv(t, d)~strata, data=testdata)
}
