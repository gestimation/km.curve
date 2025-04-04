#' #' Title Kaplan-Meier curve
#'
#' @param formula formula Model formula representing outcome and strata
#' @param data data.frame Input dataset containing survival data.
#' @param weights character Column name representing the weights. The weights must be nonnegative and it is strongly recommended that they be strictly positive, since zero weights are ambiguous, compared to use of the subset argument.
#' @param subset character Specifies a condition for subsetting the data. Defaults to NULL.
#' @param code.event integer Specifies the code of event. Defaults to 1.
#' @param code.censoring integer Specifies the code of censoring. Defaults to 0.
#' @param na.action character Specifies a missing-data filter function, applied to the model frame, after any subset argument has been used. Defaults to na.pass.
#' @param conf.int numeric The level for a two-sided confidence interval on the survival probabilities. Defaults to 0.95.
#' @param error character either the string "greenwood" for the Greenwood formula or "tsiatis" for the Tsiatis formula. Defaults to "greenwood".
#' @param conf.type character Specifies transformation used to construct the confidence interval on the probabilities. Defaults to "arcsine-square root".
#' @param use.ggsurvfit logical Draw a survival plot using ggsurvfit. Defaults to TRUE.
#' @param label.strata character Labels of strata. Defaults to NULL.
#' @param label.x character Labels of x axis. Defaults to "Survival probability".
#' @param label.y character Labels of y axis. Defaults to "Time".
#' @param lims.x vector Range of x axis. Defaults to NULL.
#' @param lims.y vector Range of y axis. Defaults to c(0, 1).
#' @param font.family character Specifies font family of the plot. Defaults to "sans".
#' @param font.size numeric Specifies font size of the plot. Defaults to 14.
#' @param legend.position character Specifies position of the legend of curves. Defaults to "top".
#' @importFrom ggsurvfit ggsurvfit
#' @importFrom Rcpp Rcpp
#'
#' @returns An object consists of Kaplan-Meier estimator and related statistics. This object is formatted to conform to the suvrfit class.
#' @export km.curve
#'
#' @examples
#' library(ggsurvfit)
#' library(Rcpp)
#' data(diabetes.complications)
#' diabetes.complications$d <- as.numeric(diabetes.complications$epsilon>0)
#' km.curve(Surv(t, d)~fruitq1, data=diabetes.complications)
km.curve <- function(formula,
                     data,
                     weights = NULL,
                     subset = NULL,
                     code.event = 1,
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
) {
  checkDependentPackages()
  out_readSurv <- readSurv(formula, data, weights, code.event, code.censoring, subset, na.action)
  sourceCpp('src/km.curve.cpp')
  out_km <- calculateKM_rcpp(out_readSurv$t, out_readSurv$d, out_readSurv$w, as.integer(out_readSurv$strata), error)
  if (!all(as.integer(out_readSurv$strata) == 1) & (is.null(label.strata))) {
    names(out_km$strata) <- levels(as.factor(out_readSurv$strata))
  } else if (!all(as.integer(out_readSurv$strata) == 1)) {
    names(out_km$strata) <- label.strata
  }
  out_ci <- calculateCI(out_km, conf.int, conf.type, conf.lower)
  if (is.null(lims.x)) {
    lims.x <- c(0, max(out_readSurv$t))
  }

  survfit_object <- list(
    time = out_km$time,
    surv = out_km$surv,
    n = out_km$n,
    n.risk = out_km$n.risk,
    n.event = out_km$n.event,
    n.censor = out_km$n.censor,
    std.err = out_km$std.err,
    upper = out_ci$upper,
    lower = out_ci$lower,
    conf.type = conf.type,
    call = match.call(),
    type = "kaplan-meier",
    method = "Kaplan-Meier"
  )
  if (any(as.integer(out_readSurv$strata) != 1)) {
    survfit_object$strata <- out_km$strata
  }

  class(survfit_object) <- c("survfit")
  if (use.ggsurvfit) {
    if (conf.type == "none" | conf.type == "n" | length(survfit_object$strata)>2) {
      survfit_object_ <- survfit_object
      survfit_object_$lower <- survfit_object_$surv
      survfit_object_$upper <- survfit_object_$surv
      out_ggsurvfit <- ggsurvfit(survfit_object_) +
        theme_classic()+
        theme(legend.position = legend.position,
              axis.title = element_text(size = (font.size+2), family = font.family),
              axis.text = element_text(size = font.size, family = font.family),
              legend.text = element_text(size = font.size, family = font.family))+
        labs(x = label.x,
             y = label.y) +
        lims(x = lims.x,
             y = lims.y) +
        theme(legend.position = "top")+
        add_risktable(risktable_stats = c("n.risk")) +
        add_censor_mark()
    } else {
      out_ggsurvfit <- ggsurvfit(survfit_object) +
        theme_classic()+
        theme(legend.position = legend.position,
              axis.title = element_text(size = (font.size+2), family = font.family),
              axis.text = element_text(size = font.size, family = font.family),
              legend.text = element_text(size = font.size, family = font.family))+
        labs(x = label.x,
             y = label.y) +
        lims(x = lims.x,
             y = lims.y) +
        add_confidence_interval() +
        add_risktable(risktable_stats = c("n.risk")) +
        add_censor_mark()
    }
    print(out_ggsurvfit)
  }
  return(survfit_object)
}
