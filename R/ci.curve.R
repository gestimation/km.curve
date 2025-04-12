#' #' Title Cumularive incidence curve
#'
#' @param formula formula Model formula representing outcome and strata
#' @param data data.frame Input dataset containing survival data.
#' @param weights character Column name representing the weights. The weights must be nonnegative and it is strongly recommended that they be strictly positive, since zero weights are ambiguous, compared to use of the subset argument.
#' @param subset character Specifies a condition for subsetting the data. Defaults to NULL.
#' @param code.event integer Specifies the code of event. Defaults to 1.
#' @param code.censoring integer Specifies the code of censoring. Defaults to 0.
#' @param na.action character Specifies a missing-data filter function, applied to the model frame, after any subset argument has been used. Defaults to na.pass.
#' @param conf.int numeric The level for a two-sided confidence interval on the survival probabilities. Defaults to 0.95.
#' @param error character Specifies standard error calculation. "greenwood" for the Greenwood formula, "tsiatis" for the Tsiatis formula or "jackknife" for the jack knife method. Defaults to "greenwood".
#' @param conf.type character Specifies transformation used to construct the confidence interval on the probabilities. Defaults to "arcsine-square root".
#' @param use.ggsurvfit logical Draw a survival plot using ggsurvfit. Defaults to TRUE.
#' @param label.strata character Labels of strata. Defaults to NULL.
#' @param label.x character Labels of x axis. Defaults to "Cumulative incidence probability".
#' @param label.y character Labels of y axis. Defaults to "Time".
#' @param lims.x vector Range of x axis. Defaults to NULL.
#' @param lims.y vector Range of y axis. Defaults to c(0, 1).
#' @param font.family character Specifies font family of the plot. Defaults to "sans".
#' @param font.size numeric Specifies font size of the plot. Defaults to 14.
#' @param legend.position character Specifies position of the legend of curves. Defaults to "top".
#' @importFrom ggsurvfit ggsurvfit
#' @importFrom Rcpp sourceCpp
#' @useDynLib km.curve, .registration = TRUE
#' @returns An object consists of Aalen-Johansen estimator and related statistics. This object is formatted to conform to the suvrfit class, so note that surv contains estimates corresponds to 1-cumulative incidence.
#' @export km.curve
#' @export ci.curve
#' @export Surv
#' @export Event
#' @export createTestData
#' @export calculateKM_rcpp
#'
#' @examples
#' library(ggsurvfit)
#' library(Rcpp)
#' testdata <- createTestData(200, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
#' ci.curve(Surv(t, epsilon)~1, data=testdata)
ci.curve <- function(formula,
                     data,
                     weights = NULL,
                     subset = NULL,
                     code.event = c(1, 2),
                     code.censoring = 0,
                     na.action = na.pass,
                     conf.int = 0.95,
                     error = "none",
                     conf.type = "none",
                     use.ggsurvfit = TRUE,
                     label.x = "Time",
                     label.y = "Cumulative incidence probability",
                     label.strata = NULL,
                     lims.x = NULL,
                     lims.y = c(0, 1),
                     font.family = "sans",
                     font.size = 14,
                     legend.position = "top"
) {
  checkDependentPackages()
  out_readSurv <- readSurv(formula, data, weights, code.event, code.censoring, subset, na.action)
  out_aj <- calculateAJ(out_readSurv)
  if (!all(as.integer(out_readSurv$strata) == 1) & (is.null(label.strata))) {
    names(out_aj$strata1) <- levels(as.factor(out_readSurv$strata))
  } else if (!all(as.integer(out_readSurv$strata) == 1)) {
    names(out_aj$strata1) <- label.strata
  }
  if (any(as.integer(out_readSurv$strata) != 1)) {
    n <- table(as.integer(out_readSurv$strata))
    rep_list <- mapply(rep, n, out_aj$strata1, SIMPLIFY = FALSE)
    n.risk <- do.call(c, rep_list) - out_aj$n.censor - out_aj$n.event1 - out_aj$n.event2
  } else {
    n <- length(out_readSurv$strata)
    n.risk <- n - out_aj$n.censor - out_aj$n.event1 - out_aj$n.event2
  }
  if (is.null(lims.x)) {
    lims.x <- c(0, max(out_readSurv$t))
  }

  survfit_object <- list(
    time = out_aj$time1,
    surv = (1-out_aj$aj1),
    n = n,
    n.risk = n.risk,
    n.event = out_aj$n.event1,
    n.censor = out_aj$n.censor,
    std.err = NULL,
    upper = NULL,
    lower = NULL,
    conf.type = conf.type,
    call = match.call(),
    type = "Aalen-Johansen",
    method = "Aalen-Johansen"
  )
  if (any(as.integer(out_readSurv$strata) != 1)) {
    survfit_object$strata <- out_aj$strata1
  }
  class(survfit_object) <- c("survfit")
  if (use.ggsurvfit) {
    if (conf.type == "none" | conf.type == "n" | length(survfit_object$strata)>2) {
      survfit_object_ <- survfit_object
      survfit_object_$std.err <- rep(1, length(survfit_object_$surv))
      survfit_object_$lower <- survfit_object_$surv
      survfit_object_$upper <- survfit_object_$surv
      out_ggsurvfit <- ggsurvfit(survfit_object_, type = "risk") +
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
      survfit_object$std.err <- rep(1, length(survfit_object$surv))
      survfit_object$lower <- survfit_object$surv
      survfit_object$upper <- survfit_object$surv
      out_ggsurvfit <- ggsurvfit(survfit_object, type = "risk") +
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

calculateAJ <- function(data) {
  #  out_km <- calculateKM_rcpp(out_readSurv$t, out_readSurv$d, out_readSurv$w, as.integer(out_readSurv$strata), error)
  out_km0 <- calculateKM_rcpp(data$t, data$d0, data$w, as.integer(data$strata), "none")
  km0 <- get_surv(data$t, out_km0$surv, out_km0$time, as.integer(data$strata), out_km0$strata)
  ip.weight <- (data$d0==0) * ifelse(km0 > 0, 1 / km0, 0)
  d1_ipw <- as.matrix(data$w*data$d1*ip.weight)

  aj1 <- NULL
  time1 <- NULL
  n.event1 <- 0
  n.event2 <- 0
  n.censor <- 0
  strata1 <- NULL
  for (level in sort(unique(as.integer(data$strata)))) {
    sub_d1_ipw <- d1_ipw[as.integer(data$strata) == level, ]
    sub_t <- data$t[as.integer(data$strata) == level]
    sub_d0 <- data$d0[as.integer(data$strata) == level]
    sub_d1 <- data$d1[as.integer(data$strata) == level]
    sub_d2 <- data$d2[as.integer(data$strata) == level]
    not_atrisk <- outer(sub_t, sub_t, ">=")
    #   not_atrisk <- matrix(as.integer(not_atrisk), nrow = nrow(not_atrisk), ncol = ncol(not_atrisk))
    sub_aj1 <- not_atrisk %*% sub_d1_ipw / length(sub_t)
    sub_n.censor <- not_atrisk %*% as.matrix(sub_d0)
    sub_n.event1 <- not_atrisk %*% as.matrix(sub_d1)
    sub_n.event2 <- not_atrisk %*% as.matrix(sub_d2)
    #    unique_aj1 <- unique(cbind(sub_t, sub_d1, sub_aj1))
    #    selected_aj1 <- unique_aj1[unique_aj1[, 2] == 1, ]
    #    sorted_aj1 <- selected_aj1[order(selected_aj1[, 1]), ]
    #    aj1 <- c(aj1, sorted_aj1[,3])
    #    time1 <- c(time1, sorted_aj1[,1])
    #    strata1 <- c(strata1, length(sorted_aj1[,2]))
    is_unique <- !duplicated(sub_t)
    unique_t <- sub_t[is_unique]
    unique_aj1 <- sub_aj1[is_unique]
    unique_n.event1 <- sub_n.event1[is_unique]
    unique_n.event2 <- sub_n.event2[is_unique]
    unique_n.censor <- sub_n.censor[is_unique]
    sorted_t <- sort(unique_t)
    sorted_aj1 <- unique_aj1[order(unique_t)]
    sorted_n.event1 <- unique_n.event1[order(unique_t)]
    sorted_n.event2 <- unique_n.event2[order(unique_t)]
    sorted_n.censor <- unique_n.censor[order(unique_t)]
    time1 <- c(time1, sorted_t)
    aj1 <- c(aj1, sorted_aj1)
    n.event1 <- c(n.event1, sorted_n.event1[-length(sorted_n.event1)])
    n.event2 <- c(n.event2, sorted_n.event2[-length(sorted_n.event2)])
    n.censor <- c(n.censor, sorted_n.censor[-length(sorted_n.censor)])
    strata1 <- c(strata1, length(sorted_t))
  }
 return(list(time1=time1, aj1=aj1,n.event1=n.event1, n.event2=n.event2, n.censor=n.censor, strata1=strata1))
}



