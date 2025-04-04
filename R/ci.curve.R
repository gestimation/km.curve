ci.curve <- function(formula,
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
  out_km0 <- calculateKM_rcpp(out_readSurv$t, out_readSurv$d0, out_readSurv$w, as.integer(out_readSurv$strata), error="none")

  km0 <- get_surv(out_readSurv$t, out_km0$surv, out_km0$time, as.integer(out_readSurv$strata), out_km0$strata)
  ip.weight <- (out_readSurv$d0==0) * ifelse(km0 > 0, 1 / km0, 0)
  atrisk <- outer(out_readSurv$t, out_readSurv$t, "<=")
  aj <- atrisk %*% as.matrix(out_readSurv$d*ip.weight)/length(out_readSurv$t)
#  if (!all(as.integer(out_readSurv$strata) == 1) & (is.null(label.strata))) {
#    names(out_km$strata) <- levels(as.factor(out_readSurv$strata))
#  } else if (!all(as.integer(out_readSurv$strata) == 1)) {
#    names(out_km$strata) <- label.strata
#  }
#  out_ci <- calculateCI(out_km, conf.int, conf.type, conf.lower)
  if (is.null(lims.x)) {
    lims.x <- c(0, max(out_readSurv$t))
  }
}
