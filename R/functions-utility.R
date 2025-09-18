checkDependentPackages <- function() {
  if (requireNamespace("ggsurvfit", quietly = TRUE) & requireNamespace("Rcpp", quietly = TRUE)) {
    suppressWarnings(library(ggsurvfit))
    suppressWarnings(library(Rcpp))
  } else {
    stop("Required packages 'ggsurvfit' and/or 'Rcpp' are not installed.")
  }
}

Surv <- function(time, event) {
  if (missing(time))
    stop("Must have a time argument")
  if (!is.numeric(time))
    stop("Time variable is not numeric")
  if (any(is.na(time)))
    warning("Invalid time variable. NA values included")
  #  if (any(time<0))
  #    warning("Invalid time variable. Non-negative values included")
  if (length(event) != length(time))
    stop("Time and event variables are different lengths")
  if (missing(event))
    stop("Must have an event argument")
  if (is.numeric(event)) {
    if (any(is.na(event)))
      warning("Invalid event variable. NA values included")
    status <- event
  } else if (is.logical(event)) {
    if (any(is.na(event)))
      warning("Invalid event variable. NA values included")
    status <- as.numeric(event)
    warning("Event variable is logical, converted to numearic")
  } else if (is.factor(event)) {
    status <- as.numeric(as.factor(event)) - 1
    warning("Event variable is a factor, converted to numearic")
  } else stop("Invalid status value, must be logical or numeric")
  if (nlevels(as.factor(event)) > 3)
    warning("Event variable should not have more than three levels")
  ss <- cbind(time = time, status = status)
  type <- "right"

  inputAttributes <- list()
  if (!is.null(attributes(time)))
    inputAttributes$time <- attributes(time)
  cname <- dimnames(ss)[[2]]
  if (length(cname) == 0) {
    cname <- c("time", "status")
  }
  dimnames(ss) <- list(NULL, cname)
  attr(ss, "type") <- type
  if (length(inputAttributes) > 0)
    attr(ss, "inputAttributes") <- inputAttributes
  class(ss) <- "Surv"
  ss
}

Event <- function(time, event) {
  if (missing(time))
    stop("A time argument is required")
  if (!is.numeric(time))
    stop("Time variable is not numeric")
  if (any(is.na(time)))
    warning("Invalid time variable. NA values included")
  #  if (any(time<0))
  #    warning("Invalid time variable. Non-negative values included")
  if (length(event) != length(time))
    stop("Time and event variables are different lengths")
  if (missing(event))
    stop("An event argument is required")
  if (is.numeric(event)) {
    if (any(is.na(event)))
      warning("Invalid event variable. NA values included")
    status <- event
  } else if (is.is.logical(event)) {
    if (any(is.na(event)))
      warning("Invalid event variable. NA values included")
    status <- as.numeric(event)
    warning("Event variable is logical, converted to numearic")
  } else if (is.factor(event)) {
    status <- as.numeric(as.factor(event)) - 1
    warning("Event variable is a factor, converted to numearic")
  } else stop("Invalid status value, must be logical or numeric")
  if (nlevels(as.factor(event)) > 3)
    warning("Event variable should not have more than three levels")
  ss <- cbind(time = time, status = status)
  type <- "right"

  inputAttributes <- list()
  if (!is.null(attributes(time)))
    inputAttributes$time <- attributes(time)
  cname <- dimnames(ss)[[2]]
  if (length(cname) == 0) {
    cname <- c("time", "status")
  }
  dimnames(ss) <- list(NULL, cname)
  attr(ss, "type") <- type
  if (length(inputAttributes) > 0)
    attr(ss, "inputAttributes") <- inputAttributes
  class(ss) <- "Event"
  ss
}

readSurv <- function(formula, data, weights, code.event, code.censoring, subset.condition, na.action) {
  data <- createAnalysisDataset(formula, data, weights, subset.condition, na.action)
  cl <- match.call()
  if (missing(formula))
    stop("A formula argument is required")
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "offset", "cluster")
  out_terms <- terms(formula, special, data = data)
  if (!is.null(attr(out_terms, "specials")$strata))
    stop("strata() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$offset))
    stop("offset() cannot appear in formula")
  if (!is.null(attr(out_terms, "specials")$cluster))
    stop("cluster() cannot appear in formula")
  mf$formula <- out_terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.extract(mf, "response")
  if (!inherits(Y, c("Event", "Surv"))) {
    stop("A 'Surv' or 'Event' object is expected")
  } else {
    t <- Y[, 1]
    if (any(t<0)) {
      stop("Invalid time variable. Expected non-negative values. ")
    }
    if (!all(Y[, 2] %in% c(code.event, code.censoring))) {
      stop("Invalid event codes. Must be 0 or 1 for survival and 0, 1 or 2 for competing risks, with 0 representing censoring, if event codes are not specified. ")
    } else {
      epsilon <- Y[, 2]
      d <- ifelse(Y[, 2] == code.censoring, 0, 1)
      d0 <- ifelse(Y[, 2] == code.censoring, 1, 0)
      d1 <- ifelse(Y[, 2] == code.event[1], 1, 0)
      d2 <- ifelse(Y[, 2] == code.event[2], 1, 0)
    }
  }
  if (is.na(all.vars(out_terms)[3])) {
    strata <- rep(1, nrow(data))
    strata_name <- NULL
  } else {
    strata_name <- all.vars(out_terms)[3]
    strata <- as.factor(data[[strata_name]])
  }
  if (is.null(weights)) {
    w <- rep(1, nrow(data))
  } else {
    w <- data[[weights]]
    if (!is.numeric(w))
      stop("Weights must be numeric")
    if (any(!is.finite(w)))
      stop("Weights must be finite")
    if (any(w < 0))
      stop("Weights must be non-negative")
    if (any(is.na(w)))
      stop("Weights contain NA values")
  }
  return(list(t = t, epsilon = epsilon, d = d, d0 = d0, d1 = d1, d2 = d2, strata = strata, strata_name = strata_name, w=w))
}

createAnalysisDataset <- function(formula, data, other.variables.analyzed=NULL, subset.condition=NULL, na.action=na.pass) {
  if (!is.null(subset.condition)) {
    analysis_dataset <- subset(data, eval(parse(text = subset.condition)))
  } else {
    analysis_dataset <- data
  }
  all_vars <- c(all.vars(formula), other.variables.analyzed)
  analysis_dataset <- analysis_dataset[, all_vars, drop = FALSE]
  return(na.action(analysis_dataset))
}

create_rr_text <- function(coefficient, cov, index, omit.conf.int=TRUE, conf.int=0.95) {
  alpha <- 1 - conf.int
  critical_value <- qnorm(1 - alpha / 2)
  coef <- coefficient[index]
  coef_se <- sqrt(diag(cov)[index])
  conf_low <- coef - critical_value * coef_se
  conf_high <- coef + critical_value * coef_se
  p_value <- floor(2 * (1 - pnorm(abs(coef) / coef_se)))
  if (omit.conf.int==TRUE) {
    if (p_value<0.01) text <- paste0("RR=", round(exp(coef), digit=2), ", p<0.01")
    else text <- paste0("RR=", round(exp(coef), digit=2), ", p=", p_value)
  } else {
    if (p_value<0.01) text <- paste0("RR=", round(exp(coef), digit=2), " (", round(exp(conf_low), digit=2), " to ", round(exp(conf_high), digit=2), ", p<0.01", ")")
    else text <- paste0("RR=", round(exp(coef), digit=2), " (", round(exp(conf_low), digit=2), " to ", round(exp(conf_high), digit=2), ", p=", p_value, ")")
  }
  return(text)
}

createTestData <- function(n, w, first_zero=FALSE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE) {
  one <- rep(1, n)
  t <- c(1:(n/2), 1:(n/2))
  epsilon <- rep(1, n)
  epsilon[2] <- 2
  epsilon[3] <- 2
  if (first_zero==TRUE) {
    epsilon[1] <- 0
    epsilon[n/2+1] <- 0
  }
  if (last_zero==TRUE) {
    epsilon[n/2] <- 0
    epsilon[n] <- 0
  }
  w <- rep(w, n)
  if (logical_strata==TRUE) {
    strata <- (t %% 2 == 0)
  } else {
    strata <- as.factor((t %% 2 == 0))
  }
  if (na_strata==TRUE) {
    strata[1] <- NA
  }
  subset <- rep(1, n)
  if (subset_present==TRUE) {
    subset[1] <- 0
  }
  d <- as.numeric(epsilon>0)
  return(data.frame(id = 1:n, t = t, epsilon = epsilon, d = d, w = w, strata = strata, subset=subset))
}

get_surv <- function(predicted.time, estimated.surv, estimated.time, predicted.strata=NULL, estimated.strata=NULL) {
  predicted.surv <- numeric(length(predicted.time))
  strata_start <- c(1, head(cumsum(estimated.strata), -1) + 1)
  strata_end <- cumsum(estimated.strata)
  if (any(is.na(predicted.time)))
    stop("Invalid predicted time variable. NA values included")

  for (i in seq_along(predicted.time)) {
    t <- predicted.time[i]
    strata_size <- estimated.strata[predicted.strata[i]]

    if (is.null(estimated.strata)|all(is.na(estimated.strata))) {
      #      time_until_t <- estimated.time[estimated.time <= t]
      time_until_t <- estimated.time[estimated.time < t]
      if (length(time_until_t) > 0) {
        time_index <- which.max(time_until_t)
        predicted.surv[i] <- estimated.surv[time_index]
      } else {
        predicted.surv[i] <- 1
      }
    } else if (strata_size > 0|all(is.na(estimated.strata))) {
      strata_indices <- strata_start[predicted.strata[i]]:strata_end[predicted.strata[i]]
      strata_time <- estimated.time[strata_indices]
      strata_surv <- estimated.surv[strata_indices]

      #      time_until_t <- strata_time[strata_time <= t]
      time_until_t <- strata_time[strata_time < t]

      if (length(time_until_t) > 0) {
        time_index <- which.max(time_until_t)
        predicted.surv[i] <- strata_surv[time_index]
      } else {
        predicted.surv[i] <- 1
      }
    } else {
      predicted.surv[i] <- NA
    }
  }
  return(predicted.surv)
}

calculateRMST <- function(out_survfit, tau){
  time <- out_survfit$time
  surv <- out_survfit$surv
  idx <- which(time <= tau)
  dt <- diff(c(0, time[idx], tau))
  surv_sub <- c(surv[idx], surv[max(idx)])
  rmst <- sum(surv_sub * dt)
  return(rmst)
}

calculateRMST <-function(out_survfit, tau, alpha = 0.05)
{
  idx = out_survfit$time <= tau
  wk.time = sort(c(out_survfit$time[idx], tau))
  wk.surv = out_survfit$surv[idx]
  wk.n.risk = out_survfit$n.risk[idx]
  wk.n.event = out_survfit$n.event[idx]
  time.diff <- diff(c(0, wk.time))
  areas <- time.diff * c(1, wk.surv)
  rmst = sum(areas)
  rmst
  wk.var <- ifelse((wk.n.risk - wk.n.event) == 0, 0, wk.n.event/(wk.n.risk*(wk.n.risk - wk.n.event)))
  wk.var = c(wk.var, 0)
  rmst.var = sum(cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
  rmst.se = sqrt(rmst.var)
  out = matrix(0, 2, 4)
  out[1, ] = c(rmst, rmst.se, rmst - qnorm(1 - alpha/2) * rmst.se,
               rmst + qnorm(1 - alpha/2) * rmst.se)
  out[2, ] = c(tau - out[1, 1], rmst.se, tau - out[1, 4], tau -
                 out[1, 3])
  rownames(out) = c("RMST", "RMTL")
  colnames(out) = c("Est.", "se", paste("lower .", round((1 - alpha) * 100, digits = 0), sep = ""), paste("upper .", round((1 - alpha) * 100, digits = 0), sep = ""))
  Z = list()
  Z$result = out
  Z$rmst = out[1, ]
  Z$rmtl = out[2, ]
  Z$tau = tau
  Z$rmst.var = rmst.var
  Z$fit = out_survfit
  return(Z)
}

check_survfit_object <- function(x, tol = 1e-10) {
  out <- list(ok = TRUE, errors = character(), warnings = character())
  add_err <- function(...)  { out$ok <<- FALSE; out$errors   <<- c(out$errors,   sprintf(...)) }
  add_warn <- function(...) {                    out$warnings <<- c(out$warnings, sprintf(...)) }

  ## 1) クラス
  if (!inherits(x, "survfit")) add_err("class(x) does not include 'survfit'")

  ## 2) 必須フィールドの存在
  req <- c("time","surv","n","n.risk","n.event","n.censor","std.err","lower","upper","conf.type")
  miss <- setdiff(req, names(x))
  if (length(miss)) add_err("missing fields: %s", paste(miss, collapse=", "))

  ## フィールドがあれば以降の整合チェック
  if (!length(miss)) {
    L <- length(x$time)
    same_len <- c(length(x$surv), length(x$n.risk), length(x$n.event),
                  length(x$n.censor), length(x$std.err), length(x$lower), length(x$upper))
    if (any(same_len != L))
      add_err("time/other vector lengths differ: time=%d vs others=%s",
              L, paste(same_len, collapse=","))

    ## 3) 数値範囲
    if (any(x$surv < -tol | x$surv > 1 + tol))
      add_err("surv values out of [0,1] range")

    if (x$conf.type %in% c("none","n")) {
      ## CI なしの場合は lower==upper==surv を期待
      if (any(abs(x$lower - x$surv) > 1e-8) || any(abs(x$upper - x$surv) > 1e-8))
        add_warn("conf.type is '%s' but lower/upper differ from surv", x$conf.type)
    } else {
      if (any(x$lower - x$surv > tol) || any(x$surv - x$upper > tol))
        add_err("CI envelope invalid: expected lower <= surv <= upper")
    }

    if (any(x$n.risk < -tol) || any(x$n.event < -tol) || any(x$n.censor < -tol))
      add_err("negative counts detected in n.risk/n.event/n.censor")

    ## 4) strata の整合
    if (!is.null(x$strata)) {
      if (!is.numeric(x$strata) && !is.integer(x$strata)) x$strata <- as.integer(x$strata)
      if (sum(x$strata) != L)
        add_err("sum(strata)=%d does not match length(time)=%d", sum(x$strata), L)
      # 各層ごとに time 単調増加、surv 単調非増加を検査
      brk <- c(0, cumsum(x$strata))
      for (s in seq_along(x$strata)) {
        idx <- (brk[s] + 1L):brk[s + 1L]
        if (any(diff(x$time[idx]) < -tol))
          add_err("time not non-decreasing within stratum %d", s)
        if (any(diff(x$surv[idx]) > tol))
          add_warn("surv not non-increasing within stratum %d (tolerance %.1e)", s, tol)
      }
    } else {
      if (any(diff(x$time) < -tol))  add_err("time not non-decreasing")
      if (any(diff(x$surv) > tol))   add_warn("surv not non-increasing (tolerance %.1e)", tol)
    }

    ## 5) n と n.risk のざっくり整合（層ごとの初期人数と合っているかの弱チェック）
    if (!is.null(x$strata)) {
      n0 <- x$n
      if (length(n0) != length(x$strata)) {
        add_warn("length(n) (%d) != number of strata (%d)", length(n0), length(x$strata))
      } else {
        brk <- c(0, cumsum(x$strata))
        for (s in seq_along(x$strata)) {
          idx <- (brk[s] + 1L):brk[s + 1L]
          if (length(idx) > 0 && abs(x$n.risk[idx[1]] - as.numeric(n0[s])) > 1e-6)
            add_warn("first n.risk in stratum %d (%g) != n[%d] (%g)", s, x$n.risk[idx[1]], s, as.numeric(n0[s]))
        }
      }
    } else if (length(x$n) == 1L && L > 0) {
      if (abs(x$n.risk[1] - as.numeric(x$n)) > 1e-6)
        add_warn("first n.risk (%g) != n (%g)", x$n.risk[1], as.numeric(x$n))
    }
  }

  if (out$ok) {
    message("✅ survfit compatibility: OK")
  } else {
    message("❌ survfit compatibility: FAILED")
  }
  if (length(out$warnings)) message("⚠️  Warnings:\n- ", paste(out$warnings, collapse = "\n- "))
  if (length(out$errors))   message("🛑 Errors:\n- ",   paste(out$errors,   collapse = "\n- "))
  invisible(out)
}
