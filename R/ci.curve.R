ci.curve <- function(formula, data, weights = NULL, subset = NULL,
                     code.event = c(1, 2), code.censoring = 0,
                     na.action = na.pass, conf.int = 0.95,
                     error = "delta",
                     conf.type = "arcsine-square root",
                     report.survfit.std.err = FALSE,
                     report.ggsurvfit = TRUE,
                     label.x = "Time",
                     label.y = "Cumulative incidence probability",
                     label.strata = NULL,
                     lims.x = NULL, lims.y = c(0, 1),
                     font.family = "sans", font.size = 14,
                     legend.position = "top") {

  checkDependentPackages()

  ### CHANGE: conf.type / error を厳密化（想定外を早期に弾く）
  conf.type <- match.arg(conf.type,
                         c("arcsine-square root","plain","log-log","none","n"))
  error     <- match.arg(error, c("delta","aalen"))

  out_readSurv <- readSurv(formula, data, weights, code.event, code.censoring, subset, na.action)
  out_aj <- calculateAJ(out_readSurv)
  out_aj <- readStrata(out_readSurv, out_aj, label.strata)

  ### CHANGE: strata 検出の判定を安定化
  str_id <- as.integer(out_readSurv$strata)
  has_strata <- (length(unique(str_id)) > 1)

  ## リスク集合サイズ n.risk（ストラタ別の初期人数から累積イベント/検閲を引く）
  if (has_strata) {
    # 各ストラタのサンプル数（1..max(str_id)の順に揃える）
    n_tab <- tabulate(str_id, nbins = max(str_id))
    stopifnot(length(n_tab) == length(out_aj$strata1))
    # strataごとの初期人数を、そのstrataのCIF点数だけ繰り返し
    n0_by_time <- rep.int(n_tab, out_aj$strata1)
    n.risk <- n0_by_time - out_aj$n.cum.censor - out_aj$n.cum.event1 - out_aj$n.cum.event2
  } else {
    n0 <- length(str_id)
    n.risk <- n0 - out_aj$n.cum.censor - out_aj$n.cum.event1 - out_aj$n.cum.event2
  }
  ### CHANGE: n.risk をガード（実装上の端数で負にならないように）
  n.risk <- pmax(as.numeric(n.risk), 0)

  ## ★★★ ストラタ別に SE を計算 ★★★
  if (has_strata) {
    se_vec <- numeric(length(out_aj$time1))

    # CIF 側のストラタ境界（各ストラタの CIF 点数）
    cif_len <- out_aj$strata1
    cif_breaks <- c(0, cumsum(cif_len))

    # KM 側のストラタ境界（各ストラタの KM ステップ数）
    km_len <- out_aj$km0_strata
    km_breaks <- c(0, cumsum(km_len))

    for (s in seq_along(cif_len)) {
      idx_cif <- (cif_breaks[s] + 1L):cif_breaks[s + 1L]
      idx_km  <- (km_breaks[s]  + 1L):km_breaks[s + 1L]
      se_vec[idx_cif] <- calculateAalenDeltaSE(
        CIF_time  = out_aj$time1[idx_cif],
        CIF_value = out_aj$aj1[idx_cif],
        n.event1  = out_aj$n.event1[idx_cif],
        n.event2  = out_aj$n.event2[idx_cif],
        n.atrisk  = n.risk[idx_cif],
        km_time   = out_aj$time0[idx_km],
        km_value  = out_aj$km0[idx_km],
        error     = error
      )
    }
    out_aj$std.err <- se_vec
  } else {
    out_aj$std.err <- calculateAalenDeltaSE(
      out_aj$time1, out_aj$aj1,
      out_aj$n.event1, out_aj$n.event2,
      n.risk,
      out_aj$time0, out_aj$km0,
      error
    )
  }

  ### CHANGE: surv 変換も数値ガード（丸め誤差で [0,1] 超え防止）
  out_aj$surv <- pmin(pmax(1 - out_aj$aj1, 0), 1)

  out_ci <- calculateCI(out_aj, conf.int, conf.type)

  ### CHANGE: surv でのSE報告時の 0 除算ガード
  if (report.survfit.std.err) {
    eps <- .Machine$double.eps^0.5
    out_aj$std.err <- out_aj$std.err / pmax(out_aj$surv, eps)
  }

  survfit_object <- list(
    time   = out_aj$time1,
    surv   = out_aj$surv,
    n      = if (has_strata) table(str_id) else length(str_id),
    n.risk = n.risk,
    n.event = out_aj$n.event1,
    n.censor = out_aj$n.censor,
    std.err = out_aj$std.err,
    upper = out_ci$upper,
    lower = out_ci$lower,
    conf.type = conf.type,
    call = match.call(),
    type = "Aalen-Johansen",
    method = "Aalen-Johansen"
  )
  if (has_strata) {
    survfit_object$strata <- out_aj$strata1
  }
  class(survfit_object) <- c("survfit")

  if (report.ggsurvfit) {
    if (is.null(lims.x)) lims.x <- c(0, max(out_readSurv$t))
    if (conf.type %in% c("none","n") || (has_strata && length(survfit_object$strata) > 2)) {
      survfit_object_ <- survfit_object
      survfit_object_$std.err <- rep(1, length(survfit_object_$surv))
      survfit_object_$lower <- survfit_object_$surv
      survfit_object_$upper <- survfit_object_$surv
      out_ggsurvfit <- ggsurvfit(survfit_object_, type = "risk") +
        theme_classic() +
        theme(legend.position = legend.position,
              axis.title = element_text(size = (font.size+2), family = font.family),
              axis.text = element_text(size = font.size, family = font.family),
              legend.text = element_text(size = font.size, family = font.family)) +
        labs(x = label.x, y = label.y) +
        lims(x = lims.x, y = lims.y) +
        theme(legend.position = "top") +
        add_risktable(risktable_stats = c("n.risk"))
    } else {
      if (length(survfit_object$time) != length(survfit_object$lower)) stop("time and upper/lower required for ggsurvfit are different lengths")
      if (length(survfit_object$time) != length(survfit_object$n.risk)) stop("time and n.risk used required ggsurvfit are different lengths")
      out_ggsurvfit <- ggsurvfit(survfit_object, type = "risk") +
        theme_classic() +
        theme(legend.position = legend.position,
              axis.title = element_text(size = (font.size+2), family = font.family),
              axis.text = element_text(size = font.size, family = font.family),
              legend.text = element_text(size = font.size, family = font.family)) +
        labs(x = label.x, y = label.y) +
        lims(x = lims.x, y = lims.y) +
        add_confidence_interval() +
        add_risktable(risktable_stats = c("n.risk"))
    }
    print(out_ggsurvfit)
  }
  return(survfit_object)
}

calculateAJ <- function(data) {
  out_km0 <- calculateKM_rcpp(data$t, data$d0, data$w, as.integer(data$strata), "none")
  # strataに応じたKMを取得
  km0_raw <- get_surv(data$t, out_km0$surv, out_km0$time, as.integer(data$strata), out_km0$strata)

  ### CHANGE: NA/Inf 防御（1/km0 のepsガード）
  eps <- .Machine$double.eps^0.5
  km0 <- pmin(pmax(km0_raw, eps), 1)  # [eps,1]に制限

  ip.weight <- (data$d0 == 0) * (1 / pmax(km0, eps))
  d1_ipw <- as.matrix(data$w * data$d1 * ip.weight)

  aj1 <- time1 <- n.cum.event1 <- n.cum.event2 <- n.cum.censor <- n.event1 <- n.event2 <- n.censor <- strata1 <- NULL
  for (level in sort(unique(as.integer(data$strata)))) {
    idx <- (as.integer(data$strata) == level)

    ### CHANGE: drop=FALSE で次元落ち防止
    sub_d1_ipw <- d1_ipw[idx, , drop = FALSE]

    sub_t  <- data$t[idx]
    sub_d0 <- data$d0[idx]
    sub_d1 <- data$d1[idx]
    sub_d2 <- data$d2[idx]

    not_atrisk <- outer(sub_t, sub_t, ">=")
    sub_aj1      <- not_atrisk %*% sub_d1_ipw / length(sub_t)
    sub_n.censor <- not_atrisk %*% as.matrix(sub_d0)
    sub_n.event1 <- not_atrisk %*% as.matrix(sub_d1)
    sub_n.event2 <- not_atrisk %*% as.matrix(sub_d2)

    is_unique <- !duplicated(sub_t)
    unique_t        <- sub_t[is_unique]
    unique_aj1      <- sub_aj1[is_unique, , drop = FALSE]  ### CHANGE: drop=FALSE
    unique_n.censor <- sub_n.censor[is_unique, , drop = FALSE]
    unique_n.event1 <- sub_n.event1[is_unique, , drop = FALSE]
    unique_n.event2 <- sub_n.event2[is_unique, , drop = FALSE]

    ord <- order(unique_t)
    sorted_t        <- unique_t[ord]
    sorted_aj1      <- unique_aj1[ord, , drop = FALSE]
    sorted_n.censor <- unique_n.censor[ord, , drop = FALSE]
    sorted_n.event1 <- unique_n.event1[ord, , drop = FALSE]
    sorted_n.event2 <- unique_n.event2[ord, , drop = FALSE]

    time1 <- c(time1, sorted_t)
    aj1   <- c(aj1,   as.numeric(sorted_aj1))

    # 直前までの累積を配置（最初の点の直前は0）
    if (nrow(sorted_n.censor) > 0) {
      n.cum.censor <- c(n.cum.censor, 0, cumsum(as.numeric(sorted_n.censor))[-nrow(sorted_n.censor)])
      n.cum.event1 <- c(n.cum.event1, 0, cumsum(as.numeric(sorted_n.event1))[-nrow(sorted_n.event1)])
      n.cum.event2 <- c(n.cum.event2, 0, cumsum(as.numeric(sorted_n.event2))[-nrow(sorted_n.event2)])
    } else {
      n.cum.censor <- c(n.cum.censor)
      n.cum.event1 <- c(n.cum.event1)
      n.cum.event2 <- c(n.cum.event2)
    }
    strata1 <- c(strata1, length(sorted_t))
  }

  # 差分で時点ごとの件数へ
  if (length(n.cum.event1) > 0) {
    n.censor <- n.event1 <- n.event2 <- numeric(length(n.cum.event1))
    n.censor[1] <- if (length(n.cum.censor)) n.cum.censor[1] else 0
    n.event1[1] <- if (length(n.cum.event1)) n.cum.event1[1] else 0
    n.event2[1] <- if (length(n.cum.event2)) n.cum.event2[1] else 0
    if (length(n.cum.event1) >= 2) {
      for (i in 2:length(n.cum.event1)) {
        n.censor[i] <- n.cum.censor[i] - n.cum.censor[i - 1L]
        n.event1[i] <- n.cum.event1[i] - n.cum.event1[i - 1L]
        n.event2[i] <- n.cum.event2[i] - n.cum.event2[i - 1L]
      }
    }
  } else {
    n.censor <- n.event1 <- n.event2 <- numeric(0)
  }

  list(time1 = time1, aj1 = aj1,
       n.event1 = n.event1, n.event2 = n.event2, n.censor = n.censor,
       n.cum.event1 = n.cum.event1, n.cum.event2 = n.cum.event2, n.cum.censor = n.cum.censor,
       strata1 = strata1,
       time0 = out_km0$time, km0 = km0,                 ### CHANGE: km0_raw → ガード済 km0 を返す
       km0_strata = out_km0$strata)
}

readStrata <- function(out_readSurv, out_aj, label.strata = NULL) {
  ### CHANGE: strata 検出の判定を安定化
  str_id <- as.integer(out_readSurv$strata)
  has_strata <- (length(unique(str_id)) > 1)

  if (has_strata && is.null(label.strata)) {
    names(out_aj$strata1) <- levels(as.factor(out_readSurv$strata))
  } else if (has_strata) {
    stopifnot(length(label.strata) == length(out_aj$strata1))
    names(out_aj$strata1) <- label.strata
  }
  return(out_aj)
}
