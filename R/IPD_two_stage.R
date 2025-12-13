#' Two-stage IPD meta-analysis with multiple imputation and forest plot
#'
#' @description
#' \code{IPD_two_stage()} fits separate regression models within each trial and
#' each imputed dataset, pools trial-specific coefficients across imputations
#' using Rubin's rules, then performs a random-effects meta-analysis across
#' trials (via \pkg{metafor}). Optionally, it computes pooled model performance
#' indices and produces a forest plot for a selected term.
#'
#' The function supports multiply imputed data (identified by \code{imp_col}),
#' clustered by trial (\code{intercept_var}), and several generalized linear
#' model families (negative binomial, Poisson, quasipoisson, Gaussian,
#' binomial). The forest plot is built with \pkg{ggplot2} and \pkg{cowplot}.
#'
#' NEW (modifications included):
#' - predictor_vars can include interaction shorthand using "*" (e.g. "x*y").
#'   These are expanded internally for formula building (so x, y, and x:y exist).
#' - Forest plot: CIs that cross the ref line (usually 1) are drawn in grey.
#' - Forest plot: if CI exceeds x_limits, the truncated side uses an arrowhead
#'   and the tick is NOT drawn on that side.
#'
#' @param data A \code{data.frame} containing the individual patient data
#'   (one row per subject per imputation).
#' @param outcome_var Character string; name of the outcome variable in
#'   \code{data}. Must be numeric for count outcomes or 0/1 for binomial.
#' @param predictor_vars Character vector of predictor terms to include in the
#'   model. Entries can be simple variable names present in \code{data} or
#'   formula-style terms such as interactions (e.g. \code{"x1:x2"},
#'   \code{"x1*x2"}).
#' @param covariables Optional character vector of additional covariate variable
#'   names to include in the model (all must be columns in \code{data}). These
#'   are checked for being "collapsed" (no variation) within trial/imputation
#'   and are dropped when necessary.
#' @param imp_col Character string indicating the name of the imputation
#'   indicator column in \code{data}. Each distinct value is treated as one
#'   imputed dataset.
#' @param followup_offset Character; either \code{"Yes"} or \code{"No"}. When
#'   \code{"Yes"}, an offset term \code{offset(log(followup_col))} is added to
#'   all models.
#' @param followup_col Character string; name of the follow-up time variable in
#'   \code{data} to be used in the offset (if \code{followup_offset = "Yes"}).
#'   Must be strictly positive where used.
#' @param intercept_var Character string; name of the trial identifier variable
#'   in \code{data}. Models are fit separately within each trial and imputation.
#' @param model_type Character; type of GLM to fit. One of
#'   \code{"nb"} (negative binomial via \code{MASS::glm.nb}),
#'   \code{"poisson"}, \code{"quasipoisson"}, \code{"gaussian"}, or
#'   \code{"binomial"}.
#' @param model_performance Logical; if \code{TRUE}, a global model including
#'   all trials (with trial fixed effects) is fit within each imputation to
#'   compute basic performance indices (AIC, BIC, RMSE, etc.), which are then
#'   pooled across imputations.
#' @param main_term Optional character string giving the name of the model term
#'   (as returned in \code{broom::tidy()$term}) to pool and display in the
#'   forest plot. If \code{NULL}, no pooling or forest plot is computed.
#' @param meta_method Character; method passed to \code{metafor::rma()} for the
#'   random-effects meta-analysis (default \code{"REML"}).
#' @param xlab Character; label for the x-axis of the forest plot.
#' @param ref_line Numeric; position of the vertical reference line in the
#'   forest plot (typically 1 for aRR / aOR).
#' @param x_limits Optional numeric vector of length 2 giving the x-axis limits
#'   (on the original scale) for the forest plot. If \code{NULL}, limits are
#'   derived from the data and possibly symmetrized on the log scale if
#'   \code{center_axis_at_one = TRUE}.
#' @param x_breaks Optional numeric vector of break points for the x-axis
#'   (original scale). If \code{NULL}, breaks are chosen automatically.
#' @param center_axis_at_one Logical; if \code{TRUE}, the x-axis is made
#'   symmetric around 1 on the log scale.
#' @param point_size_range Numeric vector of length 2 specifying the minimum and
#'   maximum point sizes for the forest plot (mapped to study weights).
#' @param show_size_legend Logical; if \code{TRUE}, a legend for point size
#'   (weights) is displayed in the forest plot.
#' @param show_headers Logical; if \code{TRUE}, column headers (e.g.
#'   \code{"aRR (95% CI)"}, \code{"N"}, \code{"Weight"}) are displayed above
#'   the right-hand text columns.
#' @param diamond_color Character; color for the pooled effect point symbol in
#'   the forest plot.
#' @param pooled_color Character; color for the pooled row label and text
#'   columns (and pooled CI line/arrow).
#' @param bold_pooled Logical; if \code{TRUE}, the pooled row label is drawn in
#'   bold.
#' @param margin_plot Numeric vector of length 4 giving the plot margins (top,
#'   right, bottom, left) for the main forest panel, passed to
#'   \code{ggplot2::margin()}.
#'
#' @export
IPD_two_stage <- function(data,
                          outcome_var,
                          predictor_vars,
                          covariables = NULL,
                          imp_col = ".imp",
                          followup_offset = c("Yes", "No"),
                          followup_col = NULL,
                          intercept_var,
                          model_type = c("nb", "poisson", "quasipoisson",
                                         "gaussian", "binomial"),
                          model_performance = FALSE,
                          ## ---- options for pooled main term & forest plot ----
                          main_term          = NULL,
                          meta_method        = "REML",
                          xlab               = "aRR",
                          ref_line           = 1,
                          x_limits           = NULL,
                          x_breaks           = NULL,
                          center_axis_at_one = TRUE,
                          point_size_range   = c(2, 6),
                          show_size_legend   = FALSE,
                          show_headers       = TRUE,
                          diamond_color      = "blue",
                          pooled_color       = "black",
                          bold_pooled        = FALSE,
                          margin_plot        = c(5.5, 2, 5.5, 5.5)
) {

  ## ------------------------------------------------------------------
  ## Internal helpers
  ## ------------------------------------------------------------------

  is_collapsed <- function(vec) {
    x <- vec[!is.na(vec)]
    if (length(x) == 0L) TRUE else length(unique(x)) < 2L
  }

  auto_log_breaks <- function(lims) {
    lb <- lims[1]; ub <- lims[2]
    std <- c(0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9,
             1, 1.1, 1.25, 1.4, 1.5, 1.6, 1.75,
             2, 2.5, 3, 4, 5, 7, 10, 15, 20)
    br <- std[std >= lb & std <= ub]
    if (length(br) == 0) br <- exp(pretty(log(c(lb, ub)), n = 5))
    if (!any(abs(br - 1) < 1e-8)) br <- sort(unique(c(br, 1)))
    br
  }

  compute_centered_limits <- function(lims, from_data_lower, from_data_upper, center = TRUE) {
    if (!center) {
      if (is.null(lims)) return(c(from_data_lower, from_data_upper))
      return(lims)
    }
    if (is.null(lims)) {
      log_lower <- log(from_data_lower); log_upper <- log(from_data_upper)
    } else {
      stopifnot(length(lims) == 2, all(is.finite(lims)), all(lims > 0))
      log_lower <- log(lims[1]); log_upper <- log(lims[2])
    }
    range_max <- max(abs(log_lower), abs(log_upper))
    exp(c(-range_max, range_max))
  }

  ## --- NEW: expand predictor_vars that include "*" into the full set of terms ---
  expand_star_terms <- function(term_vec) {
    out <- character(0)

    for (t in term_vec) {
      if (grepl("\\*", t) && !grepl("I\\(", t)) {
        # build terms(~ <t>) and extract labels (drops intercept)
        f <- stats::as.formula(paste0("~", t))
        tl <- attr(stats::terms(f), "term.labels")
        out <- c(out, tl)
      } else {
        out <- c(out, t)
      }
    }

    unique(out)
  }

  predictor_vars_expanded <- expand_star_terms(predictor_vars)

  pool_IPD_two_stage_inner <- function(coefs_per_trial, term, method = "REML") {
    coefs_term <- coefs_per_trial %>%
      dplyr::filter(.data$term == !!term) %>%
      dplyr::filter(!is.na(.data$trial), .data$trial != "")

    if (nrow(coefs_term) == 0L) {
      stop("No coefficients found for term = ", term,
           ". If this is a factor, you may need the level-specific name (e.g. varYes).")
    }

    meta_per_imp <- list()
    for (imp_val in sort(unique(coefs_term$imp))) {
      dat_imp <- coefs_term %>%
        dplyr::filter(.data$imp == imp_val) %>%
        dplyr::select(.data$trial, .data$estimate, .data$std.error)

      meta_fit <- metafor::rma(
        yi     = dat_imp$estimate,
        sei    = dat_imp$std.error,
        method = method
      )
      meta_per_imp[[paste0("imp_", imp_val)]] <- meta_fit
    }

    Q_i <- vapply(meta_per_imp, function(x) x$b[1, 1], numeric(1))
    U_i <- vapply(meta_per_imp, function(x) x$se^2,   numeric(1))
    m   <- length(Q_i)

    Q_bar <- mean(Q_i)
    U_bar <- mean(U_i)
    B     <- stats::var(Q_i)
    if (!is.finite(B)) B <- 0
    T_var <- U_bar + (1 + 1 / m) * B
    se_pooled <- sqrt(T_var)

    df <- (m - 1) * (1 + U_bar / ((1 + 1/m) * B))^2
    if (!is.finite(df)) df <- m - 1

    t_val <- Q_bar / se_pooled
    p_val <- 2 * stats::pt(abs(t_val), df = df, lower.tail = FALSE)
    ci_log <- Q_bar + c(-1, 1) * stats::qt(0.975, df = df) * se_pooled

    list(
      term                = term,
      method              = method,
      m_imputations       = m,
      pooled_log_estimate = Q_bar,
      pooled_log_se       = se_pooled,
      df                  = df,
      p_value             = p_val,
      ci_log              = ci_log,
      pooled_estimate_exp = exp(Q_bar),
      ci_exp              = exp(ci_log)
    )
  }

  pool_IPD_two_stage_performance_inner <- function(perf_list) {
    perf_df <- dplyr::bind_rows(
      lapply(names(perf_list), function(nm) {
        if (is.null(perf_list[[nm]])) return(NULL)
        tibble::tibble(
          imp    = as.integer(gsub("imp_", "", nm)),
          AIC    = perf_list[[nm]]$AIC,
          BIC    = perf_list[[nm]]$BIC,
          logLik = perf_list[[nm]]$logLik,
          RMSE   = perf_list[[nm]]$RMSE,
          Sigma  = perf_list[[nm]]$Sigma,
          Score_log = perf_list[[nm]]$Score_log
        )
      })
    )

    pooled <- perf_df %>%
      dplyr::summarise(
        m = dplyr::n(),
        AIC_mean  = mean(.data$AIC, na.rm = TRUE),
        AIC_sd    = stats::sd(.data$AIC, na.rm = TRUE),
        BIC_mean  = mean(.data$BIC, na.rm = TRUE),
        BIC_sd    = stats::sd(.data$BIC, na.rm = TRUE),
        logLik_mean = mean(.data$logLik, na.rm = TRUE),
        logLik_sd   = stats::sd(.data$logLik, na.rm = TRUE),
        RMSE_mean   = mean(.data$RMSE, na.rm = TRUE),
        RMSE_sd     = stats::sd(.data$RMSE, na.rm = TRUE),
        Sigma_mean  = mean(.data$Sigma, na.rm = TRUE),
        Sigma_sd    = stats::sd(.data$Sigma, na.rm = TRUE),
        Score_log_mean = mean(.data$Score_log, na.rm = TRUE),
        Score_log_sd   = stats::sd(.data$Score_log, na.rm = TRUE)
      )

    list(per_imputation = perf_df, pooled_summary = pooled)
  }

  forest_IPD_two_stage_inner <- function(coefs_per_trial,
                                         fit_log_tbl,
                                         trial_levels,
                                         term,
                                         xlab,
                                         ref_line,
                                         x_limits,
                                         x_breaks,
                                         center_axis_at_one,
                                         point_size_range,
                                         show_size_legend,
                                         show_headers,
                                         diamond_color,
                                         pooled_color,
                                         bold_pooled,
                                         meta_method,
                                         margin_plot) {

    coefs_term <- coefs_per_trial %>%
      dplyr::filter(.data$term == !!term) %>%
      dplyr::filter(!is.na(.data$trial), .data$trial != "")

    if (nrow(coefs_term) == 0L) stop("No coefficients found for term = ", term)

    ## Rubin pooling per trial across imputations
    trial_effects <- coefs_term %>%
      dplyr::group_by(.data$trial) %>%
      dplyr::summarise(
        M_imp   = dplyr::n(),
        Q_bar   = mean(.data$estimate, na.rm = TRUE),
        U_bar   = mean(.data$std.error^2, na.rm = TRUE),
        B_raw   = stats::var(.data$estimate, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        B     = ifelse(.data$M_imp > 1, .data$B_raw, 0),
        T_var = .data$U_bar + (1 + 1 / .data$M_imp) * .data$B,
        SE_mi = sqrt(.data$T_var)
      ) %>%
      dplyr::filter(is.finite(.data$SE_mi), .data$SE_mi > 0)

    if (nrow(trial_effects) == 0L) stop("All trial-specific SEs are non-finite for term = ", term)

    ## N per trial
    N_by_trial <- fit_log_tbl %>%
      dplyr::group_by(.data$trial) %>%
      dplyr::summarise(N = suppressWarnings(max(.data$n, na.rm = TRUE)), .groups = "drop")

    trial_effects <- trial_effects %>%
      dplyr::left_join(N_by_trial, by = "trial")

    ## meta across trials
    meta_fit <- metafor::rma(
      yi     = trial_effects$Q_bar,
      sei    = trial_effects$SE_mi,
      method = meta_method
    )

    pooled_log <- as.numeric(meta_fit$b)
    pooled_se  <- as.numeric(meta_fit$se)
    pooled_ci  <- c(pooled_log - 1.96 * pooled_se, pooled_log + 1.96 * pooled_se)

    wt <- tryCatch(as.numeric(stats::weights(meta_fit)), error = function(e) numeric(0))
    if (length(wt) == nrow(trial_effects) && all(is.finite(wt))) {
      Weight_percent <- 100 * wt / sum(wt)
    } else {
      wt2 <- 1 / (trial_effects$SE_mi^2)
      Weight_percent <- 100 * wt2 / sum(wt2)
    }

    format_num <- function(x) ifelse(x >= 10, sprintf("%.0f", x), sprintf("%.1f", x))

    df_plot <- trial_effects %>%
      dplyr::mutate(
        Trial  = as.character(.data$trial),
        RR     = exp(.data$Q_bar),
        RR_lo  = exp(.data$Q_bar - 1.96 * .data$SE_mi),
        RR_hi  = exp(.data$Q_bar + 1.96 * .data$SE_mi),
        Weight_percent = Weight_percent,
        aRR_CI = sprintf("%.2f [%.2f–%.2f]", .data$RR, .data$RR_lo, .data$RR_hi),
        N_col  = ifelse(is.finite(.data$N), format(.data$N, big.mark = ","), ""),
        W_col  = paste0(format_num(.data$Weight_percent), "%")
      ) %>%
      dplyr::select(.data$Trial, .data$RR, .data$RR_lo, .data$RR_hi,
                    .data$aRR_CI, .data$N_col, .data$W_col,
                    .data$Weight_percent)

    ## pooled row
    pooled_row <- df_plot[1, , drop = FALSE]
    pooled_row[] <- NA
    pooled_row$Trial          <- "Pooled (REML)"
    pooled_row$RR             <- exp(pooled_log)
    pooled_row$RR_lo          <- exp(pooled_ci[1])
    pooled_row$RR_hi          <- exp(pooled_ci[2])
    pooled_row$aRR_CI         <- sprintf("%.2f [%.2f–%.2f]",
                                         exp(pooled_log), exp(pooled_ci[1]), exp(pooled_ci[2]))
    pooled_row$N_col          <- format(sum(trial_effects$N, na.rm = TRUE), big.mark = ",")
    pooled_row$W_col          <- ""
    pooled_row$Weight_percent <- max(Weight_percent, na.rm = TRUE)

    df_plot <- dplyr::bind_rows(df_plot, pooled_row)

    ## x limits
    if (is.null(x_limits)) {
      lower <- min(df_plot$RR_lo, na.rm = TRUE)
      upper <- max(df_plot$RR_hi, na.rm = TRUE)
    } else {
      lower <- x_limits[1]; upper <- x_limits[2]
    }

    final_limits <- compute_centered_limits(
      lims            = c(lower, upper),
      from_data_lower = min(df_plot$RR_lo, na.rm = TRUE),
      from_data_upper = max(df_plot$RR_hi, na.rm = TRUE),
      center          = center_axis_at_one
    )

    ## truncation flags + plotting values
    df_plot <- df_plot %>%
      dplyr::mutate(
        lo_trunc   = .data$RR_lo < final_limits[1],
        hi_trunc   = .data$RR_hi > final_limits[2],
        RR_lo_plot = pmax(.data$RR_lo, final_limits[1]),
        RR_hi_plot = pmin(.data$RR_hi, final_limits[2]),
        RR_plot    = pmin(pmax(.data$RR, final_limits[1]), final_limits[2]),
        is_pooled  = (.data$Trial == "Pooled (REML)")
      )

    ## breaks
    if (!is.null(x_breaks)) {
      x_breaks_final <- x_breaks[x_breaks >= final_limits[1] & x_breaks <= final_limits[2]]
      if (length(x_breaks_final) == 0) x_breaks_final <- auto_log_breaks(final_limits)
    } else {
      x_breaks_final <- auto_log_breaks(final_limits)
    }

    ## y order
    trials_ord <- trial_levels[trial_levels %in% df_plot$Trial]
    levels_vec <- c(trials_ord, "Pooled (REML)")
    levels_vec <- levels_vec[!duplicated(levels_vec)]
    y_levels_rev <- rev(levels_vec)
    df_plot$Trial <- factor(df_plot$Trial, levels = y_levels_rev)

    if (is.null(margin_plot) || length(margin_plot) != 4L) {
      margin_plot <- c(5.5, 2, 5.5, 5.5)
    }

    ## ---- NEW: CI significance colouring + arrows without ticks ----
    df_plot <- df_plot %>%
      dplyr::mutate(
        crosses_ref = (.data$RR_lo <= ref_line & .data$RR_hi >= ref_line),
        ci_col = dplyr::case_when(
          .data$is_pooled ~ pooled_color,
          .data$crosses_ref ~ "grey60",
          TRUE ~ "black"
        )
      )

    arrow_len <- grid::unit(0.18, "cm")
    tick_h    <- 0.20

    p_forest <- ggplot2::ggplot(df_plot, ggplot2::aes(y = .data$Trial)) +
      # CI segment (already truncated)
      ggplot2::geom_segment(
        ggplot2::aes(
          x    = .data$RR_lo_plot,
          xend = .data$RR_hi_plot,
          yend = .data$Trial,
          colour = .data$ci_col
        ),
        linewidth = 0.6,
        lineend   = "butt"
      ) +
      # left tick if NOT truncated on left
      ggplot2::geom_segment(
        data = df_plot %>% dplyr::filter(!.data$lo_trunc),
        ggplot2::aes(
          x    = .data$RR_lo_plot,
          xend = .data$RR_lo_plot,
          y    = as.numeric(.data$Trial) - tick_h/2,
          yend = as.numeric(.data$Trial) + tick_h/2,
          colour = .data$ci_col
        ),
        inherit.aes = FALSE,
        linewidth = 0.6,
        lineend   = "butt"
      ) +
      # right tick if NOT truncated on right
      ggplot2::geom_segment(
        data = df_plot %>% dplyr::filter(!.data$hi_trunc),
        ggplot2::aes(
          x    = .data$RR_hi_plot,
          xend = .data$RR_hi_plot,
          y    = as.numeric(.data$Trial) - tick_h/2,
          yend = as.numeric(.data$Trial) + tick_h/2,
          colour = .data$ci_col
        ),
        inherit.aes = FALSE,
        linewidth = 0.6,
        lineend   = "butt"
      ) +
      # right arrow if truncated (no tick)
      ggplot2::geom_segment(
        data = df_plot %>% dplyr::filter(.data$hi_trunc),
        ggplot2::aes(
          x    = final_limits[2] / 1.12,
          xend = final_limits[2],
          y    = .data$Trial,
          yend = .data$Trial,
          colour = .data$ci_col
        ),
        inherit.aes = FALSE,
        arrow = ggplot2::arrow(type = "closed", length = arrow_len),
        linewidth = 0.6,
        lineend   = "butt"
      ) +
      # left arrow if truncated (no tick)
      ggplot2::geom_segment(
        data = df_plot %>% dplyr::filter(.data$lo_trunc),
        ggplot2::aes(
          x    = final_limits[1] * 1.12,
          xend = final_limits[1],
          y    = .data$Trial,
          yend = .data$Trial,
          colour = .data$ci_col
        ),
        inherit.aes = FALSE,
        arrow = ggplot2::arrow(type = "closed", length = arrow_len),
        linewidth = 0.6,
        lineend   = "butt"
      ) +
      # points (kept default; sizes reflect weights)
      ggplot2::geom_point(
        ggplot2::aes(x = .data$RR_plot, size = .data$Weight_percent, shape = .data$is_pooled),
        show.legend = show_size_legend
      ) +
      ggplot2::scale_size_continuous(range = point_size_range, name = "Weight (%)") +
      ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
      # pooled diamond (explicit colour)
      ggplot2::geom_point(
        data = df_plot %>% dplyr::filter(.data$is_pooled),
        ggplot2::aes(y = .data$Trial, x = .data$RR_plot),
        inherit.aes = FALSE,
        colour = diamond_color,
        shape  = 18,
        size   = max(point_size_range) * 1.1
      ) +
      ggplot2::geom_vline(xintercept = ref_line, linetype = "dashed", linewidth = 0.5) +
      ggplot2::scale_colour_identity() +
      ggplot2::scale_x_log10(breaks = x_breaks_final, limits = final_limits) +
      ggplot2::scale_y_discrete(limits = y_levels_rev) +
      ggplot2::labs(x = xlab, y = NULL) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.margin      = ggplot2::margin(margin_plot[1], margin_plot[2],
                                           margin_plot[3], margin_plot[4]),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.y      = ggplot2::element_blank(),
        axis.ticks.y     = ggplot2::element_blank()
      )

    ## left labels
    df_left <- data.frame(
      Trial = df_plot$Trial,
      lab   = as.character(df_plot$Trial),
      col   = ifelse(as.character(df_plot$Trial) == "Pooled (REML)", pooled_color, "black"),
      face  = ifelse(as.character(df_plot$Trial) == "Pooled (REML)" & bold_pooled, "bold", "plain")
    )

    p_left <- ggplot2::ggplot(
      df_left,
      ggplot2::aes(y = .data$Trial, x = 1, label = .data$lab, colour = .data$col)
    ) +
      ggplot2::geom_text(hjust = 1, size = 3.2, fontface = df_left$face) +
      ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0.02, 0)) +
      ggplot2::scale_colour_identity() +
      ggplot2::labs(x = NULL, y = NULL, title = "Trials") +
      ggplot2::theme_void(base_size = 12) +
      ggplot2::theme(
        plot.margin = ggplot2::margin(5.5, 0, 5.5, 5.5),
        plot.title  = ggplot2::element_text(hjust = 1, size = 11, face = "bold",
                                            margin = ggplot2::margin(b = 4))
      )

    make_text_col <- function(labels, title, color_vec, x_pos = 0.02) {
      df_txt <- data.frame(Trial = df_plot$Trial, lab = labels, col = color_vec)
      ggplot2::ggplot(df_txt,
                      ggplot2::aes(y = .data$Trial, x = x_pos, label = .data$lab, colour = .data$col)) +
        ggplot2::geom_text(hjust = 0, size = 3.2) +
        ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0.02, 0)) +
        ggplot2::scale_colour_identity() +
        ggplot2::labs(x = NULL, y = NULL, title = if (show_headers) title else NULL) +
        ggplot2::theme_void(base_size = 12) +
        ggplot2::theme(
          plot.margin = ggplot2::margin(5.5, 1, 5.5, 0),
          plot.title  = ggplot2::element_text(hjust = 0, size = 11, face = "bold",
                                              margin = ggplot2::margin(b = 4)),
          axis.text.y = ggplot2::element_blank()
        ) +
        ggplot2::coord_cartesian(clip = "off")
    }

    pooled_mask <- as.character(df_plot$Trial) == "Pooled (REML)"
    col_vec <- ifelse(pooled_mask, pooled_color, "black")

    p_col1 <- make_text_col(df_plot$aRR_CI, "aRR (95% CI)", col_vec)
    p_col2 <- make_text_col(df_plot$N_col,  "N",            col_vec)
    p_col3 <- make_text_col(df_plot$W_col,  "Weight",       col_vec)

    full_plot <- cowplot::plot_grid(
      p_left, p_forest, p_col1, p_col2, p_col3,
      nrow       = 1,
      rel_widths = c(0.9, 2.2, 0.6, 0.25, 0.5),
      align      = "hv",
      axis       = "tb"
    )

    heterogeneity_df <- data.frame(
      tau2 = meta_fit$tau2,
      I2   = meta_fit$I2,
      H2   = meta_fit$H2,
      Q    = meta_fit$QE,
      df_Q = meta_fit$k - 1,
      k    = meta_fit$k
    )

    list(plot = full_plot, pooled_log = pooled_log, pooled_ci = pooled_ci,
         heterogeneity = heterogeneity_df)
  }

  ## ------------------------------------------------------------------
  ## 1) Basic checks & set-up
  ## ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  if (!imp_col %in% names(data)) stop("imp_col '", imp_col, "' not found in data.")
  if (!outcome_var %in% names(data)) stop("outcome_var '", outcome_var, "' not found in data.")
  if (!intercept_var %in% names(data)) stop("intercept_var '", intercept_var, "' not found in data.")

  model_type      <- match.arg(model_type)
  followup_offset <- match.arg(followup_offset)

  if (followup_offset == "Yes") {
    if (is.null(followup_col) || !followup_col %in% names(data)) {
      stop("followup_offset = 'Yes' but followup_col is missing or not in data.")
    }
  }

  all_covars <- if (is.null(covariables)) character(0) else covariables

  rhs_template <- c(predictor_vars, all_covars)
  rhs_string   <- paste(rhs_template, collapse = " + ")
  if (followup_offset == "Yes") rhs_string <- paste0(rhs_string, " + offset(log(", followup_col, "))")
  form_global <- stats::as.formula(paste0(outcome_var, " ~ ", rhs_string))

  trial_levels <- if (is.factor(data[[intercept_var]])) levels(data[[intercept_var]])
  else unique(as.character(data[[intercept_var]]))

  imp_values   <- sort(unique(data[[imp_col]]))
  trial_values <- trial_levels

  coefs_list     <- list()
  models_per_imp <- vector("list", length(imp_values))
  names(models_per_imp) <- paste0("imp_", imp_values)

  fit_log   <- list()
  log_index <- 0L

  fit_one_model <- function(dat_trial, model_type, form) {
    w_msgs <- character(0)
    res <- tryCatch({
      m <- withCallingHandlers({
        if (model_type == "nb") {
          MASS::glm.nb(formula = form, data = dat_trial)
        } else if (model_type == "poisson") {
          stats::glm(formula = form, data = dat_trial, family = stats::poisson())
        } else if (model_type == "quasipoisson") {
          stats::glm(formula = form, data = dat_trial, family = stats::quasipoisson())
        } else if (model_type == "gaussian") {
          stats::glm(formula = form, data = dat_trial, family = stats::gaussian())
        } else if (model_type == "binomial") {
          stats::glm(formula = form, data = dat_trial, family = stats::binomial())
        } else {
          stop("Unsupported model_type: ", model_type)
        }
      }, warning = function(w) {
        w_msgs <<- c(w_msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
      list(ok = TRUE, model = m, error = NULL, warning = w_msgs)
    }, error = function(e) {
      list(ok = FALSE, model = NULL, error = conditionMessage(e), warning = w_msgs)
    })
    res
  }

  ## ------------------------------------------------------------------
  ## 2) Loop over imputations & trials
  ## ------------------------------------------------------------------
  for (imp in imp_values) {
    dat_imp <- data[data[[imp_col]] == imp, , drop = FALSE]

    models_this_imp <- vector("list", length(trial_values))
    names(models_this_imp) <- trial_values

    for (trial in trial_values) {
      dat_trial <- dat_imp[as.character(dat_imp[[intercept_var]]) == trial, , drop = FALSE]

      if (nrow(dat_trial) == 0L || length(unique(dat_trial[[outcome_var]])) < 2L) {
        log_index <- log_index + 1L
        fit_log[[log_index]] <- tibble::tibble(
          trial          = as.character(trial),
          imp            = imp,
          n              = nrow(dat_trial),
          ok             = FALSE,
          error          = "No data or no variation in outcome for this trial/imp.",
          warning        = NA_character_,
          dropped_covars = NA_character_
        )
        next
      }

      collapsed_covars <- all_covars[
        vapply(all_covars, function(v) {
          if (v %in% names(dat_trial)) is_collapsed(dat_trial[[v]]) else FALSE
        }, logical(1))
      ]

      # check collapse on expanded main effects when possible (columns only)
      predictor_collapsed <- vapply(predictor_vars_expanded, function(v) {
        if (v %in% names(dat_trial)) is_collapsed(dat_trial[[v]]) else FALSE
      }, logical(1))

      if (any(predictor_collapsed)) {
        log_index <- log_index + 1L
        fit_log[[log_index]] <- tibble::tibble(
          trial          = as.character(trial),
          imp            = imp,
          n              = nrow(dat_trial),
          ok             = FALSE,
          error          = paste0(
            "Predictor(s) with no variation in this trial/imp: ",
            paste(predictor_vars_expanded[predictor_collapsed], collapse = ", ")
          ),
          warning        = NA_character_,
          dropped_covars = if (length(collapsed_covars) > 0L) paste(collapsed_covars, collapse = "; ")
          else NA_character_
        )
        next
      }

      rhs_terms_trial <- c(predictor_vars, setdiff(all_covars, collapsed_covars))
      rhs_trial <- paste(rhs_terms_trial, collapse = " + ")
      if (followup_offset == "Yes") rhs_trial <- paste0(rhs_trial, " + offset(log(", followup_col, "))")

      form_trial <- stats::as.formula(paste0(outcome_var, " ~ ", rhs_trial))
      fit_res <- fit_one_model(dat_trial, model_type, form_trial)

      log_index <- log_index + 1L
      fit_log[[log_index]] <- tibble::tibble(
        trial          = as.character(trial),
        imp            = imp,
        n              = nrow(dat_trial),
        ok             = fit_res$ok,
        error          = if (!is.null(fit_res$error)) fit_res$error else NA_character_,
        warning        = if (length(fit_res$warning) > 0L) paste(unique(fit_res$warning), collapse = " | ")
        else NA_character_,
        dropped_covars = if (length(collapsed_covars) > 0L) paste(collapsed_covars, collapse = "; ")
        else NA_character_
      )

      if (!fit_res$ok || is.null(fit_res$model)) next

      models_this_imp[[as.character(trial)]] <- fit_res$model

      tt <- broom::tidy(fit_res$model) %>%
        dplyr::mutate(trial = as.character(trial), imp = imp, .before = 1)

      coefs_list[[length(coefs_list) + 1L]] <- tt
    }

    models_per_imp[[paste0("imp_", imp)]] <- models_this_imp
  }

  coefs_per_trial <- dplyr::bind_rows(coefs_list)
  fit_log_tbl     <- dplyr::bind_rows(fit_log)

  ## ------------------------------------------------------------------
  ## 3) Performance per imputation (optional)
  ## ------------------------------------------------------------------
  performance_per_imp <- NULL
  if (isTRUE(model_performance)) {
    perf_list <- vector("list", length(imp_values))
    names(perf_list) <- paste0("imp_", imp_values)

    for (imp in imp_values) {
      dat_imp <- data[data[[imp_col]] == imp, , drop = FALSE]

      collapsed_covars_imp <- all_covars[
        vapply(all_covars, function(v) {
          if (v %in% names(dat_imp)) is_collapsed(dat_imp[[v]]) else FALSE
        }, logical(1))
      ]

      rhs_perf_terms <- c(
        predictor_vars,
        setdiff(all_covars, collapsed_covars_imp),
        paste0("as.factor(", intercept_var, ")")
      )
      rhs_perf <- paste(rhs_perf_terms, collapse = " + ")
      if (followup_offset == "Yes") rhs_perf <- paste0(rhs_perf, " + offset(log(", followup_col, "))")

      form_perf <- stats::as.formula(paste0(outcome_var, " ~ ", rhs_perf))

      perf_fit <- tryCatch({
        if (model_type == "nb") {
          MASS::glm.nb(formula = form_perf, data = dat_imp)
        } else if (model_type == "poisson") {
          stats::glm(formula = form_perf, data = dat_imp, family = stats::poisson())
        } else if (model_type == "quasipoisson") {
          stats::glm(formula = form_perf, data = dat_imp, family = stats::quasipoisson())
        } else if (model_type == "gaussian") {
          stats::glm(formula = form_perf, data = dat_imp, family = stats::gaussian())
        } else if (model_type == "binomial") {
          stats::glm(formula = form_perf, data = dat_imp, family = stats::binomial())
        } else {
          stop("Unsupported model_type for performance.")
        }
      }, error = function(e) NULL)

      if (!is.null(perf_fit)) {
        y_obs  <- dat_imp[[outcome_var]]
        mu_hat <- stats::predict(perf_fit, type = "response")
        rmse <- sqrt(mean((y_obs - mu_hat)^2, na.rm = TRUE))

        sigma <- if (inherits(perf_fit, "negbin")) sqrt(1) else sqrt(summary(perf_fit)$dispersion)
        loglik    <- as.numeric(stats::logLik(perf_fit))
        n_obs     <- length(y_obs)
        score_log <- loglik / n_obs

        perf_list[[paste0("imp_", imp)]] <- list(
          AIC    = stats::AIC(perf_fit),
          BIC    = stats::BIC(perf_fit),
          logLik = loglik,
          RMSE   = rmse,
          Sigma  = sigma,
          Score_log = score_log
        )
      } else {
        perf_list[[paste0("imp_", imp)]] <- NULL
      }
    }

    performance_per_imp <- perf_list
  }

  pooled_perf <- if (!is.null(performance_per_imp)) pool_IPD_two_stage_performance_inner(performance_per_imp) else NULL

  ## ------------------------------------------------------------------
  ## 4) Final object
  ## ------------------------------------------------------------------
  out <- list(
    call                = match.call(),
    formula             = form_global,
    outcome_var         = outcome_var,
    predictor_vars      = predictor_vars,
    covariables         = all_covars,
    intercept_var       = intercept_var,
    imp_col             = imp_col,
    model_type          = model_type,
    followup_offset     = followup_offset,
    followup_col        = followup_col,
    coefs_per_trial     = coefs_per_trial,
    models_per_imp      = models_per_imp,
    performance_per_imp = performance_per_imp,
    fit_log             = fit_log_tbl,
    trial_levels        = trial_levels,
    Pooled_performance  = pooled_perf,
    perf_pool           = pooled_perf
  )

  ## ------------------------------------------------------------------
  ## 5) Optional: pooled main term + forest plot
  ## ------------------------------------------------------------------
  if (!is.null(main_term)) {
    out$pool_main <- pool_IPD_two_stage_inner(coefs_per_trial, main_term, meta_method)

    forest_res <- forest_IPD_two_stage_inner(
      coefs_per_trial    = coefs_per_trial,
      fit_log_tbl        = fit_log_tbl,
      trial_levels       = trial_levels,
      term               = main_term,
      xlab               = xlab,
      ref_line           = ref_line,
      x_limits           = x_limits,
      x_breaks           = x_breaks,
      center_axis_at_one = center_axis_at_one,
      point_size_range   = point_size_range,
      show_size_legend   = show_size_legend,
      show_headers       = show_headers,
      diamond_color      = diamond_color,
      pooled_color       = pooled_color,
      bold_pooled        = bold_pooled,
      meta_method        = meta_method,
      margin_plot        = margin_plot
    )

    out$forest_plot        <- forest_res$plot
    out$pooled_log_rr      <- forest_res$pooled_log
    out$pooled_ci_log      <- forest_res$pooled_ci
    out$heterogeneity_main <- forest_res$heterogeneity
  }

  class(out) <- c("IPD_two_stage", "list")
  out
}
