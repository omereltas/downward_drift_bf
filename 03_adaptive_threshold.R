
if (!exists("compute_drifted_prior_down")) {
  message("01_sim_engine.R yukleniyor...")
  source("01_sim_engine.R")
}

suppressPackageStartupMessages({
  library(BayesFactor)
  library(doParallel)
  library(foreach)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

for (d in c("results", "results/calib_checkpoints", "figures", "tables")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- 1. Type I error -----------------------------------------------
estimate_typeI_down <- function(bf_upper, prior_r0, drift_magnitude,
                                drift_speed, drift_type,
                                n_min = 20, n_max = 500,
                                r_min = 0.1, K = 2000, seed = 42) {
  set.seed(seed)
  h1_decisions <- 0L

  for (k in seq_len(K)) {
    x <- rnorm(n_max, mean = 0, sd = 1)

    for (n in n_min:n_max) {
      r_now <- compute_drifted_prior_down(n, prior_r0, drift_magnitude,
                                          drift_speed, drift_type, r_min)
      bf_val <- tryCatch({
        as.numeric(exp(
          ttestBF(x = x[1:n], rscale = r_now)@bayesFactor$bf[1]
        ))
      }, error = function(e) NA_real_)

      if (is.na(bf_val) || !is.finite(bf_val)) next
      if (bf_val >= bf_upper)   { h1_decisions <- h1_decisions + 1L; break }
      if (bf_val <= 1/bf_upper) break
    }
  }
  h1_decisions / K
}

# ---- 2.  ------------------------------------
find_calibrated_threshold_down <- function(prior_r0, drift_magnitude,
                                           drift_speed, drift_type,
                                           target_alpha   = 0.05,
                                           tol            = 0.004,
                                           early_stop_tol = 0.01,
                                           early_stop_n   = 3,
                                           K              = 2000,
                                           lo             = 3,
                                           hi             = 300,
                                           max_iter       = 30,
                                           r_min          = 0.1,
                                           seed_base      = 99999) {
  achieved    <- NA_real_
  mid         <- NA_real_
  mid_history <- numeric(max_iter)
  stop_reason <- "max_iter"

  for (iter in seq_len(max_iter)) {
    mid <- (lo + hi) / 2
    mid_history[iter] <- mid

    t1 <- estimate_typeI_down(
      bf_upper        = mid,
      prior_r0        = prior_r0,
      drift_magnitude = drift_magnitude,
      drift_speed     = drift_speed,
      drift_type      = drift_type,
      r_min           = r_min,
      K               = K,
      seed            = seed_base + iter
    )
    achieved <- t1

    if (abs(t1 - target_alpha) <= tol) {
      stop_reason <- "converged"; break
    }

    if (iter >= early_stop_n) {
      recent <- mid_history[(iter - early_stop_n + 1):iter]
      if (max(recent) - min(recent) < early_stop_tol) {
        stop_reason <- "early_stop"; break
      }
    }

    if (t1 > target_alpha) lo <- mid else hi <- mid
  }

  list(bf_threshold   = mid,
       achieved_typeI = achieved,
       n_iter         = iter,
       stop_reason    = stop_reason)
}

# ---- 3. Calibration grid ---------------------------------------------------
calib_grid <- expand.grid(
  prior_r0        = c(0.3, 0.5, 0.707),
  drift_magnitude = c(0, 0.1, 0.3, 0.6),
  drift_speed     = c(10, 50),
  drift_type      = c("gradual", "sudden"),
  stringsAsFactors = FALSE
)

calib_grid$drift_type[calib_grid$drift_magnitude == 0]  <- "none"
calib_grid$drift_speed[calib_grid$drift_magnitude == 0] <- Inf
calib_grid <- unique(calib_grid)

calib_grid <- calib_grid[
  calib_grid$drift_magnitude == 0 |
  (calib_grid$prior_r0 - calib_grid$drift_magnitude >= 0.1),
]

calib_grid$row_id <- seq_len(nrow(calib_grid))
message(sprintf("Kalibrasyon kosul sayisi: %d", nrow(calib_grid)))

# ---- 4. Tamamlanan satirlari bul (checkpoint'lerden) ------------------------
# Her worker kendi checkpoint'ini results/calib_checkpoints/calib_XXXX.rds
# olarak yazacak. Boylece crash durumunda sadece tamamlanmayanlar yeniden calisir.

find_done_ids <- function() {
  cp_files <- list.files("results/calib_checkpoints",
                          pattern = "^calib_\\d+\\.rds$",
                          full.names = FALSE)
  as.integer(gsub("calib_(\\d+)\\.rds", "\\1", cp_files))
}

done_ids <- find_done_ids()
todo_ids <- setdiff(calib_grid$row_id, done_ids)

if (length(todo_ids) == 0) {
  message("Tum kalibrasyon kosullari tamamlanmis.")
} else {
  message(sprintf("%d kosuldan %d tanesi kaldi.", nrow(calib_grid), length(todo_ids)))

  n_cores <- max(1L, parallel::detectCores() - 1L)
  n_cores <- min(n_cores, length(todo_ids))
  message(sprintf("Paralel cekirdek: %d", n_cores))

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)

  fns_to_export <- c("estimate_typeI_down",
                     "find_calibrated_threshold_down",
                     "compute_drifted_prior_down")

  foreach(
    i              = todo_ids,
    .packages      = "BayesFactor",
    .export        = fns_to_export,
    .errorhandling = "pass",
    .verbose       = FALSE
  ) %dopar% {

    row <- calib_grid[calib_grid$row_id == i, ]

    res <- find_calibrated_threshold_down(
      prior_r0        = row$prior_r0,
      drift_magnitude = row$drift_magnitude,
      drift_speed     = row$drift_speed,
      drift_type      = row$drift_type,
      K               = 2000,
      seed_base       = 99999 + i * 100
    )

    result_df <- data.frame(
      row_id          = row$row_id,
      prior_r0        = row$prior_r0,
      drift_magnitude = row$drift_magnitude,
      drift_speed     = as.numeric(row$drift_speed),
      drift_type      = row$drift_type,
      bf_threshold    = res$bf_threshold,
      achieved_typeI  = res$achieved_typeI,
      n_iter          = res$n_iter,
      stop_reason     = res$stop_reason,
      stringsAsFactors = FALSE
    )

    # Her worker kendi checkpoint'ini yazar
    cp_path <- file.path("results", "calib_checkpoints",
                         sprintf("calib_%04d.rds", i))
    saveRDS(result_df, cp_path)
    NULL
  }

  message("Paralel kalibrasyon tamamlandi.")
}

# ---- 5.  -------------------------------------------
cp_files <- list.files("results/calib_checkpoints",
                        pattern = "^calib_\\d+\\.rds$",
                        full.names = TRUE)

if (length(cp_files) == 0) stop("Hic kalibrasyon checkpoint'i bulunamadi.")

calib_df <- do.call(rbind, lapply(cp_files, readRDS))
calib_df  <- calib_df[order(calib_df$row_id), ]
calib_df$drift_speed <- as.numeric(calib_df$drift_speed)

# Durdurma nedeni ozeti
message("Durdurma nedeni ozeti:")
print(table(calib_df$stop_reason))
message(sprintf("Iterasyon: ort=%.1f, min=%d, max=%d",
                mean(calib_df$n_iter),
                min(calib_df$n_iter),
                max(calib_df$n_iter)))

# ---- 6. Baseline and CF  ----------------------------------------------
baseline <- calib_df %>%
  filter(drift_magnitude == 0) %>%
  select(prior_r0, baseline_bf = bf_threshold)

calib_df <- calib_df %>%
  left_join(baseline, by = "prior_r0") %>%
  mutate(correction_factor = bf_threshold / baseline_bf)

write.csv(calib_df, "tables/table3_calibrated_thresholds_downward.csv",
          row.names = FALSE)
message("Tablo 3 -> tables/table3_calibrated_thresholds_downward.csv")

# ---- 7.  -----------------------------------------------
fit_df <- calib_df %>%
  filter(drift_magnitude > 0, is.finite(drift_speed)) %>%
  mutate(
    log_correction = log(correction_factor),
    log_speed      = log(drift_speed),
    is_sudden      = as.integer(drift_type == "sudden"),
    r0_c           = prior_r0 - mean(prior_r0)
  )

if (nrow(fit_df) >= 4) {
  lm_fit <- lm(
    log_correction ~ drift_magnitude + log_speed + is_sudden +
      drift_magnitude:is_sudden + r0_c + drift_magnitude:r0_c,
    data = fit_df
  )
  cat("\n=== Duzeltme Formulu (OLS fit) ===\n")
  print(summary(lm_fit))
  cat("\n=== VIF (multikolinearlik kontrolu) ===\n")
  print(car::vif(lm_fit))
  saveRDS(lm_fit, "results/correction_model_downward.rds")
}

# ---- 8. Figure 3: Calibration line ----------------------------------------
plot_df <- calib_df %>%
  filter(drift_magnitude > 0, is.finite(drift_speed)) %>%
  mutate(speed_label = factor(as.character(drift_speed), levels = c("10","50")))

baseline_lines <- calib_df %>%
  filter(drift_magnitude == 0) %>%
  distinct(prior_r0, baseline_bf)

if (nrow(plot_df) > 0) {
  p4 <- ggplot(plot_df,
               aes(x        = drift_magnitude,
                   y        = bf_threshold,
                   color    = factor(prior_r0),
                   linetype = drift_type,
                   shape    = speed_label,
                   group    = interaction(prior_r0, drift_type, speed_label))) +
    geom_hline(data        = baseline_lines,
               aes(yintercept = baseline_bf, color = factor(prior_r0)),
               linetype = "dashed", linewidth = 0.5, alpha = 0.6,
               inherit.aes = FALSE) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 3) +
    facet_wrap(~ speed_label, labeller = label_both) +
    scale_color_manual(values = c("#2C7BB6", "#E66100", "#009E73"),
                       name   = "Prior r\u2080") +
    scale_linetype_manual(values = c("gradual" = "solid",
                                     "sudden"  = "dotted",
                                     "none"    = "dashed"),
                          name   = "Drift type") +
    scale_shape_manual(values = c("10" = 16, "50" = 17),
                       name   = "Drift speed") +
    scale_y_continuous(breaks = pretty_breaks(6)) +
    labs(
      title    = expression("Calibrated BF Threshold (BF"^"*"*
                            ") \u2014 Downward Drift"),
      subtitle = expression(
        "Minimum BF threshold to maintain Type I error at "*alpha*" = .05"),
      x        = expression(
        "Drift magnitude ("*Delta*italic(r)*", downward)"),
      y        = expression("Calibrated threshold (BF"^"*"*")"),
      caption  = "Dashed lines: no-drift baseline threshold per prior scale"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position  = "bottom",
          strip.background = element_rect(fill = "gray95"),
          panel.grid.minor = element_blank())

  ggsave("figures/fig3_calibration_downward.pdf", p4, width = 8, height = 5.5)
  ggsave("figures/fig3_calibration_downward.png", p4, width = 8, height = 5.5,
         dpi = 300)
  message("Figure 3 kaydedildi.")
}

message("=== 03_adaptive_threshold.R tamamlandi ===")


