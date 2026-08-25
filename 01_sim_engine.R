
required_pkgs <- c("BayesFactor", "doParallel", "foreach", "dplyr", "ggplot2", "scales")

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Eksik paketler yukleniyor: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(BayesFactor)
  library(doParallel)
  library(foreach)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ---- 1. Design ---------------------------------------------------------
for (d in c("results", "results/checkpoints", "figures", "tables")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- 2. Simulation parameters ---------------------------------------------
SIM_PARAMS <- list(
  true_delta      = c(0, 0.2, 0.5),       # 0.8 cikarildi (tavan etkisi)
  prior_r0        = c(0.3, 0.5, 0.707),
  drift_magnitude = c(0, 0.1, 0.3, 0.6),  # asagi yonlu: r azalacak
  drift_speed     = c(10, 50),             # Inf cikarildi (baseline ile ayni)
  drift_type      = c("gradual", "sudden"),
  n_min           = 20,
  n_max           = 500,
  bf_upper        = 10,
  bf_lower        = 0.1,
  K               = 2000,
  seed_base       = 20240602,
  r_min           = 0.1    # prior olcegi alt siniri
)

# ---- 3.  -----------------------------------------------------
build_grid <- function(p) {
  grid <- expand.grid(
    true_delta      = p$true_delta,
    prior_r0        = p$prior_r0,
    drift_magnitude = p$drift_magnitude,
    drift_speed     = p$drift_speed,
    drift_type      = p$drift_type,
    stringsAsFactors = FALSE
  )

  # drift = 0 icin speed ve type tekillesir
  grid$drift_type[grid$drift_magnitude == 0]  <- "none"
  grid$drift_speed[grid$drift_magnitude == 0] <- Inf
  grid <- unique(grid)

  # GECERSIZ KOMBINASYONLARI CIKAR:
  # r0 - drift_magnitude >= r_min olmali
  # Gradual drift icin: en fazla drift_magnitude kadar daralir
  # Sudden drift icin: ayni kural gecerli
  grid <- grid[
    grid$drift_magnitude == 0 |
    (grid$prior_r0 - grid$drift_magnitude >= p$r_min),
  ]

  grid$condition_id <- seq_len(nrow(grid))
  grid
}

GRID <- build_grid(SIM_PARAMS)
message(sprintf("Toplam kosul sayisi: %d", nrow(GRID)))
message("Izgara ozeti:")
print(table(GRID$prior_r0, GRID$drift_magnitude))

# ---- 4.  -----------------------------------

compute_drifted_prior_down <- function(n_current, r0, magnitude, speed, type,
                                       r_min = 0.1) {
  if (magnitude == 0 || type == "none" || is.infinite(speed)) return(r0)

  N_MAX <- 500L

  if (type == "gradual") {
    max_steps  <- floor(N_MAX / speed)
    steps_done <- floor(n_current / speed)
    drift      <- min(steps_done, max_steps) * (magnitude / max(max_steps, 1))
  } else {
    # sudden: speed. gozlemden sonra tum drift birdenbire iner
    drift <- ifelse(n_current >= speed, magnitude, 0)
  }

  new_r <- r0 - drift                  # ASAGI: cikarma
  max(r_min, min(new_r, 2.0))          # [r_min, 2.0] araliginda tut
}

# ---- 5.  -------------------------------------------------
simulate_condition <- function(true_delta, prior_r0,
                               drift_magnitude, drift_speed, drift_type,
                               K, n_min, n_max, bf_upper, bf_lower,
                               r_min, seed) {
  set.seed(seed)

  out_decision <- character(K)
  out_n_stop   <- integer(K)
  out_final_bf <- numeric(K)
  out_final_r  <- numeric(K)

  for (k in seq_len(K)) {
    decision <- "inconclusive"
    n_stop   <- n_max
    final_bf <- NA_real_
    final_r  <- prior_r0

    x <- rnorm(n_max, mean = true_delta, sd = 1)

    for (n in n_min:n_max) {
      r_now <- compute_drifted_prior_down(
        n_current = n,
        r0        = prior_r0,
        magnitude = drift_magnitude,
        speed     = drift_speed,
        type      = drift_type,
        r_min     = r_min
      )

      bf_val <- tryCatch({
        bf_obj <- ttestBF(x = x[1:n], rscale = r_now)
        as.numeric(exp(bf_obj@bayesFactor$bf[1]))
      }, error = function(e) NA_real_)

      if (is.na(bf_val) || !is.finite(bf_val)) next

      if (bf_val >= bf_upper) {
        decision <- "H1"; n_stop <- n; final_bf <- bf_val; final_r <- r_now
        break
      }
      if (bf_val <= bf_lower) {
        decision <- "H0"; n_stop <- n; final_bf <- bf_val; final_r <- r_now
        break
      }
    }

    if (decision == "inconclusive") {
      final_r  <- compute_drifted_prior_down(
        n_max, prior_r0, drift_magnitude, drift_speed, drift_type, r_min
      )
      final_bf <- tryCatch({
        bf_obj <- ttestBF(x = x, rscale = final_r)
        as.numeric(exp(bf_obj@bayesFactor$bf[1]))
      }, error = function(e) NA_real_)
    }

    out_decision[k] <- decision
    out_n_stop[k]   <- n_stop
    out_final_bf[k] <- final_bf
    out_final_r[k]  <- final_r
  }

  data.frame(
    decision = out_decision,
    n_stop   = out_n_stop,
    final_bf = out_final_bf,
    final_r  = out_final_r,
    stringsAsFactors = FALSE
  )
}

# ---- 6.  (parallel + checkpoint) -----------------------
run_all_simulations <- function(grid, params, n_cores = NULL) {

  if (is.null(n_cores)) {
    n_cores <- max(1L, parallel::detectCores() - 1L)
  }
  n_cores <- min(n_cores, nrow(grid))
  message(sprintf("Paralel cekirdek: %d", n_cores))

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)

  done_files <- list.files("results/checkpoints",
                            pattern = "^checkpoint_\\d+\\.rds$",
                            full.names = FALSE)
  done_ids <- as.integer(gsub("checkpoint_(\\d+)\\.rds", "\\1", done_files))
  done_ids <- done_ids[!is.na(done_ids)]
  todo_ids <- setdiff(grid$condition_id, done_ids)

  if (length(todo_ids) == 0) {
    message("Tum kosullar tamamlanmis. Birlestiriliyor...")
    return(merge_checkpoints())
  }

  message(sprintf("%d kosuldan %d tanesi kaldi.", nrow(grid), length(todo_ids)))
  todo_grid <- grid[grid$condition_id %in% todo_ids, ]

  fns_to_export <- c("simulate_condition", "compute_drifted_prior_down")

  foreach(
    i              = seq_len(nrow(todo_grid)),
    .packages      = "BayesFactor",
    .export        = fns_to_export,
    .errorhandling = "pass",
    .verbose       = FALSE
  ) %dopar% {

    row  <- todo_grid[i, ]
    cid  <- row$condition_id
    seed <- params$seed_base + cid

    sim_df <- simulate_condition(
      true_delta      = row$true_delta,
      prior_r0        = row$prior_r0,
      drift_magnitude = row$drift_magnitude,
      drift_speed     = row$drift_speed,
      drift_type      = row$drift_type,
      K               = params$K,
      n_min           = params$n_min,
      n_max           = params$n_max,
      bf_upper        = params$bf_upper,
      bf_lower        = params$bf_lower,
      r_min           = params$r_min,
      seed            = seed
    )

    sim_df$condition_id    <- cid
    sim_df$true_delta      <- row$true_delta
    sim_df$prior_r0        <- row$prior_r0
    sim_df$drift_magnitude <- row$drift_magnitude
    sim_df$drift_speed     <- row$drift_speed
    sim_df$drift_type      <- row$drift_type

    cp_path <- file.path("results", "checkpoints",
                         sprintf("checkpoint_%04d.rds", cid))
    saveRDS(sim_df, cp_path)
    NULL
  }

  message("Tum kosullar tamamlandi. Checkpointler birlestiriliyor...")
  merge_checkpoints()
}

# ---- 7. ---------------------------------------------
merge_checkpoints <- function() {
  cp_files <- list.files("results/checkpoints",
                          pattern = "^checkpoint_\\d+\\.rds$",
                          full.names = TRUE)
  if (length(cp_files) == 0) stop("Hic checkpoint bulunamadi.")
  combined <- do.call(rbind, lapply(cp_files, readRDS))
  saveRDS(combined, "results/all_simulations.rds")
  message(sprintf("Toplam %d satir -> results/all_simulations.rds", nrow(combined)))
  combined
}

# ---- 8. -----------------------------------------------
compute_error_rates <- function(sim_data) {
  sim_data %>%
    group_by(condition_id, true_delta, prior_r0,
             drift_magnitude, drift_speed, drift_type) %>%
    summarise(
      K             = n(),
      type_I_error  = ifelse(first(true_delta) == 0,
                             mean(decision == "H1"), NA_real_),
      type_II_error = ifelse(first(true_delta) > 0,
                             mean(decision == "H0"), NA_real_),
      power         = ifelse(first(true_delta) > 0,
                             mean(decision == "H1"), NA_real_),
      inconclusive  = mean(decision == "inconclusive"),
      median_n      = median(n_stop),
      mean_n        = mean(n_stop),
      .groups       = "drop"
    )
}

message("01_sim_engine.R yuklendi.")
message(sprintf("Izgara: %d kosul, K = %d", nrow(GRID), SIM_PARAMS$K))
message("Calistirmak icin:")
message("  sim_data    <- run_all_simulations(GRID, SIM_PARAMS)")
message("  error_rates <- compute_error_rates(sim_data)")
message("  saveRDS(error_rates, 'results/error_rates.rds')")
