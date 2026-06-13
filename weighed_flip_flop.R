############################################################
## Masked flip-flop simulation for nu = d = 2
## Row/column missingness example
############################################################

rm(list = ls())

global_start <- Sys.time()

############################################################
## User controls
############################################################

master_seed <- 20260613

## Monte Carlo sample sizes and repetitions.
## Increase n_mc for smoother graphs.
n_grid <- c(25, 50, 100, 250, 500, 1000, 2500, 5000)
n_mc   <- 200

## One large run for iteration-level diagnostics.
n_demo <- 50000

maxit <- 1000
tol   <- 1e-12

graphs_dir <- file.path(getwd(), "graphs")
if (!dir.exists(graphs_dir)) {
  dir.create(graphs_dir, recursive = TRUE)
}

############################################################
## True matrix-normal scale matrices
##
## X ~ MN_{2 x 2}(0, Psi_true, Sigma_true)
## vec(X) has covariance Sigma_true %x% Psi_true.
############################################################

Psi_true <- matrix(c(1.50, 0.60,
                     0.60, 2.00), nrow = 2, byrow = TRUE)

Sigma_true <- matrix(c(2.20, -0.70,
                       -0.70,  1.40), nrow = 2, byrow = TRUE)

############################################################
## Masks
############################################################

I2 <- diag(2)

masks <- list(
  row_mask = list(
    P = diag(c(0, 1)),
    Q = I2,
    label = "Row mask: P = diag(0,1), Q = I2"
  ),
  column_mask = list(
    P = I2,
    Q = diag(c(0, 1)),
    label = "Column mask: P = I2, Q = diag(0,1)"
  )
)

############################################################
## Utility functions
############################################################

tr_mat <- function(A) {
  sum(diag(A))
}

sym_mat <- function(A) {
  0.5 * (A + t(A))
}

frob_norm <- function(A) {
  sqrt(sum(A * A))
}

mp_inv <- function(A, tol = sqrt(.Machine$double.eps)) {
  A <- sym_mat(A)
  ee <- eigen(A, symmetric = TRUE)
  cutoff <- tol * max(1, max(abs(ee$values)))
  inv_vals <- ifelse(ee$values > cutoff, 1 / ee$values, 0)
  ee$vectors %*% diag(inv_vals, nrow = length(inv_vals)) %*% t(ee$vectors)
}

make_targets <- function(Psi_true, Sigma_true, P, Q) {
  ## Scale convention:
  ## tr(Q Sigma Q) = tr(Q).
  ##
  ## Because (Psi, Sigma) are identifiable only up to inverse scaling,
  ## the true target under this convention is
  ## Sigma_target = Sigma_true / alpha,
  ## Psi_target   = Psi_true * alpha,
  ## where alpha = tr(Q Sigma_true Q) / tr(Q).
  alpha <- tr_mat(Q %*% Sigma_true %*% Q) / tr_mat(Q)
  
  Sigma_target_obs <- Q %*% (Sigma_true / alpha) %*% Q
  Psi_target_obs   <- P %*% (Psi_true * alpha) %*% P
  
  list(
    alpha = alpha,
    Sigma_obs = sym_mat(Sigma_target_obs),
    Psi_obs = sym_mat(Psi_target_obs),
    Kron_obs = kronecker(sym_mat(Sigma_target_obs), sym_mat(Psi_target_obs))
  )
}

calc_errors <- function(Psi_hat, Sigma_hat, targets, P, Q) {
  Sigma_hat_obs <- Q %*% Sigma_hat %*% Q
  Psi_hat_obs   <- P %*% Psi_hat %*% P
  
  Kron_hat_obs <- kronecker(Sigma_hat_obs, Psi_hat_obs)
  
  sigma_den <- max(frob_norm(targets$Sigma_obs), .Machine$double.eps)
  psi_den   <- max(frob_norm(targets$Psi_obs), .Machine$double.eps)
  kron_den  <- max(frob_norm(targets$Kron_obs), .Machine$double.eps)
  
  c(
    sigma_relerr = frob_norm(Sigma_hat_obs - targets$Sigma_obs) / sigma_den,
    psi_relerr   = frob_norm(Psi_hat_obs - targets$Psi_obs) / psi_den,
    kron_relerr  = frob_norm(Kron_hat_obs - targets$Kron_obs) / kron_den
  )
}

############################################################
## Matrix-normal simulation with masking
############################################################

simulate_masked_array <- function(n, Psi, Sigma, P, Q) {
  nu <- nrow(Psi)
  d  <- nrow(Sigma)
  
  L_psi   <- t(chol(Psi))
  L_sigma <- t(chol(Sigma))
  
  Y_array <- array(0, dim = c(nu, d, n))
  
  for (k in seq_len(n)) {
    Z <- matrix(rnorm(nu * d), nrow = nu, ncol = d)
    X <- L_psi %*% Z %*% t(L_sigma)
    Y_array[, , k] <- P %*% X %*% Q
  }
  
  Y_array
}

############################################################
## Masked flip-flop algorithm
############################################################

masked_flipflop <- function(Y_array,
                            P,
                            Q,
                            Psi_init = diag(dim(Y_array)[1]),
                            Sigma_init = diag(dim(Y_array)[2]),
                            maxit = 1000,
                            tol = 1e-12,
                            targets = NULL,
                            store_history = TRUE) {
  nu <- dim(Y_array)[1]
  d  <- dim(Y_array)[2]
  n  <- dim(Y_array)[3]
  
  trP <- tr_mat(P)
  trQ <- tr_mat(Q)
  
  Psi   <- Psi_init
  Sigma <- Sigma_init
  
  history <- NULL
  
  if (store_history && !is.null(targets)) {
    err0 <- calc_errors(Psi, Sigma, targets, P, Q)
    history <- data.frame(
      iter = 0,
      delta = NA_real_,
      sigma_relerr = err0["sigma_relerr"],
      psi_relerr = err0["psi_relerr"],
      kron_relerr = err0["kron_relerr"]
    )
  }
  
  converged <- FALSE
  
  for (it in seq_len(maxit)) {
    ########################################################
    ## Sigma update:
    ## Sigma_raw = 1 / (n tr(P)) sum_k Y_k' (P Psi P)^+ Y_k
    ########################################################
    
    A_plus <- mp_inv(P %*% Psi %*% P)
    
    Sigma_raw <- matrix(0, nrow = d, ncol = d)
    
    for (k in seq_len(n)) {
      Yk <- Y_array[, , k]
      Sigma_raw <- Sigma_raw + t(Yk) %*% A_plus %*% Yk
    }
    
    Sigma_raw <- sym_mat(Sigma_raw / (n * trP))
    Sigma_raw <- Q %*% Sigma_raw %*% Q
    
    ########################################################
    ## Psi update:
    ## Psi_raw = 1 / (n tr(Q)) sum_k Y_k (Q Sigma_raw Q)^+ Y_k'
    ########################################################
    
    C_plus <- mp_inv(Q %*% Sigma_raw %*% Q)
    
    Psi_raw <- matrix(0, nrow = nu, ncol = nu)
    
    for (k in seq_len(n)) {
      Yk <- Y_array[, , k]
      Psi_raw <- Psi_raw + Yk %*% C_plus %*% t(Yk)
    }
    
    Psi_raw <- sym_mat(Psi_raw / (n * trQ))
    Psi_raw <- P %*% Psi_raw %*% P
    
    ########################################################
    ## Identifiability rescaling:
    ## tr(Q Sigma Q) = tr(Q), inverse scaling on P Psi P.
    ########################################################
    
    alpha_hat <- tr_mat(Q %*% Sigma_raw %*% Q) / trQ
    
    if (!is.finite(alpha_hat) || alpha_hat <= 0) {
      stop("Non-positive or non-finite scaling factor encountered.")
    }
    
    Sigma_next <- sym_mat(Sigma_raw / alpha_hat)
    Psi_next   <- sym_mat(Psi_raw * alpha_hat)
    
    Sigma_next <- Q %*% Sigma_next %*% Q
    Psi_next   <- P %*% Psi_next %*% P
    
    ########################################################
    ## Convergence check on identifiable blocks only
    ########################################################
    
    delta_sigma <- frob_norm(Q %*% (Sigma_next - Sigma) %*% Q) /
      max(1, frob_norm(Q %*% Sigma %*% Q))
    
    delta_psi <- frob_norm(P %*% (Psi_next - Psi) %*% P) /
      max(1, frob_norm(P %*% Psi %*% P))
    
    delta <- max(delta_sigma, delta_psi)
    
    Sigma <- Sigma_next
    Psi   <- Psi_next
    
    if (store_history && !is.null(targets)) {
      err <- calc_errors(Psi, Sigma, targets, P, Q)
      history <- rbind(
        history,
        data.frame(
          iter = it,
          delta = delta,
          sigma_relerr = err["sigma_relerr"],
          psi_relerr = err["psi_relerr"],
          kron_relerr = err["kron_relerr"]
        )
      )
    }
    
    if (delta < tol) {
      converged <- TRUE
      break
    }
  }
  
  list(
    Psi_hat = Psi,
    Sigma_hat = Sigma,
    iterations = it,
    converged = converged,
    history = history
  )
}

############################################################
## Plotting functions
############################################################

save_iteration_plot <- function(history, main, file_name) {
  out_file <- file.path(graphs_dir, file_name)
  
  ymat <- cbind(
    pmax(history$sigma_relerr, .Machine$double.eps),
    pmax(history$psi_relerr, .Machine$double.eps),
    pmax(history$kron_relerr, .Machine$double.eps)
  )
  
  png(out_file, width = 1000, height = 700)
  matplot(
    history$iter,
    ymat,
    type = "b",
    pch = 1:3,
    lty = 1:3,
    log = "y",
    xlab = "Flip-flop iteration",
    ylab = "Relative error to true identifiable target",
    main = main
  )
  legend(
    "topright",
    legend = c("Sigma block", "Psi block", "Observed Kronecker covariance"),
    pch = 1:3,
    lty = 1:3,
    bty = "n"
  )
  dev.off()
  
  invisible(out_file)
}

summarise_metric <- function(dat, metric) {
  split_dat <- split(dat, list(dat$case_name, dat$n), drop = TRUE)
  
  out <- do.call(
    rbind,
    lapply(split_dat, function(dd) {
      x <- dd[[metric]]
      data.frame(
        case_name = dd$case_name[1],
        n = dd$n[1],
        metric = metric,
        mean = mean(x),
        median = median(x),
        q10 = unname(quantile(x, 0.10)),
        q90 = unname(quantile(x, 0.90)),
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(out) <- NULL
  out
}

save_mc_metric_plot <- function(summary_df, metric, ylab, file_name) {
  dd <- summary_df[summary_df$metric == metric, ]
  
  out_file <- file.path(graphs_dir, file_name)
  
  png(out_file, width = 1000, height = 700)
  
  plot(
    NA,
    xlim = range(dd$n),
    ylim = range(pmax(dd$median, .Machine$double.eps)),
    log = "xy",
    xlab = "Sample size n",
    ylab = ylab,
    main = paste("Statistical convergence:", metric)
  )
  
  case_names <- unique(dd$case_name)
  
  for (j in seq_along(case_names)) {
    tmp <- dd[dd$case_name == case_names[j], ]
    tmp <- tmp[order(tmp$n), ]
    
    lines(
      tmp$n,
      pmax(tmp$median, .Machine$double.eps),
      type = "b",
      pch = j,
      lty = j
    )
  }
  
  legend(
    "topright",
    legend = case_names,
    pch = seq_along(case_names),
    lty = seq_along(case_names),
    bty = "n"
  )
  
  dev.off()
  
  invisible(out_file)
}

save_iteration_count_plot <- function(dat, file_name) {
  split_dat <- split(dat, list(dat$case_name, dat$n), drop = TRUE)
  
  iter_summary <- do.call(
    rbind,
    lapply(split_dat, function(dd) {
      data.frame(
        case_name = dd$case_name[1],
        n = dd$n[1],
        mean_iterations = mean(dd$iterations),
        median_iterations = median(dd$iterations),
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(iter_summary) <- NULL
  
  out_file <- file.path(graphs_dir, file_name)
  
  png(out_file, width = 1000, height = 700)
  
  plot(
    NA,
    xlim = range(iter_summary$n),
    ylim = range(iter_summary$mean_iterations),
    log = "x",
    xlab = "Sample size n",
    ylab = "Mean number of flip-flop iterations",
    main = "Iteration counts"
  )
  
  case_names <- unique(iter_summary$case_name)
  
  for (j in seq_along(case_names)) {
    tmp <- iter_summary[iter_summary$case_name == case_names[j], ]
    tmp <- tmp[order(tmp$n), ]
    
    lines(
      tmp$n,
      tmp$mean_iterations,
      type = "b",
      pch = j,
      lty = j
    )
  }
  
  legend(
    "topright",
    legend = case_names,
    pch = seq_along(case_names),
    lty = seq_along(case_names),
    bty = "n"
  )
  
  dev.off()
  
  invisible(out_file)
}

############################################################
## Print true identifiable targets
############################################################

cat("\nTrue scale matrices used to generate data:\n")
cat("\nPsi_true:\n")
print(Psi_true)
cat("\nSigma_true:\n")
print(Sigma_true)

for (case_name in names(masks)) {
  P <- masks[[case_name]]$P
  Q <- masks[[case_name]]$Q
  targets <- make_targets(Psi_true, Sigma_true, P, Q)
  
  cat("\n============================================================\n")
  cat(masks[[case_name]]$label, "\n")
  cat("Scale alpha used for true target:", targets$alpha, "\n")
  cat("\nIdentifiable Sigma target, Q Sigma Q, under convention:\n")
  print(round(targets$Sigma_obs, 6))
  cat("\nIdentifiable Psi target, P Psi P, under convention:\n")
  print(round(targets$Psi_obs, 6))
}

############################################################
## Large single-dataset diagnostics
############################################################

set.seed(master_seed)

for (case_name in names(masks)) {
  P <- masks[[case_name]]$P
  Q <- masks[[case_name]]$Q
  
  targets <- make_targets(Psi_true, Sigma_true, P, Q)
  
  Y_demo <- simulate_masked_array(
    n = n_demo,
    Psi = Psi_true,
    Sigma = Sigma_true,
    P = P,
    Q = Q
  )
  
  fit_demo <- masked_flipflop(
    Y_array = Y_demo,
    P = P,
    Q = Q,
    Psi_init = diag(2),
    Sigma_init = diag(2),
    maxit = maxit,
    tol = tol,
    targets = targets,
    store_history = TRUE
  )
  
  cat("\n============================================================\n")
  cat("Large-demo fit:", case_name, "\n")
  cat("Converged:", fit_demo$converged, "\n")
  cat("Iterations:", fit_demo$iterations, "\n")
  
  cat("\nEstimated identifiable Sigma block:\n")
  print(round(Q %*% fit_demo$Sigma_hat %*% Q, 6))
  
  cat("\nTrue identifiable Sigma block:\n")
  print(round(targets$Sigma_obs, 6))
  
  cat("\nEstimated identifiable Psi block:\n")
  print(round(P %*% fit_demo$Psi_hat %*% P, 6))
  
  cat("\nTrue identifiable Psi block:\n")
  print(round(targets$Psi_obs, 6))
  
  final_errors <- calc_errors(
    Psi_hat = fit_demo$Psi_hat,
    Sigma_hat = fit_demo$Sigma_hat,
    targets = targets,
    P = P,
    Q = Q
  )
  
  cat("\nFinal relative errors:\n")
  print(final_errors)
  
  save_iteration_plot(
    history = fit_demo$history,
    main = paste("Iteration diagnostics:", case_name),
    file_name = paste0(case_name, "_iteration_errors.png")
  )
}

############################################################
## Parallel Monte Carlo convergence study
##
## Windows note:
## makeCluster() creates PSOCK workers, so every function and object
## needed by run_one_task() is explicitly exported below.
############################################################

set.seed(master_seed)

tasks <- expand.grid(
  case_name = names(masks),
  n = n_grid,
  rep = seq_len(n_mc),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

tasks$seed <- sample.int(.Machine$integer.max, size = nrow(tasks))

task_list <- lapply(seq_len(nrow(tasks)), function(i) {
  as.list(tasks[i, ])
})

run_one_task <- function(task) {
  set.seed(task$seed)
  
  case_name <- task$case_name
  n <- task$n
  
  P <- masks[[case_name]]$P
  Q <- masks[[case_name]]$Q
  
  targets <- make_targets(Psi_true, Sigma_true, P, Q)
  
  Y <- simulate_masked_array(
    n = n,
    Psi = Psi_true,
    Sigma = Sigma_true,
    P = P,
    Q = Q
  )
  
  fit <- masked_flipflop(
    Y_array = Y,
    P = P,
    Q = Q,
    Psi_init = diag(2),
    Sigma_init = diag(2),
    maxit = maxit,
    tol = tol,
    targets = targets,
    store_history = FALSE
  )
  
  errs <- calc_errors(
    Psi_hat = fit$Psi_hat,
    Sigma_hat = fit$Sigma_hat,
    targets = targets,
    P = P,
    Q = Q
  )
  
  data.frame(
    case_name = case_name,
    n = n,
    rep = task$rep,
    converged = fit$converged,
    iterations = fit$iterations,
    sigma_relerr = errs["sigma_relerr"],
    psi_relerr = errs["psi_relerr"],
    kron_relerr = errs["kron_relerr"],
    stringsAsFactors = FALSE
  )
}

n_cores_available <- parallel::detectCores(logical = TRUE)
if (is.na(n_cores_available)) {
  n_cores_available <- 2
}

n_workers <- max(1, n_cores_available - 1)

cat("\n============================================================\n")
cat("Starting parallel Monte Carlo with", n_workers, "workers.\n")
cat("Total tasks:", length(task_list), "\n")

cl <- parallel::makeCluster(n_workers)

export_names <- c(
  "tr_mat",
  "sym_mat",
  "frob_norm",
  "mp_inv",
  "make_targets",
  "calc_errors",
  "simulate_masked_array",
  "masked_flipflop",
  "run_one_task",
  "Psi_true",
  "Sigma_true",
  "masks",
  "maxit",
  "tol"
)

parallel::clusterExport(
  cl = cl,
  varlist = export_names,
  envir = environment()
)

mc_list <- tryCatch(
  {
    parallel::parLapplyLB(cl, task_list, run_one_task)
  },
  finally = {
    parallel::stopCluster(cl)
  }
)

mc_results <- do.call(rbind, mc_list)
rownames(mc_results) <- NULL

############################################################
## Save Monte Carlo results and graphs
############################################################

mc_results_file <- file.path(graphs_dir, "masked_flipflop_mc_results.csv")
write.csv(mc_results, mc_results_file, row.names = FALSE)

summary_sigma <- summarise_metric(mc_results, "sigma_relerr")
summary_psi   <- summarise_metric(mc_results, "psi_relerr")
summary_kron  <- summarise_metric(mc_results, "kron_relerr")

mc_summary <- rbind(summary_sigma, summary_psi, summary_kron)

mc_summary_file <- file.path(graphs_dir, "masked_flipflop_mc_summary.csv")
write.csv(mc_summary, mc_summary_file, row.names = FALSE)

save_mc_metric_plot(
  summary_df = mc_summary,
  metric = "sigma_relerr",
  ylab = "Median relative error for identifiable Sigma block",
  file_name = "statistical_convergence_sigma.png"
)

save_mc_metric_plot(
  summary_df = mc_summary,
  metric = "psi_relerr",
  ylab = "Median relative error for identifiable Psi block",
  file_name = "statistical_convergence_psi.png"
)

save_mc_metric_plot(
  summary_df = mc_summary,
  metric = "kron_relerr",
  ylab = "Median relative error for observed Kronecker covariance",
  file_name = "statistical_convergence_observed_kron.png"
)

save_iteration_count_plot(
  dat = mc_results,
  file_name = "iteration_counts.png"
)

############################################################
## Console summary
############################################################

cat("\n============================================================\n")
cat("Monte Carlo convergence summary:\n")
cat("Results saved to:\n")
cat(mc_results_file, "\n")
cat(mc_summary_file, "\n")

cat("\nGraphs saved in:\n")
cat(graphs_dir, "\n\n")

cat("Convergence rate by case:\n")
print(
  aggregate(
    converged ~ case_name + n,
    data = mc_results,
    FUN = mean
  )
)

cat("\nMedian final relative errors by case and n:\n")
print(
  mc_summary[mc_summary$metric %in% c("sigma_relerr", "psi_relerr", "kron_relerr"),
             c("case_name", "n", "metric", "median", "q10", "q90")]
)

global_end <- Sys.time()
elapsed_minutes <- as.numeric(difftime(global_end, global_start, units = "mins"))

cat("\n============================================================\n")
cat(sprintf("Total elapsed time: %.3f minutes\n", elapsed_minutes))
cat("Done.\n")