## ----------------------------------------------
## Parallel weighted flip–flop (nu = d = 2) with masks & diagnostics
## Streams data; parallelizes per-replicate sums over cores-1 (PSOCK).
## ----------------------------------------------
set.seed(123)
suppressWarnings(suppressMessages({
  library(parallel)
}))

## ===== Helpers =====
symmetrize <- function(A) (A + t(A))/2
fro <- function(A) sqrt(sum(A*A))
mat_sqrt <- function(A) {
  ev <- eigen(symmetrize(A), symmetric = TRUE)
  Q <- ev$vectors; vals <- pmax(ev$values, 0)
  Q %*% diag(sqrt(vals), nrow=length(vals)) %*% t(Q)
}
mat_invsqrt <- function(A) {
  ev <- eigen(symmetrize(A), symmetric = TRUE)
  Q <- ev$vectors; vals <- pmax(ev$values, 1e-12)
  Q %*% diag(1/sqrt(vals), nrow=length(vals)) %*% t(Q)
}
psd_project <- function(A, eps = 1e-12) {
  ev <- eigen(symmetrize(A), symmetric = TRUE)
  vals <- pmax(ev$values, eps)
  ev$vectors %*% diag(vals, nrow=length(vals)) %*% t(ev$vectors)
}
tr <- function(A) sum(diag(A))

logdet_spd <- function(A) {
  ev <- eigen(symmetrize(A), symmetric = TRUE, only.values = TRUE)$values
  sum(log(pmax(ev, 1e-15)))
}

## ===== True parameters (2x2, SPD) =====
Psi_true   <- matrix(c(1, 0.5,
                       0.5, 1.5), 2, 2, byrow=TRUE)
Sigma_true <- matrix(c(1.2, 0.4,
                       0.4, 0.9), 2, 2, byrow=TRUE)
nu <- 2; d <- 2

## Precompute square-roots for generation
Psi_true_sqrt   <- mat_sqrt(Psi_true)
Sigma_true_sqrt <- mat_sqrt(Sigma_true)

## ===== Simulation settings =====
n <- 1e7  ## total number of iid matrices; increase as desired

## Partition n across workers deterministically
cores_use <- max(1, detectCores() - 1L)
cl <- makeCluster(cores_use)
on.exit(try(stopCluster(cl), silent=TRUE), add=TRUE)

block_sizes <- rep(n %/% cores_use, cores_use)
if ((n %% cores_use) > 0) block_sizes[seq_len(n %% cores_use)] <- block_sizes[seq_len(n %% cores_use)] + 1
block_seeds <- 4321L + seq_len(cores_use)  ## fixed per block -> identical dataset each pass

## Export static items to workers
clusterExport(cl, c("nu","d","Psi_true_sqrt","Sigma_true_sqrt",
                    "symmetrize","psd_project","mat_sqrt","mat_invsqrt"),
              envir=environment())

## ===== Masks (simple missingness) =====
use_row_mask    <- TRUE    # TRUE -> mask a row in Ψ-equation / Σ-equation weights
use_column_mask <- FALSE    # TRUE -> mask a column in Σ-equation / Ψ-equation weights
mask_row_index  <- 1        # which row to mask if use_row_mask=TRUE (1 or 2)
mask_col_index  <- 2        # which column to mask if use_column_mask=TRUE (1 or 2)

P <- diag(nu)  # row mask
Q <- diag(d)   # column mask
if (use_row_mask)    P[mask_row_index, mask_row_index] <- 0
if (use_column_mask) Q[mask_col_index, mask_col_index] <- 0

W <- P  # weight in Σ-equation (row-side)
U <- Q  # weight in Ψ-equation (column-side)
if (tr(W) <= 0 || tr(U) <= 0) stop("Mask trace is zero: cannot update with empty information.")

## ===== Identifiability: we compare to trace-normalized Σ =====
c_truth        <- d / tr(Sigma_true)
Sigma_true_sc  <- Sigma_true *  c_truth
Psi_true_sc    <- Psi_true   /  c_truth

obs_rows <- which(diag(P) == 1)
obs_cols <- which(diag(Q) == 1)
subblock <- function(A, idx) A[idx, idx, drop=FALSE]

## ===== Utilities: parallel block kernels (each uses its own seed) =====
## Each kernel loops over its assigned replicates, generates X = Psi_sqrt Z Sigma_sqrt,
## applies the masks, and accumulates the desired sum.

## Σ-update kernel: returns sum_k t(X_eff) %*% Y %*% X_eff
block_sum_S <- function(arg) {
  seed <- arg$seed; n_i <- arg$n_i
  set.seed(seed)
  S <- matrix(0, d, d)
  for (kk in seq_len(n_i)) {
    Z  <- matrix(rnorm(nu*d), nu, d)
    Xk <- Psi_true_sqrt %*% Z %*% Sigma_true_sqrt
    X_eff <- Xk %*% Q                          # honor column mask in Σ-update
    S <- S + t(X_eff) %*% Y %*% X_eff
  }
  S
}

## Ψ-update kernel: returns sum_k X_eff %*% Zmask %*% t(X_eff)
block_sum_P <- function(arg) {
  seed <- arg$seed; n_i <- arg$n_i
  set.seed(seed)
  Psum <- matrix(0, nu, nu)
  for (kk in seq_len(n_i)) {
    Z  <- matrix(rnorm(nu*d), nu, d)
    Xk <- Psi_true_sqrt %*% Z %*% Sigma_true_sqrt
    X_eff <- P %*% Xk                           # honor row mask in Ψ-update
    Psum <- Psum + X_eff %*% Zmask %*% t(X_eff)
  }
  Psum
}

## Residual Σ kernel at (Ψ*,Σ*): same as Σ kernel with Y_res
block_sum_S_res <- function(arg) {
  seed <- arg$seed; n_i <- arg$n_i
  set.seed(seed)
  S <- matrix(0, d, d)
  for (kk in seq_len(n_i)) {
    Z  <- matrix(rnorm(nu*d), nu, d)
    Xk <- Psi_true_sqrt %*% Z %*% Sigma_true_sqrt
    X_eff <- Xk %*% Q
    S <- S + t(X_eff) %*% Y_res %*% X_eff
  }
  S
}

## Residual Ψ kernel at (Ψ*,Σ*): same as Ψ kernel with Z_res
block_sum_P_res <- function(arg) {
  seed <- arg$seed; n_i <- arg$n_i
  set.seed(seed)
  Psum <- matrix(0, nu, nu)
  for (kk in seq_len(n_i)) {
    Z  <- matrix(rnorm(nu*d), nu, d)
    Xk <- Psi_true_sqrt %*% Z %*% Sigma_true_sqrt
    X_eff <- P %*% Xk
    Psum <- Psum + X_eff %*% Z_res %*% t(X_eff)
  }
  Psum
}

## Log-likelihood quadratic term kernel at (Ψ,Σ): sum_k tr(Σ^{-1} X^T Ψ^{-1} X)
block_sum_quad <- function(arg) {
  seed <- arg$seed; n_i <- arg$n_i
  set.seed(seed)
  acc <- 0
  for (kk in seq_len(n_i)) {
    Z  <- matrix(rnorm(nu*d), nu, d)
    Xk <- Psi_true_sqrt %*% Z %*% Sigma_true_sqrt
    acc <- acc + sum(diag(Si_inv %*% t(Xk) %*% Ps_inv %*% Xk))
  }
  acc
}

## Fixed block argument list (same dataset each call)
block_args <- mapply(function(seed, n_i) list(seed=seed, n_i=n_i),
                     block_seeds, block_sizes, SIMPLIFY = FALSE)

## ===== Flip–flop iterations with diagnostics (parallel) =====

# start timer
time_start <- Sys.time()

Tmax <- 50
tol  <- 1e-10
ridge <- 0

Psi_t   <- diag(nu)
Sigma_t <- diag(d)

hist_S_err    <- numeric(Tmax)
hist_P_err    <- numeric(Tmax)
hist_S_delta  <- numeric(Tmax)
hist_P_delta  <- numeric(Tmax)
hist_Kron_err <- numeric(Tmax)
hist_loglik   <- numeric(Tmax)
hist_res_S    <- numeric(Tmax)
hist_res_P    <- numeric(Tmax)

for (t in 1:Tmax) {
  
  ## ---------- Σ update (parallel) ----------
  Y <- mat_invsqrt(Psi_t) %*% W %*% mat_invsqrt(Psi_t)          # Ψ^{-1/2} W Ψ^{-1/2}
  clusterExport(cl, c("Y","Q","block_sum_S"), envir=environment())
  S_list <- parLapply(cl, block_args, block_sum_S)
  S_sum  <- Reduce("+", S_list)
  Sigma_tp1 <- symmetrize(S_sum / (n * tr(W)))
  if (ridge > 0) Sigma_tp1 <- Sigma_tp1 + ridge * diag(d)
  Sigma_tp1 <- psd_project(Sigma_tp1)
  
  ## ---------- Ψ update (parallel) ----------
  Zmask <- mat_invsqrt(Sigma_tp1) %*% U %*% mat_invsqrt(Sigma_tp1)
  clusterExport(cl, c("Zmask","P","block_sum_P"), envir=environment())
  P_list <- parLapply(cl, block_args, block_sum_P)
  P_sum  <- Reduce("+", P_list)
  Psi_tp1 <- symmetrize(P_sum / (n * tr(U)))
  if (ridge > 0) Psi_tp1 <- Psi_tp1 + ridge * diag(nu)
  Psi_tp1 <- psd_project(Psi_tp1)
  
  ## ---------- Identifiability rescaling ----------
  cscale <- d / tr(Sigma_tp1)
  Sigma_tp1 <- Sigma_tp1 * cscale
  Psi_tp1   <- Psi_tp1   / cscale
  
  ## ---------- Residuals of the estimating equations (parallel) ----------
  Y_res <- mat_invsqrt(Psi_tp1) %*% W %*% mat_invsqrt(Psi_tp1)
  clusterExport(cl, c("Y_res","block_sum_S_res"), envir=environment())
  Sres_list <- parLapply(cl, block_args, block_sum_S_res)
  S_res_sum <- Reduce("+", Sres_list)
  R_S <- symmetrize(S_res_sum / (n * tr(W)) - Sigma_tp1)
  hist_res_S[t] <- fro(R_S)
  
  Z_res <- mat_invsqrt(Sigma_tp1) %*% U %*% mat_invsqrt(Sigma_tp1)
  clusterExport(cl, c("Z_res","block_sum_P_res"), envir=environment())
  Pres_list <- parLapply(cl, block_args, block_sum_P_res)
  P_res_sum <- Reduce("+", Pres_list)
  R_P <- symmetrize(P_res_sum / (n * tr(U)) - Psi_tp1)
  hist_res_P[t] <- fro(R_P)
  
  ## ---------- Errors vs truth on identifiable blocks ----------
  if (use_column_mask) {
    hist_S_err[t] <- fro(subblock(Sigma_tp1, obs_cols) - subblock(Sigma_true_sc, obs_cols))
  } else {
    hist_S_err[t] <- fro(Sigma_tp1 - Sigma_true_sc)
  }
  if (use_row_mask) {
    hist_P_err[t] <- fro(subblock(Psi_tp1, obs_rows) - subblock(Psi_true_sc, obs_rows))
  } else {
    hist_P_err[t] <- fro(Psi_tp1 - Psi_true_sc)
  }
  
  ## Kronecker error on identifiable blocks
  Psi_eval_est   <- if (use_row_mask) subblock(Psi_tp1,   obs_rows) else Psi_tp1
  Psi_eval_true  <- if (use_row_mask) subblock(Psi_true,  obs_rows) else Psi_true
  Sigma_eval_est <- if (use_column_mask) subblock(Sigma_tp1, obs_cols) else Sigma_tp1
  Sigma_eval_true<- if (use_column_mask) subblock(Sigma_true,obs_cols) else Sigma_true
  hist_Kron_err[t] <- fro(kronecker(Psi_eval_est, Sigma_eval_est) -
                            kronecker(Psi_eval_true, Sigma_eval_true))
  
  ## Update deltas
  hist_S_delta[t] <- if (t == 1) NA else fro(Sigma_tp1 - Sigma_t)
  hist_P_delta[t] <- if (t == 1) NA else fro(Psi_tp1   - Psi_t)
  
  ## ---------- Log-likelihood (parallel; scale-invariant) ----------
  Ps_inv <- solve(psd_project(Psi_tp1))
  Si_inv <- solve(psd_project(Sigma_tp1))
  clusterExport(cl, c("Ps_inv","Si_inv","block_sum_quad"), envir=environment())
  quad_list <- parLapply(cl, block_args, block_sum_quad)
  quad_sum  <- Reduce("+", quad_list)
  term_det  <- -(n*d/2)*logdet_spd(Psi_tp1) - (n*nu/2)*logdet_spd(Sigma_tp1)
  hist_loglik[t] <- term_det - 0.5 * quad_sum
  
  ## Print a compact iteration log
  cat(sprintf(
    "iter %3d:  ||ΔΣ||=%.2e  ||ΔΨ||=%.2e  ||resΣ||=%.2e  ||resΨ||=%.2e  L=%.6e\n",
    t, ifelse(t==1, NA, hist_S_delta[t]), ifelse(t==1, NA, hist_P_delta[t]),
    hist_res_S[t], hist_res_P[t], hist_loglik[t]
  ))
  
  ## Update for next iteration and early stop
  Sigma_t <- Sigma_tp1; Psi_t <- Psi_tp1
  if (t >= 2 && max(hist_S_delta[t], hist_P_delta[t]) < tol) {
    cat("Converged by update deltas.\n")
    hist_S_err    <- hist_S_err[1:t]
    hist_P_err    <- hist_P_err[1:t]
    hist_S_delta  <- hist_S_delta[1:t]
    hist_P_delta  <- hist_P_delta[1:t]
    hist_Kron_err <- hist_Kron_err[1:t]
    hist_loglik   <- hist_loglik[1:t]
    hist_res_S    <- hist_res_S[1:t]
    hist_res_P    <- hist_res_P[1:t]
    break
  }
}

## Final report
cat("\nFinal estimates (trace-normalized Σ):\n")
print(Sigma_t); print(Psi_t)

if (use_column_mask) {
  cat(sprintf("\nΣ error on observed columns: ||Q Σ Q - Q Σ*_{sc} Q||_F = %.4e\n",
              fro(subblock(Sigma_t, obs_cols) - subblock(Sigma_true_sc, obs_cols))))
} else {
  cat(sprintf("\nΣ error (full): ||Σ - Σ*_{sc}||_F = %.4e\n", fro(Sigma_t - Sigma_true_sc)))
}
if (use_row_mask) {
  cat(sprintf("Ψ error on observed rows: ||P Ψ P - P Ψ*_{sc} P||_F = %.4e\n",
              fro(subblock(Psi_t, obs_rows) - subblock(Psi_true_sc, obs_rows))))
} else {
  cat(sprintf("Ψ error (full): ||Ψ - Ψ*_{sc}||_F = %.4e\n", fro(Psi_t - Psi_true_sc)))
}
cat(sprintf("Kronecker (on identifiable blocks): ||·||_F = %.4e\n", hist_Kron_err[length(hist_Kron_err)]))
cat(sprintf("Moment residuals at the fixed point:  ||resΣ||_F=%.4e  ||resΨ||_F=%.4e\n",
            hist_res_S[length(hist_res_S)], hist_res_P[length(hist_res_P)]))

# end timer
time_end <- Sys.time()
time_spent <- as.numeric(difftime(time_end, time_start, units = "mins"))
cat("Time spent:", round(time_spent, 2), "minutes\n")

## ===== Save plots to disk =====
out_dir <- "C:/Users/fred1/Dropbox/GORS_2025_projects/1_GORS_Stein_matrix_normal/parameter_estimation/graphs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

png(file.path(out_dir, "flipflop_errors.png"), width=1400, height=700, res=150)
par(mfrow=c(1,2))
plot(hist_S_err, type="l", lwd=2, xlab="iteration", ylab="Frobenius error",
     main=if (use_column_mask) "||Q Σ^{(t)} Q - Q Σ*_{sc} Q||_F" else "||Σ^{(t)} - Σ*_{sc}||_F")
plot(hist_P_err, type="l", lwd=2, xlab="iteration", ylab="Frobenius error",
     main=if (use_row_mask) "||P Ψ^{(t)} P - P Ψ*_{sc} P||_F" else "||Ψ^{(t)} - Ψ*_{sc}||_F")
dev.off()

png(file.path(out_dir, "flipflop_deltas_residuals.png"), width=1400, height=700, res=150)
par(mfrow=c(1,2))
plot(hist_S_delta, type="l", lwd=2, xlab="iteration", ylab="Frobenius delta",
     main="||Σ^{(t)} - Σ^{(t-1)}||_F")
lines(hist_P_delta, type="l", lwd=2, lty=2, col="gray40")
legend("topright", c("Sigma delta","Psi delta"), lty=c(1,2), lwd=2, col=c("black","gray40"))
plot(hist_res_S, type="l", lwd=2, xlab="iteration", ylab="Frobenius residual",
     main="moment residuals at iterate")
lines(hist_res_P, type="l", lwd=2, lty=2, col="gray40")
legend("topright", c("res Σ","res Ψ"), lty=c(1,2), lwd=2, col=c("black","gray40"))
dev.off()

png(file.path(out_dir, "flipflop_loglik_kronecker.png"), width=1400, height=700, res=150)
par(mfrow=c(1,2))
plot(hist_loglik, type="l", lwd=2, xlab="iteration", ylab="log-likelihood",
     main="Matrix-normal log-likelihood (scale-invariant)")
plot(hist_Kron_err, type="l", lwd=2, xlab="iteration", ylab="Frobenius error",
     main="Kronecker error on identifiable blocks")
dev.off()
