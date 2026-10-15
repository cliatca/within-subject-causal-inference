library(sandwich)
library(dplyr)
library(tidyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args) >= 1) as.integer(args[1]) else 1000L
out <- if (length(args) >= 2) args[2] else "results/draft_6_simulations"
cores <- if (length(args) >= 3) as.integer(args[3]) else 4L
stopifnot(n_rep > 0, cores > 0)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

grid1 <- expand.grid(beta1 = c(0, 1), delta_B0 = c(0, 3),
                     delta_C0 = c(0, 2), delta_B1 = c(0, 1),
                     delta_C1 = c(0, 1), T = c(10, 20, 50, 100, 150))
grid2 <- expand.grid(beta1 = c(0, 1), delta0 = c(0, 2, 3, 4),
                     delta1 = c(0, 1), T = c(15, 30, 60))
p <- c(.4, .3, .3)
Sigma <- diag(p[-1]) - tcrossprod(p[-1])
M <- solve(Sigma)
L1 <- kronecker(diag(2), M)

infer <- function(fit, L, truth, hac_adjust) {
  k <- length(truth)
  ans <- matrix(NA_real_, 2 * k, 4,
                dimnames = list(NULL, c("estimate", "se", "covered", "rejected")))
  if (fit$rank < length(coef(fit)) || df.residual(fit) <= 0) return(ans)
  estimate <- drop(L %*% coef(fit))
  ans[, "estimate"] <- rep(estimate, 2)
  for (j in 1:2) {
    if (j == 2 && any(hatvalues(fit) >= 1 - 1e-10)) next
    V <- tryCatch(if (j == 1) {
      NeweyWest(fit, lag = floor(nobs(fit)^.25), prewhite = FALSE,
                adjust = hac_adjust)
    } else vcovHC(fit, type = "HC3"), error = function(e) NULL)
    if (is.null(V)) next
    variance <- diag(L %*% V %*% t(L))
    se <- sqrt(pmax(variance, 0))
    se[!is.finite(se) | se <= 0 | variance < 0] <- NA_real_
    margin <- qt(.975, df.residual(fit)) * se
    rows <- (j - 1) * k + seq_len(k)
    ans[rows, "se"] <- se
    ans[rows, "covered"] <- abs(estimate - truth) <= margin
    ans[rows, "rejected"] <- abs(estimate) > margin
  }
  ans
}

mean_valid <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

summarize_runs <- function(a, truth, effects, model) {
  k <- length(truth)
  bind_rows(lapply(seq_len(2 * k), function(j) {
    h <- (j - 1) %% k + 1
    estimate <- a[, j, "estimate"]
    bias <- mean_valid(estimate) - truth[h]
    rejection <- mean_valid(a[, j, "rejected"])
    data.frame(model, se_type = if (j <= k) "HAC" else "HC3",
               effect = effects[h], truth = truth[h], n_rep,
               n_estimated = sum(is.finite(estimate)),
               n_inference = sum(is.finite(a[, j, "covered"])),
               estimate = mean_valid(estimate), bias,
               relative_bias = if (truth[h] == 0) NA_real_ else bias / truth[h],
               rmse = sqrt(mean_valid((estimate - truth[h])^2)),
               coverage = mean_valid(a[, j, "covered"]),
               type1_error = if (truth[h] == 0) rejection else NA_real_,
               power = if (truth[h] != 0) rejection else NA_real_)
  }))
}

run_three_arm <- function(i) {
  set.seed(2026 + i)
  g <- grid1[i, ]
  T <- g$T
  truth <- unlist(g[c("delta_B0", "delta_C0", "delta_B1", "delta_C1")],
                  use.names = FALSE)
  a <- array(NA_real_, c(n_rep, 8, 4),
             dimnames = list(NULL, NULL, c("estimate", "se", "covered", "rejected")))
  for (r in seq_len(n_rep)) {
    W <- sample(0:2, T, replace = TRUE, prob = p)
    D <- cbind(as.numeric(W == 1), as.numeric(W == 2))
    D_lag <- rbind(c(0, 0), D[-T, ])
    Y <- drop(g$beta1 * log(seq_len(T)) + cbind(D, D_lag) %*% truth + rnorm(T))
    Z <- sweep(D, 2, p[-1]) %*% M
    X <- cbind(Z[-1, ], Z[-T, ])
    fit <- lm(Y[-1] ~ X - 1)
    a[r, , ] <- infer(fit, L1, truth, hac_adjust = TRUE)
  }
  s <- summarize_runs(a, truth, c("Immediate B vs A", "Immediate C vs A",
                                  "Lag-1 B vs A", "Lag-1 C vs A"), "Multiple-arm")
  bind_cols(g[rep(1, nrow(s)), ], s)
}

run_restricted <- function(i) {
  set.seed(2026 + nrow(grid1) + i)
  g <- grid2[i, ]
  T <- g$T
  truth <- c(g$delta0, g$delta1)
  a <- b <- array(NA_real_, c(n_rep, 4, 4),
                  dimnames = list(NULL, NULL, c("estimate", "se", "covered", "rejected")))
  for (r in seq_len(n_rep)) {
    W <- integer(T)
    p_t <- rep(.5, T)
    for (t in seq_len(T)) {
      if (t > 2 && W[t - 1] == W[t - 2]) p_t[t] <- 1 - W[t - 1]
      W[t] <- if (p_t[t] == .5) rbinom(1, 1, .5) else p_t[t]
    }
    Y <- g$beta1 * log(seq_len(T)) + g$delta0 * W +
      g$delta1 * c(0, W[-T]) + rnorm(T)
    v_t <- p_t * (1 - p_t)
    Z <- numeric(T)
    eligible <- v_t > 0
    Z[eligible] <- (W[eligible] - p_t[eligible]) / v_t[eligible]
    Z_lag <- c(0, Z[-T])
    v_lag <- c(0, v_t[-T])
    c0 <- mean(1 / v_t[eligible])
    c1 <- mean(1 / v_lag[v_lag > 0])
    fit0 <- lm(Y ~ Z + Z_lag - 1)
    fit1 <- lm(Y ~ Z + Z_lag + p_t)
    L0 <- diag(c(c0, c1))
    L <- cbind(0, L0, 0)
    a[r, , ] <- infer(fit0, L0, truth, hac_adjust = FALSE)
    b[r, , ] <- infer(fit1, L, truth, hac_adjust = FALSE)
  }
  effects <- c("Immediate", "Lag-1")
  s <- bind_rows(summarize_runs(a, truth, effects, "Unadjusted"),
                 summarize_runs(b, truth, effects, "Probability-adjusted"))
  bind_cols(g[rep(1, nrow(s)), ], s)
}

run_grid <- function(n, fun) {
  if (.Platform$OS.type == "windows" || cores == 1) {
    bind_rows(lapply(seq_len(n), fun))
  } else {
    bind_rows(parallel::mclapply(seq_len(n), fun, mc.cores = cores,
                                mc.set.seed = FALSE))
  }
}

message("Simulation I: ", nrow(grid1), " conditions, ", n_rep, " replications each")
s1 <- run_grid(nrow(grid1), run_three_arm)
write.csv(s1, file.path(out, "simulation_I.csv"), row.names = FALSE)
message("Simulation II: ", nrow(grid2), " conditions, ", n_rep, " replications each")
s2 <- run_grid(nrow(grid2), run_restricted)
write.csv(s2, file.path(out, "simulation_II.csv"), row.names = FALSE)

point_plot <- function(d, metric, filename, label) {
  d <- d %>% filter(se_type == "HAC", is.finite(.data[[metric]])) %>%
    group_by(model, effect, T) %>%
    summarise(value = mean(.data[[metric]]),
              lo = quantile(.data[[metric]], .25),
              hi = quantile(.data[[metric]], .75), .groups = "drop")
  p <- ggplot(d, aes(T, value, group = model, shape = model, linetype = model)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 1) +
    geom_line() + geom_point(size = 2) + facet_wrap(~effect, ncol = 2) +
    scale_x_continuous(breaks = sort(unique(d$T))) +
    labs(x = "Number of occasions", y = label, shape = NULL, linetype = NULL) +
    theme_classic(base_size = 12) + theme(legend.position = "bottom")
  if (metric == "relative_bias") p <- p + geom_hline(yintercept = 0, linetype = 3)
  ggsave(file.path(out, filename), p, width = 9, height = 6)
}

inference_plot <- function(d, filename) {
  d <- d %>% group_by(se_type, effect, T) %>%
    summarise(across(c(coverage, type1_error, power), mean_valid), .groups = "drop") %>%
    pivot_longer(c(coverage, type1_error, power), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, c("coverage", "type1_error", "power"),
                            c("Coverage", "Type I error", "Power")))
  reference <- data.frame(metric = factor(c("Coverage", "Type I error"),
                                          levels = levels(d$metric)), value = c(.95, .05))
  p <- ggplot(d, aes(T, value, group = se_type, shape = se_type, linetype = se_type)) +
    geom_hline(data = reference, aes(yintercept = value), inherit.aes = FALSE,
               linetype = 3, colour = "grey50") +
    geom_line() + geom_point(size = 2) + facet_grid(metric ~ effect, scales = "free_y") +
    scale_x_continuous(breaks = sort(unique(d$T))) +
    labs(x = "Number of occasions", y = NULL, shape = NULL, linetype = NULL) +
    theme_classic(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(out, filename), p, width = 11, height = 7)
}

point_plot(s1, "relative_bias", "figure_2.pdf", "Relative bias")
point_plot(s1, "rmse", "figure_3.pdf", "RMSE")
inference_plot(s1, "figure_4.pdf")
point_plot(s2, "relative_bias", "figure_5.pdf", "Relative bias")
point_plot(filter(s2, model == "Probability-adjusted"), "rmse", "figure_6.pdf", "RMSE")
inference_plot(filter(s2, model == "Probability-adjusted"), "figure_7.pdf")
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
message("Results saved to ", out)
