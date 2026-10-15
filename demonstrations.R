library(sandwich)
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else "demo/suggestions.csv"
out <- if (length(args) >= 2) args[2] else "results/draft_6_demonstrations"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

effect_table <- function(fit, L, V, effects, se_type) {
  estimate <- drop(L %*% coef(fit))
  se <- sqrt(diag(L %*% V %*% t(L)))
  margin <- qt(.975, df.residual(fit)) * se
  data.frame(effect = effects, se_type, estimate, se,
             lower = estimate - margin, upper = estimate + margin,
             p_value = 2 * pt(-abs(estimate / se), df.residual(fit)))
}

to01 <- function(x) as.integer(tolower(trimws(as.character(x))) %in%
                                c("1", "true", "t", "yes", "y"))
hs <- read.csv(data_file) %>%
  filter(user.index == 8, to01(avail) == 1, to01(connect) == 1) %>%
  mutate(W = case_when(to01(send) == 0 ~ 0L,
                        to01(send.active) == 1 ~ 1L,
                        to01(send.sedentary) == 1 ~ 2L,
                        TRUE ~ NA_integer_)) %>%
  filter(!is.na(W)) %>%
  arrange(decision.index.nogap, decision.index) %>%
  mutate(W_lag = lag(W), Y = gfsteps30) %>%
  filter(!is.na(W_lag), is.finite(Y))

p <- c(.4, .3, .3)
Sigma <- diag(p[-1]) - tcrossprod(p[-1])
M <- solve(Sigma)
D <- with(hs, cbind(as.numeric(W == 1), as.numeric(W == 2)))
D_lag <- with(hs, cbind(as.numeric(W_lag == 1), as.numeric(W_lag == 2)))
Z <- sweep(D, 2, p[-1]) %*% M
Z_lag <- sweep(D_lag, 2, p[-1]) %*% M
X <- cbind(Z, Z_lag)
fit <- lm(hs$Y ~ X - 1)
L <- kronecker(diag(2), M)
effects <- c("Immediate: walking vs no suggestion",
             "Immediate: anti-sedentary vs no suggestion",
             "Lag-1: walking vs no suggestion",
             "Lag-1: anti-sedentary vs no suggestion")
heartsteps <- bind_rows(
  effect_table(fit, L, vcovHC(fit, type = "HC3"), effects, "HC3"),
  effect_table(fit, L, NeweyWest(fit, lag = floor(nobs(fit)^.25),
                                prewhite = FALSE, adjust = TRUE), effects, "HAC"))
cat("HeartSteps participant 8; analysis occasions:", nrow(hs), "\n")
print(heartsteps, row.names = FALSE, digits = 6)
write.csv(heartsteps, file.path(out, "table_3.csv"), row.names = FALSE)

jacob <- data.frame(
  session = 1:16,
  Y = c(26.7, 33.1, 16.6, 26.9, 20.1, 13.2, 16.7, 3.1,
        6.6, 0, 36.4, 13.1, 9.8, 9.8, 0, 0),
  W = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1))
jacob$p_t <- .5
for (t in 3:nrow(jacob)) {
  if (jacob$W[t - 1] == jacob$W[t - 2]) jacob$p_t[t] <- 1 - jacob$W[t - 1]
}
jacob <- jacob %>%
  mutate(v_t = p_t * (1 - p_t), Z = 0,
         assignment = ifelse(v_t > 0, "Randomized", "Deterministic"),
         condition = ifelse(W == 0, "Control (A)", "Treatment (B)"))
eligible <- jacob$v_t > 0
jacob$Z[eligible] <- with(jacob[eligible, ], (W - p_t) / v_t)
jacob$Z_lag <- lag(jacob$Z, default = 0)
v_lag <- lag(jacob$v_t, default = 0)
c0 <- mean(1 / jacob$v_t[eligible])
c1 <- mean(1 / v_lag[v_lag > 0])
fit <- lm(Y ~ Z + Z_lag + p_t, data = jacob)
L <- cbind(0, diag(c(c0, c1)), 0)
jacob_effects <- effect_table(
  fit, L, NeweyWest(fit, lag = floor(nobs(fit)^.25), prewhite = FALSE, adjust = FALSE),
  c("Immediate", "Lag-1"), "HAC")
print(jacob_effects, row.names = FALSE, digits = 6)
write.csv(jacob_effects, file.path(out, "jacob_effects.csv"), row.names = FALSE)
write.csv(select(jacob, session, assignment, condition, Z),
           file.path(out, "table_4.csv"), row.names = FALSE)

p <- ggplot(jacob, aes(session, Y, shape = condition, linetype = condition)) +
  geom_line() + geom_point(size = 2.5) +
  scale_x_continuous(breaks = 1:16) +
  labs(x = "Session", y = "Intervals with problem behavior (%)",
       shape = NULL, linetype = NULL) +
  theme_classic(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(out, "figure_8.pdf"), p, width = 8, height = 5)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
