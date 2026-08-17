### List of estimator functions

source("./scripts/utils.R")

## Difference-in-means estimator
# Data: Tibble
# Alpha: Allowable type I error
estimate_DM <- function(data, alpha = 0.05, true_value = 1) {
  data %>% group_by(A) %>% 
    summarise(
      mean_outcome = mean(Y, na.rm = TRUE),
      var_norm = var(Y) / n()
    ) %>% 
    summarise(
      estimator = "DM",
      Qform = NA,
      gform = NA,
      est_init = NA,
      est = mean_outcome[2] - mean_outcome[1],
      se = sqrt(var_norm[2] + var_norm[1]),
      ci_lower = est - qnorm(1-alpha) * se,
      ci_upper = est + qnorm(1-alpha) * se,
      bias = est - true_value,
      bias.to.se = bias / se
    ) %>% 
    as_tibble()
}

## Ordinary Least Squares estimator
# Data: Tibble
# Alpha: Allowable type I error
estimate_OLS <- function(data, alpha = 0.05, true_value = 1,
                         Qform = as.formula(paste("Y ~ A +", paste(covariates, collapse = " + ")))
                         )
  {
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  # OLS formula Y ~ A + W1 + ... + Wk
  #formula_OLS <- as.formula(paste("Y ~ A +", paste(covariates, collapse = " + "))) #, "+ A:(", paste(covariates, collapse = " + "), ")"))
  model <- lm(Qform, data = data)
  coefs_A <- coef(summary(model))["A",]
  return(
    tibble(
      estimator = "OR",
      Qform = paste("Y ~",as.character(Qform)[3]),
      gform = NA,
      est_init = NA,
      est = coefs_A["Estimate"],
      se = coefs_A["Std. Error"],
      ci_lower = confint(model, "A", level = 1-alpha)[1],
      ci_upper = confint(model, "A", level = 1-alpha)[2],
      bias = est - true_value,
      bias.to.se = bias / se
    )
  )
}


## Ordinary Least Squares substitution estimator
# Data: Tibble
# Alpha: Allowable type I error
estimate_OLS_sub <- function(data, alpha = 0.05, true_value = 1,
                             Qform = as.formula(paste("Y ~ A +", paste(covariates, collapse = " + ")))
                             )
  {
  
  estimate <- OLS_sub(data, Qform) 
  boot_res <- boot(data,
                   statistic = function(data, indices){ OLS_stat(data, indices, Qform)},
                   R = 1000)

  return(
    tibble(
      estimator = "OR_sub",
      Qform = paste("Y ~",as.character(Qform)[3]),
      gform = NA,
      est_init = NA,
      est = estimate,
      se = sd(boot_res$t),
      ci_lower = quantile(boot_res$t, probs = alpha/2, na.rm = TRUE),
      ci_upper = quantile(boot_res$t, probs = 1-alpha/2, na.rm = TRUE),
      bias = est - true_value,
      bias.to.se = bias / se
    )
  )
}


## Inverse Propensity Weighting estimator (survey package)
# Data: Tibble
# Alpha: Allowable type I error
# gform: Propensity score formula, e.g., A ~ W1*W2
estimate_IPW <- function(data, alpha = 0.05, true_value = 1,
                         gform = as.formula(paste("A ~", paste(
                             colnames(data)[grepl("^W", colnames(data))], collapse = " + "))
                             )
                         ) {
  # covariates <- colnames(data)[grepl("^W", colnames(data))]
  propensity_model <- glm2(gform,
                           family = binomial(link = "logit"),
                           data = data)
  
  data$weights <- ifelse(data$A == 1, 
                         1 / predict(propensity_model, type = "response", newdata = data),
                         1 / (1 - predict(propensity_model, type = "response", newdata = data))
                         )
  
  ipw_design <- svydesign(ids = ~1, weights = ~weights, data = data)
  
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  formula_IPW <- as.formula("Y ~ A")
  #formula_IPW <- as.formula(paste("Y ~ A +", paste(
  #  colnames(data)[grepl("^W", colnames(data))], collapse = " + ")))
  ipw_result <- svyglm(formula_IPW, design = ipw_design)
  
  coef_ipw <- coef(ipw_result)
  se_ipw <- sqrt(diag(vcov(ipw_result)))
  ci_ipw <- coef_ipw["A"] + c(-qnorm(1-alpha), qnorm(1-alpha)) * se_ipw["A"]

  return(
    tibble(
      estimator = "IPW",
      Qform = NA,
      gform = paste("A ~",as.character(gform)[3]),
      est_init = NA,
      est = coef_ipw["A"],
      se = se_ipw["A"],
      ci_lower = ci_ipw[1],
      ci_upper = ci_ipw[2],
      bias = est - true_value,
      bias.to.se = bias / se
    )
  )
}


## Inverse Propensity Weighting estimator (by hand)
# Data: Tibble
# Alpha: Allowable type I error
estimate_IPW_hand <- function(data, alpha = 0.05, true_value = 1,
                              gform = as.formula(paste("A ~", paste(
                                colnames(data)[grepl("^W", colnames(data))], collapse = " + "))
                              )
) {
  # covariates <- colnames(data)[grepl("^W", colnames(data))]
  propensity_model <- glm2(gform,
                           family = binomial(link = "logit"),
                           data = data)
  
  data$weights <- ifelse(data$A == 1, 
                         1 / predict(propensity_model, type = "response", newdata = data),
                         1 / (1 - predict(propensity_model, type = "response", newdata = data))
  )
  
  est_hand <- data %>% summarise(mean(Y * (2*A-1) * weights)) %>% as.numeric()
  
  return(
    tibble(
      estimator = "IPW_hand",
      est = est_hand,
      se = NA,
      ci_lower = NA,
      ci_upper = NA,
      bias = est - true_value,
      bias.to.se = bias / se
    )
  )
}


## TMLE estimator
# Data: Tibble
# Alpha: Allowable type I error
# Qform
# gform
# TMLE specs: Weighted, single_H, bounded_continuous_outcome
estimate_TMLE <- function(data, alpha = 0.05, true_value = 1,
                          Qform = as.formula(paste("Y ~ A +", paste(covariates, collapse = " + "))),
                          gform = as.formula(paste("A ~", paste(
                            colnames(data)[grepl("^W", colnames(data))], collapse = " + "))),
                          weight = TRUE, use_single_H = FALSE, transforms = TRUE) {
  
  if(transforms){
    a <- min(data$Y)
    b <- max(data$Y)
    data$Y <- (data$Y - a) / (b - a)
  }
  
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  
  # Outcome regression
  Q_init <- glm(Qform, family = gaussian, data = data)
  data$Q0W <- predict(Q_init, newdata = data.frame(A = 0, data[covariates]))
  data$Q1W <- predict(Q_init, newdata = data.frame(A = 1, data[covariates]))
  data$Q <- predict(Q_init, type = "response")
  
  # Treatment mechanism
  g_init <- glm(gform, family = binomial, data = data)
  data$g1W <- predict(g_init, type = "response")
  
  data$k_Q1 <- if (weight) 1 else 1 / data$g1W
  data$k_Q0 <- if (weight) 1 else 1 / (1 - data$g1W)
  
  data$H <-if (weight) data$A - (1 - data$A) else (data$A / data$g1W) - ((1 - data$A) / (1 - data$g1W))
  data$weights_single_H <- if (weight) ifelse(data$A == 1, 1 / data$g1W, 1 / (1 - data$g1W)) else 1
  data$H1 <- if(weight) data$A else data$A / data$g1W
  data$H0 <- if(weight) (1  - data$A) else ((1 - data$A) / (1 - data$g1W))
  data$weights_H1 <- if (weight) 1 / data$g1W else 1
  data$weights_H0 <- if (weight) 1 / (1 - data$g1W) else 1
  
  if (use_single_H) {
    if (transforms) {
      second_stage <- glm2(data$Y ~ -1 + H + offset(logit(data$Q)), 
                           family = quasibinomial, data = data, weights = weights_single_H)
      epsilon <- coef(second_stage)[1]
      Q_star1 <- expit(logit(data$Q1W) + epsilon * data$k_Q1)
      Q_star0 <- expit(logit(data$Q0W) - epsilon * data$k_Q0)
    } else {
      second_stage <- glm2(data$Y ~ -1 + H + offset(data$Q), 
                           family = gaussian, data = data, weights = weights_single_H)
      epsilon <- coef(second_stage)[1]
      Q_star1 <- data$Q1W + epsilon * data$k_Q1
      Q_star0 <- data$Q0W - epsilon * data$k_Q0
    }
  } else {
    if (transforms) {
      second_stage_1 <- glm2(data$Y ~ -1 + H1 + offset(logit(data$Q)), 
                             family = quasibinomial, data = data, weights = weights_H1)
      second_stage_0 <- glm2(data$Y ~ -1 + H0 + offset(logit(data$Q)), 
                             family = quasibinomial, data = data, weights = weights_H0)
      epsilon_1 <- coef(second_stage_1)[1]
      epsilon_0 <- coef(second_stage_0)[1]
      Q_star1 <- expit(logit(data$Q1W) + epsilon_1 * data$k_Q1)
      Q_star0 <- expit(logit(data$Q0W) + epsilon_0 * data$k_Q0)
    } else {
      second_stage_1 <- glm2(data$Y ~ -1 + H1 + offset(data$Q), 
                             family = gaussian, data = data, weights = weights_H1)
      second_stage_0 <- glm2(data$Y ~ -1 + H0 + offset(data$Q), 
                             family = gaussian, data = data, weights = weights_H0)
      epsilon_1 <- coef(second_stage_1)[1]
      epsilon_0 <- coef(second_stage_0)[1]
      Q_star1 <- data$Q1W + epsilon_1 * data$k_Q1
      Q_star0 <- data$Q0W + epsilon_0 * data$k_Q0
    }
  }
  
  if (transforms) {
    ATE_pre <- mean(data$Q1W - data$Q0W) * (b - a)
    ATE <- mean(Q_star1 - Q_star0) * (b - a)
  }else{
    ATE_pre <- mean(data$Q1W - data$Q0W)
    ATE <- mean(Q_star1 - Q_star0)
  }
  
  if  (transforms) {
    ATE_n <- ATE / (b - a)
  }else{
    ATE_n <- ATE 
  }
  
  data$H_n <- (data$A / data$g1W) - ((1 - data$A) / (1 - data$g1W))
  data$H_1 <- data$A / data$g1W
  data$H_0 <- (1 - data$A) / (1 - data$g1W)
  
  if(use_single_H){
    IF <- data$H_n * (data$Y - (data$A * Q_star1 + (1 - data$A) * Q_star0)) + (Q_star1 - Q_star0) - ATE_n
    IF_1 <- mean(data$H_n * (data$Y - (data$A * Q_star1 + (1 - data$A) * Q_star0)))
    IF_0 <- mean((Q_star1 - Q_star0) - ATE_n)
  }else{
    IF <- data$H_1* (data$Y - Q_star1) - data$H_0 * (data$Y - Q_star0) + (Q_star1 - Q_star0) - ATE_n
    IF_1 <- mean(data$H_1* (data$Y - Q_star1) - data$H_0 * (data$Y - Q_star0))
    IF_0 <- mean((Q_star1 - Q_star0) - ATE_n)
  }
  
  n <- nrow(data)
  SE <- if(transforms) sqrt(var(IF) / n) * (b - a) else sqrt(var(IF) / n)
  z <- qnorm(1 - alpha / 2)
  v <- if(transforms) (var(IF) / n) * (b - a)^2 else var(IF) / n #var(IF)/n THIS WAS A TYPO
  CI_lower <- ATE - z * SE
  CI_upper <- ATE + z * SE
  #MSE <- bias^2 + v
  
  #return(list(estimate = ATE, var = v, bias = bias, se = SE, mse = MSE, ci_lower = CI_lower, ci_upper = CI_upper))
  
  return(
    tibble(
      estimator = "TMLE",
      Qform = paste("Y ~",as.character(Qform)[3]),
      gform = paste("A ~",as.character(gform)[3]),
      est_init = ATE_pre,
      est = ATE,
      se = SE,
      ci_lower = CI_lower,
      ci_upper = CI_upper,
      bias = est - true_value,
      bias.to.se = bias / se
    )
  )
}

