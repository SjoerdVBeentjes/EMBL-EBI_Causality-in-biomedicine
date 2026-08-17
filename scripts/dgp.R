### Data Generating Processes

## Study 1: low-dimensional, highly correlated covariates
# Simulation 1a
sim1a <- function(n=1000,
                  mu=c(0,1),
                  sigma=matrix(c(1,0.5,0.5,1), nrow = 2)) {
  # Covariates
  W <- mvrnorm(n, mu = mu, Sigma = sigma)
  W1 <- W[,1]
  W2 <- W[,2]
  
  # Treatment mechanism
  g0 <- expit(0.5 - 1.5 * W1 + 0.5 * W2)
  A <- rbinom(n, 1, g0)
  
  # Outcome mechanism
  Q0 <- 1 + A - 0.7 * W1 + 0.3*exp(-W1*W2)
  epsilon <- rnorm(n, mean = 0, sd = 1) 
  Y <- Q0 + epsilon
  
  data1 <- tibble(W1 = W1, W2 = W2, A = A, Y = Y)
  return(data1)
}


# Simulation 1b
sim1b <- function(n=1000,
                  mu=c(0,1),
                  sigma=matrix(c(1,0.5,0.5,1), nrow = 2)) {
  # Covariates
  W <- mvrnorm(n, mu = mu, Sigma = sigma)
  W1 <- W[,1]
  W2 <- W[,2]
  
  # Treatment mechanism
  g0 <- expit(0.5 - 1.5 * W1 * W2)
  A <- rbinom(n, 1, g0)
  
  # Outcome mechanism
  Q0 <- 1 + A - 0.5 * W1 - 0.25*W1*W2
  epsilon <- rnorm(n, mean = 0, sd = 1) 
  Y <- Q0 + epsilon
  
  data1 <- tibble(W1 = W1, W2 = W2, A = A, Y = Y)
  return(data1)
}

# Simulation 1c: Random treatment assignment
sim1c <- function(n=1000,
                  mu=c(0,1),
                  sigma=matrix(c(1,0.5,0.5,1), nrow = 2),
                  prob_treatment = 0.5) {
  # Covariates
  W <- mvrnorm(n, mu = mu, Sigma = sigma)
  W1 <- W[,1]
  W2 <- W[,2]
  
  # Treatment mechanism
  #g0 <- expit(0.5 - 1.5 * W1 + 0.5 * W2)
  A <- rbinom(n, 1, prob_treatment)
  
  # Outcome mechanism
  Q0 <- 1 + A - 0.7 * W1 + 0.3*exp(-W1*W2)
  epsilon <- rnorm(n, mean = 0, sd = 1) 
  Y <- Q0 + epsilon
  
  data1 <- tibble(W1 = W1, W2 = W2, A = A, Y = Y)
  return(data1)
}


## Study 2: highly correlated covariates
# Simulation 2a
sim2 <- function(n=1000){
  # Covariates
  W1 <- rbinom(n, 1, 0.5)
  W2 <- rbinom(n, 1, 0.5)
  W3 <- rbinom(n, 1, 0.5)
  W4 <- rbinom(n, 1, 0.2 + 0.5 * W1)
  W5 <- rbinom(n, 1, 0.05 + 0.3 * W1 + 0.1 * W2 + 0.05 * W3 + 0.4 * W4)
  W6 <- rbinom(n, 1, 0.2 + 0.6 * W5)
  W7 <- rbinom(n, 1, 0.5 + 0.2 * W3)
  W8 <- rbinom(n, 1, 0.1 + 0.2 * W2 + 0.3 * W6 + 0.1 * W7)
  
  # Treatment mechanism
  logit_g0 <- -0.05 + 0.1*W1 + 0.2*W2 + 0.2*W3 - 0.02*W4 - 0.6*W5 - 0.2*W6 - 0.1*W7
  g0 <- expit(-logit_g0) 
  A <- rbinom(n, 1, g0)
  
  # Outcome mechanism
  Q0 <- 10 + A + W1 + W2 + W4 + 2*W6+ W7
  epsilon <- rnorm(n, mean = 0, sd = 1) 
  Y <- Q0 + epsilon
  
  data <- tibble(W1, W2, W3, W4, W5, W6, W7, W8, A, Y)
  data$g0 <- g0
  
  return(data)
}
