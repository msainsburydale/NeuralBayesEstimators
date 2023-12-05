library("NeuralEstimators")
library("JuliaConnectoR")
library("dplyr")
library("ggplot2")
library("reshape2")

RNGversion("3.6.0")
set.seed(1)

# ---- Construct the curves ----

sigma <- 0.05
n <- 100 # number of observations in each replicate
x <- seq(-1, 1, length.out = n)
X <- cbind(1, x)
error <- function() rnorm(n, sd = sigma)


y1 <- exp(x) + error()
y2 <- exp(-x)  + error()
y3 <- x^3 + error()
y4 <- -x^3 + error()
df <- data.frame(x, y1, y2, y3, y4)
df_long <- melt(df, id.vars = "x", value.name = "y", variable.name = "curve")

ggplot(df_long) +
  geom_point(aes(x, y), size = 0.5) +
  facet_wrap(curve~., nrow = 1) +
  theme_bw() +
  theme(
    legend.position = "none",
  )

# ---- Compute the OLS estimates for each curve ----

(OLS1 <- coef(lm(y1~x)))
(OLS2 <- coef(lm(y2~x)))
(OLS3 <- coef(lm(y3~x)))
(OLS4 <- coef(lm(y4~x)))


# Rather than OLS, use the posterior median of beta
# For simplicity, just assume standard multivariate normal prior
posteriorprecision <- function(sigma, S, X) (sigma^-2) * t(X)%*%X #+ solve(S)
posteriormean <- function(sigma, S, X, y) c((sigma^-2) * solve(posteriorprecision(sigma, S, X)) %*% t(X) %*% y)

(OLS1 <- posteriormean(sigma, S, X, y1))
(OLS2 <- posteriormean(sigma, S, X, y2))
(OLS3 <- posteriormean(sigma, S, X, y3))
(OLS4 <- posteriormean(sigma, S, X, y4))


# ---- Construct the prior ----

prior <- function(K) {
  b0 <- rnorm(K)
  b1 <- rnorm(K)
  theta <- matrix(c(b0, b1), byrow = TRUE, ncol = K)
  return(theta)
}
theta_train <- prior(50000)
theta_val   <- prior(5000)


# ---- Train the neural estimator ----

m <- 1

simulate <- function(theta_set) {
  lapply(1:ncol(theta_set), function(k) {
    theta <- theta_set[, k]
    Z <- theta[1] + theta[2] * x + error()
    dim(Z) <- c(n, 1, m)
    Z
  })
}


Z_train <- simulate(theta_train)
Z_val   <- simulate(theta_val)

estimator <- juliaLet('

  w = 128   # number of neurons in each layer
  p = 2    # number of parameters in the statistical model

  using Flux
  psi = Chain(Dense(n, w, relu), Dense(w, w, relu), Dense(w, w, relu))
  phi = Chain(Dense(w, w, relu), Dense(w, w, relu), Dense(w, p), Flux.flatten)

  using NeuralEstimators
  estimator = DeepSet(psi, phi)
', n = as.integer(n))


estimator <- train(
  estimator,
  theta_train = theta_train,
  theta_val   = theta_val,
  Z_train = Z_train,
  Z_val   = Z_val,
  epochs = 50L
)


# ---- Obtain the neural estimates ----


neuralestimate <- function(estimator, y) {
  dim(y) <- c(n, 1, m)
  estimate(estimator, list(y))
}

(neural1 <- neuralestimate(estimator, y1))
(neural2 <- neuralestimate(estimator, y2))
(neural3 <- neuralestimate(estimator, y3))
(neural4 <- neuralestimate(estimator, y4))

pred <- function(thetahat, x) thetahat[1] + thetahat[2] * x

tmp <- df_long
df_long$neural <- c(sapply(list(neural1, neural2, neural3, neural4), pred, x))
df_long$OLS <- c(sapply(list(OLS1, OLS2, OLS3, OLS4), pred, x))
df_long$y <- NULL
df_long <- melt(df_long, id.vars = c("x", "curve"))

figure <- ggplot(df_long) +
    geom_point(data = tmp, aes(x, y), size = 0.25) +
  geom_line(aes(x, value, colour = variable)) +
  facet_wrap(curve~., nrow = 1) +
  scale_colour_discrete(labels = c(
    "y" = expression(f(x)),
    "neural" = expression(bold(X) * hat(bold(beta))[neural]),
    "OLS" = expression(bold(X) * hat(beta)[analytic])
  )) +
  labs(x = expression(italic(x)),
       y = expression(italic(Z)),
       colour = "") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.text.align = 0
  )


ggsave(
  figure,
  file = "misspecification.pdf",
  width = 8, height = 3, path = file.path("img", "Univariate"), device = "pdf"
)
