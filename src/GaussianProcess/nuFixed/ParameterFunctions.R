suppressMessages({
library("fields")
library("parallel")
library("abind")
library("testthat")
library("reticulate")
})

# d <- document::document(file_name = "~/Dropbox/SpatialDeepSets/src/GaussianProcess/prep_for_Julia.R", check_package = FALSE, output_directory = "./")

#' Compute the Cholesky factors of a Matern covariance matrix for fixed
#' smoothness parameter and an arbitrary number of range parameters.
#'
#' Under the working model as described by Gerber and Nychka (2021), the
#' covariance matrix, \eqn{\Sigma}, of the latent spatial process
#' depends only on a smoothness parameter, \eqn{\nu}, and a range parameter,
#' \eqn{\theta}. This function computes the Cholesky factors of \eqn{\Sigma}
#' for a fixed smoothness parameter and an arbitrary number of range parameters.
#'
#' @param nu Positive scalar representing the Matern smoothness parameter, \eqn{\nu}.
#' @param thetaGrid Vector of positive scalars representing Matern range parameters, \eqn{\theta}.
#' @param S Matrix of spatial locations.
compute_Cholesky_Factors <- function(nu, thetaGrid, S) {

  chols <- lapply(thetaGrid, function(theta) {
    sigma <- stationary.cov(
      S, S,
      Covariance = "Matern",
      theta = theta,
      smoothness = nu
    )

    L <- t(chol(sigma))

    return(L)
  })

  ## Bind the list of Cholesky factors into a 3-dimensional array
  chols <- do.call("abind", c(chols, along = 3)) # array of size (N^2, N^2, nTheta)
  dimnames(chols) <- NULL

  return(chols)
}

#' Prepare the parameter configurations (\eqn{\lambda}, \eqn{\theta}, \eqn{\nu}) and the
#' associated Cholesky factors for a fixed smoothness parameter, \eqn{\nu}.
#'
#' @inheritParams compute_Cholesky_Factors
#' @param lambdaGrid The complete grid of values of \eqn{\lambda} used for
#' each value of \eqn{\lambda} in \code{thetaGrid}; see details.
#' @param nLambda The number of values of \eqn{\lambda} to use
#' with each value of \eqn{\theta} in \code{thetaGrid}; applicable only if
#' \code{is.null(lambdaGrid)}.
#' @param dfRange Range of EDF that \code{lambdaGrid} should cover; see details.
#' @param rangeLambda Range of values that \eqn{\lambda} should take; see details.
#'
#' @details If \code{lambdaGrid} is \code{NULL}, it will be constructed internally by
#' a method that accounts for the EDF (if \code{dfRange} was provided) or it
#' will be constructed to be uniform on the log scale (if \code{dfRange)}
#' and \code{rangeLambda} were provided). Note that, if the EDF method is used,
#' the resulting \code{lambdaGrid} will vary for each value in \code{thetaGrid}.
#'
#' @return A \code{list} of Cholesky factors, parameter configurations, and
#' other information. The parameters are ordered such that \eqn{\lambda} runs
#' faster than \eqn{\theta}.
prepare <- function(S,
                    thetaGrid,
                    nu,
                    lambdaGrid = NULL,
                    dfRange = NULL,
                    nLambda = NULL,
                    rangeLambda = NULL) {

  nTheta <- length(thetaGrid)

  ## Create a matrix of size (nLambda, nTheta) of log(lambda) values:
  if (!is.null(lambdaGrid)) {

    ## The user has supplied the values of lambda. These values will be the same
    ## for all values of theta.
    nLambda <- length(lambdaGrid)
    logLambdaMat <- log(matrix(rep(lambdaGrid, nTheta), ncol = nTheta))

  } else if (!is.null(dfRange)) {

    ## dfRange was supplied, so we want the EDF design: This requires nLambda
    if (is.null(nLambda))
      stop("Since dfRange is not NULL, nLambda must be given")

    logLambdaMat <- makeDFGrid(
      nLambda = nLambda, thetaGrid = thetaGrid, dfRange = dfRange, S = S, nu = nu
    )

  } else {
    ## At a minimum, the user needs to give the number of lambda values and the
    ## range the values should take:
    if (is.null(nLambda) || is.null(rangeLambda))
      stop("You did not specify lambdaGrid, so nLambda and rangeLambda must be given")

    # Define log-scale-uniform lambdaGrid
    logLambdaGrid <- seq(log(rangeLambda[1]), log(rangeLambda[2]), length.out = nLambda)
    logLambdaMat  <- matrix(rep(logLambdaGrid, nTheta), ncol = nTheta)
  }

  xi <- cbind(
    logLambda = c(logLambdaMat),
    logTheta  = rep(log(thetaGrid), each = nLambda),
    logNu     = rep(log(nu), nTheta * nLambda)
  )

  chols <- compute_Cholesky_Factors(nu = nu, thetaGrid = thetaGrid, S = S)

  return(list(chols = chols, xi = xi))
}


#' Wrapper for \code{\link{prepare}()} that facilitates varying the smoothness parameter.
#'
#' @inheritParams prepare
#' @param N Integer indicating the side length of the assumed-square
#' spatial domain. This argument is ignored if \code{S} is not \code{NULL}.
#' @param S matrix of spatial locations (default \code{NULL})
#' @param path Save path.
#' @param nNu,rangeNu,nuGrid The number, range of values, or complete grid of
#' grid of values for \eqn{\nu}. Some can be left \code{NULL}; see details.
#' @param nTheta,rangeTheta,thetaGrid The number, range of values, or complete
#' grid of values for \eqn{\theta}. Some can be left \code{NULL}; see details.
#' @param ... Arguments to pass on to \code{\link{prepare}()} involving the
#' parameter \eqn{\lambda}.
#' @param completeGrid Logical: If \code{TRUE}, all combinations of \eqn{\theta}
#' and \eqn{\nu} will be used; otherwise, only element-wise pairs of \eqn{\theta}
#' and \eqn{\nu} are used and \eqn{n_\theta} must equal \eqn{n_\nu}.
#' @param completeLambdaGrid Logical: If \code{TRUE}, supplied values of
#' \eqn{\lambda} will be used with each \eqn{\theta}-\eqn{\nu} pair. This is
#' \code{TRUE} by default, as the Cholesky factors are
#' independent of \eqn{\lambda}, and so it is computationally inexpensive to
#' generate parameter configurations for a large number of \eqn{\lambda} values.
#' Only applicable if \code{lambdaGrid = NULL} and
#' \code{completeGrid = FALSE}.
#'
#' @details If a grid of parameters is not provided (i.e., \code{xGrid} is
#' \code{NULL}), then, using \code{rangeX} and \code{nX}, a grid of parameters
#' is generated that is uniform on the log scale.
#'
#' @return A \code{list} of Cholesky factors, parameter configurations, and
#' other information. The parameters are ordered such that \eqn{\lambda} runs
#' faster than
#' \eqn{\theta}, which in turn runs faster than \eqn{\nu}.
prepare_wrapper <- function(N = 16, S = NULL,
                            path = NULL,
                            completeGrid = TRUE, completeLambdaGrid = TRUE,
                            nNu = NULL, rangeNu = NULL, nuGrid = NULL,
                            nTheta = NULL, rangeTheta = NULL, thetaGrid = NULL,
                            saveLambda = TRUE,
                            ...) {

  if (is.null(S)) {
    #cat("The spatial domain is [0, N] x [0, N], where N =", N,"\b.\n")
    S <- make.surface.grid(list(x = 1:N, y = 1:N))
  }

  ## Matrix of spatial distances between each location in S.
  ## rdist() is used within stationary.cov(), so we also use it for consistency.
  D <- rdist(S)

  if (is.null(path)) cat("No path provided: Computed objects will not be saved.")

  if (is.null(nuGrid)) {
    if (is.null(nNu) || is.null(rangeNu))
      stop("nuGrid is NULL, so nNu and rangeNu must be given")

    # cat("Defining log-uniform nuGrid based on nNu and rangeNu.\n")
    logNuGrid <- seq(log(rangeNu[1]), log(rangeNu[2]), length.out = nNu)
    nuGrid <- exp(logNuGrid)
  } else {
    nNu <- length(nuGrid)
    logNuGrid <- log(nuGrid)
  }

  if (is.null(thetaGrid)) {
    if (is.null(nTheta) || is.null(rangeTheta))
      stop("thetaGrid is NULL, so nTheta and rangeTheta must be given")

    # cat("Defining log-uniform thetaGrid based on nTheta and rangeTheta.\n")
    logThetaGrid <- seq(log(rangeTheta[1]), log(rangeTheta[2]), length.out = nTheta)
    thetaGrid <- exp(logThetaGrid)
  } else {
    nTheta <- length(thetaGrid)
  }

  if (completeGrid) {

    # cat("Since completeGrid = TRUE, all combinations of theta and nu will be used.\n")
    nChols <- nTheta * nNu
    cat("Computing", nChols, "Cholesky factors...\n")

    time <- system.time(
      param_configs <- lapply(
        nuGrid,
        function(nu) prepare(thetaGrid = thetaGrid, nu = nu, S = S, ...)
      )
    )

  } else if (!completeGrid) {

    if (nTheta != nNu) {
      stop("Since completeGrid = FALSE, we are not generating a complete-grid of parameters and length(thetaGrid) must equal length(nuGrid).")
    } else {
      nChols <- nTheta # also equal to nNu
    }

    # cat("Since completeGrid = FALSE, we will generate parameter configurations for the pairs (thetaGrid[i], nuGrid[i]), i = 1, ...\n")
    cat(paste("Computing", nChols, "Cholesky factors... "))

    if (completeLambdaGrid) {
      # Apply prepare() for a single theta value at a time; entire lambdaGrid
      # passed into prepare() at each iteration
      time <- system.time(
        param_configs <- lapply(
          1:nChols,
          function(i) prepare(thetaGrid = thetaGrid[i], nu = nuGrid[i], S = S, ...)
        )
      )

    } else if (!completeLambdaGrid) {
      # Apply prepare() for a single theta, lambda value at a time

      lambdaGrid <- list(...)$lambdaGrid

      if (length(lambdaGrid) != nChols) stop("Length of lambdaGrid should equal the length of thetaGrid and nuGrid")

      time <- system.time(
        param_configs <- lapply(
          1:nChols,
          function(i) do.call(prepare, list(thetaGrid = thetaGrid[i], lambdaGrid = lambdaGrid[i], nu = nuGrid[i], S = S))
        )
      )
    }
  }

  cat(paste("Cholesky factors computed in", round(time[["elapsed"]], 4), "seconds, with each factor taking", round(time[["elapsed"]]/nChols, 4), "seconds.\n"))

  # Extract and combine the cholesky factors
  chols <- lapply(param_configs, function(x) x$chols)
  chols <- abind(chols, along = 3)

  # extract and combine the true parameter values
  xi <- lapply(param_configs, function(x) x$xi)
  xi <- abind(xi, along = 1)

  # Convert the parameters from the original scale to the log-scale, and
  # convert lambda to the square-root of lambda (the standard deviation)
  xi <- exp(xi)
  xi[, 1] <- sqrt(xi[, 1])

  colnames(xi) <- gsub("log", "", colnames(xi))
  colnames(xi) <- gsub("Lambda", "sigma", colnames(xi))
  colnames(xi) <- gsub("Theta", "rho", colnames(xi))
  colnames(xi) <- gsub("Nu", "nu", colnames(xi))

  if (!saveLambda) xi <- xi[, -1]

  ## Save objects
  if (!is.null(path)) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    path2 <- paste0(path, "dummy") # just in case path ends in a "/"
    save(file = paste(dirname(path2), "S.rda", sep = "/"), S)
    save(file = paste(dirname(path2), "D.rda", sep = "/"), D)
    save(file = paste0(path, "_chols", ".rda"), chols)
    save(file = paste0(path, "_xi", ".rda"), xi)
    save(file = paste0(path, "_logNuGrid", ".rda"), logNuGrid)
  }

  return(list(chols = chols, xi = xi, nNu = nNu))
}

## helper function to account for the effective degrees of freedom (EDF) when
## creating parameter configurations
# For
# example, if \code{dfRange = c(1, 256)}, then
# \code{lambdaGrid} is constructed such that the model resulting from
# \eqn{(\lambda_i, \theta, \nu)}, \eqn{i = 1, ..., n_\lambda}, have
# an EDF 1 and 256. Note that the resulting \code{lambdaGrid}
# will vary for each value in \code{thetaGrid}.
makeDFGrid <- function(S, nLambda, thetaGrid, dfRange, nu) {

  ## This function works by defining a range of EDF that we want to use: For
  ## example, dfRange = c(1, 256) means that we want to find nLambda
  ## values of lambda such that (lambda_i, theta), i = 1, ..., nLambda have an
  ## EDF of between 1 and 256. It works by defining a large range of possible
  ## lambda values (300 values between 1e-6 and 2000), computing the EDF for each
  ## of those values, and then interpolating to get the values of lambda that
  ## approximately yield the desired EDF.

  nL     <- 300
  lRange <- c(1e-6, 2000) # range of lambda values that will be used to try to match EDF with
  lGrid  <- seq(log(lRange[1]), log(lRange[2]), length.out = nL)

  nTheta <- length(thetaGrid)
  df     <- matrix(NA, nL, nTheta)
  for (j in 1:nTheta) {
    df[, j] <- findDF(S, theta = thetaGrid[j], lGrid, nu)
  }

  ## The values for EDF that we want to use
  dfGrid <- seq(dfRange[1], dfRange[2], length.out = nLambda)

  ## For each value of theta, find the value of lambda associated with each EDF.
  ## Interpolation between EDF values is done using splint().
  logLambdaMat <- matrix(NA, nLambda, nTheta)
  for (j in 1:nTheta) {
    logLambdaMat[, j] <- splint(log(df[, j]), lGrid, log(dfGrid))
  }

  return(logLambdaMat)
}


findDF <- function(S, theta, lGrid, nu) {

  sigma <- stationary.cov(S, S, Covariance = "Matern", theta = theta, smoothness = nu)
  d <- eigen(sigma, symmetric = TRUE)$values

  nLambda <- length(lGrid)
  alphaAll <- exp(lGrid) / (1 + exp(lGrid))

  out <- rep(NA, nLambda)
  for (k in 1:nLambda) {
    alpha <- alphaAll[k]
    out[k] <- sum((1 - alpha) * d / ((1 - alpha) * d + alpha))
  }
  return(out)
}


# ---- Parameter sampling ----

# Latin Hypercube Sampling
lhs_design <- function(n_samples, range_params) {

  n_params <- ncol(range_params)
  param_names <- colnames(range_params)
  A <- lhs::randomLHS(n_samples, n_params)

  ## Make a list of transforms. The parameters will be simulated to be uniform
  ## on the original scale if their range contains negative values, and uniform
  ## on the log scale otherwise (using the log-uniform distribution, also
  ## known as the reciprocal distribution).
  neg_params <- which(apply(range_params, 2, function(x) any(x < 0)))
  transforms <- lapply(1:n_params, function(i) {
    if (i %in% neg_params) {
      function(x, range) qunif(x, min = range[1], max = range[2])
    } else {
      function(x, range) KScorrect::qlunif(x, min = range[1], max = range[2])
    }
  })

  B <- sapply(1:n_params, function(i) transforms[[i]](A[, i], range_params[, i]))

  colnames(B) <- param_names

  return(B)
}


sample_uniformly <- function(n_samples, range_params) {

  n_params <- ncol(range_params)
  param_names <- colnames(range_params)

  ## Make a list of transforms. The parameters will be simulated to be uniform
  ## on the original scale if their range contains negative values, and uniform
  ## on the log scale otherwise (using the log-uniform distribution, also
  ## known as the reciprocal distribution).
  neg_params <- which(apply(range_params, 2, function(x) any(x < 0)))
  transforms <- lapply(1:n_params, function(i) {
    if (i %in% neg_params) {
      function(x, range) runif(x, min = range[1], max = range[2])
    } else {
      function(x, range) KScorrect::rlunif(x, min = range[1], max = range[2])
    }
  })


  B <- sapply(1:n_params, function(i) transforms[[i]](n_samples, range_params[, i]))

  colnames(B) <- param_names

  return(B)
}
