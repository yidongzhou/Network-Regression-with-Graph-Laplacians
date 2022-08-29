#' @title Global Network Regression
#' @description Global regression for network data quantified by graph Laplacians with Euclidean predictors.
#' @param gl a list (length n) of m by m matrices (resp. an m by m by n array) where
#'   \code{M[[i]]} (resp. \code{M[, , i]}) holds the i-th graph Laplacian
#'   of dimension m by m.
#' @param x an n by p matrix or data frame of predictors.
#'   It can be a vector of length n if p = 1.
#' @param xOut an nOut by p matrix or data frame of output predictor levels.
#'   It can be a vector of length p if m = 1 or a vector of length nOut if p = 1
#'   or a scalar if both p and nOut are equal to 0. Default is \code{NULL}.
#' @param optns a list of options control parameters specified by
#'   \code{list(name = value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{metric}{choice of metric. 'frobenius' and 'power' are supported,
#'   which corresponds to Frobenius metric and Euclidean power metric,
#'   respectively. One can use 'f' and 'p' for short.
#'   Default is Frobenius metric.}
#' \item{alpha}{the power for Euclidean power metric.
#'   Default is 1 which corresponds to Frobenius metric.}
#' \item{digits}{the integer indicating the number of decimal places (round)
#'   to be kept in the output. Default is \code{NULL}, which means no round operation.}
#' }
#' @return A \code{nr} object --- a list containing the following fields:
#' \item{fit}{a list of estimated graph Laplacians at \code{x}.}
#' \item{predict}{a list of predicted graph Laplacians at \code{xOut}.
#'   Included if \code{xOut} is not \code{NULL}.}
#' \item{RSquare}{Fréchet coefficient of determination.}
#' \item{AdjRSquare}{adjusted Fréchet coefficient of determination.}
#' \item{residuals}{Frobenius distance between the true and
#'   fitted graph Laplacians.}
#' \item{gl}{the graph Laplacians used.}
#' \item{x}{the predictor used.}
#' \item{xOut}{the output predictor level used.}
#' \item{optns}{the control options used.}
#' @examples
#' set.seed(1)
#' n <- 100
#' m <- 10
#' d <- m * (m - 1) / 2
#' xOut <- seq(0.1, 0.9, length.out = 9)
#' x <- runif(n, min = 0, max = 1)
#' gl <- list()
#' for (i in 1:n) {
#'   glVec <- -rbeta(d, shape1 = x[i], shape2 = 1 - x[i])
#'   temp <- matrix(0, nrow = m, ncol = m)
#'   temp[lower.tri(temp)] <- glVec
#'   temp <- temp + t(temp)
#'   diag(temp) <- -colSums(temp)
#'   gl[[i]] <- temp
#' }
#' res <- gnr(gl, x, xOut)
#' res$predict
#' @references
#' \itemize{
#' \item \cite{Zhou, Y. and Müller, H.G., 2022. Network Regression with Graph Laplacians. arXiv preprint arXiv:2109.02981.}
#' \item \cite{Petersen, A. and Müller, H.-G. (2019). Fréchet regression for random objects with Euclidean predictors. The Annals of Statistics, 47(2), 691--719.}
#' }
#' @export

gnr <- function(gl = NULL, x = NULL, xOut = NULL, optns = list()) {
  if (is.null(gl) | is.null(x)) {
    stop("requires the input of both gl and x")
  }
  if (is.null(optns$metric)) {
    optns$metric <- "frobenius"
  }
  if (!(optns$metric %in% c("frobenius", "power"))) {
    stop("metric choice not supported")
  }
  if (is.null(optns$alpha)) {
    optns$alpha <- 1
  }
  if (optns$alpha < 0) {
    stop("alpha must be non-negative")
  }
  if (is.null(optns$digits)) {
    optns$digits <- NA
  }
  if (!is.matrix(x)) {
    if (is.data.frame(x) | is.vector(x)) {
      x <- as.matrix(x)
    } else {
      stop("x must be a matrix or a data frame or a vector")
    }
  }
  n <- nrow(x) # number of observations
  p <- ncol(x) # number of covariates
  if (!is.list(gl)) {
    if (is.array(gl)) {
      gl <- lapply(seq(dim(gl)[3]), function(i) gl[, , i])
    } else {
      stop("gl must be a list or an array")
    }
  }
  if (length(gl) != n) {
    stop("the number of rows in x must be the same as the number of graph Laplacians in gl")
  }
  if (!is.null(xOut)) {
    if (!is.matrix(xOut)) {
      if (is.data.frame(xOut)) {
        xOut <- as.matrix(xOut)
      } else if (is.vector(xOut)) {
        if (p == 1) {
          xOut <- as.matrix(xOut)
        } else {
          xOut <- t(xOut)
        }
      } else {
        stop("xOut must be a matrix or a data frame or a vector")
      }
    }
    if (ncol(xOut) != p) {
      stop("x and xOut must have the same number of columns")
    }
    nOut <- nrow(xOut) # number of predictions
  } else {
    nOut <- 0
  }
  nodes <- colnames(gl[[1]])
  m <- ncol(gl[[1]]) # number of nodes
  if (is.null(nodes)) nodes <- 1:m
  glVec <- matrix(unlist(gl), ncol = m^2, byrow = TRUE) # n by m^2
  if (substr(optns$metric, 1, 1) == "p") {
    glAlpha <- lapply(gl, function(gli) {
      eigenDecom <- eigen(gli)
      Lambda <- pmax(Re(eigenDecom$values), 0) # exclude 0i
      U <- eigenDecom$vectors
      U %*% diag(Lambda^optns$alpha) %*% t(U)
    })
    glAlphaVec <- matrix(unlist(glAlpha), ncol = m^2, byrow = TRUE) # n by m^2
  }

  # initialization of OSQP solver
  W <- 2^32 # bound on the weights of edges in the graph
  nConsts <- m^2 # number of constraints
  l <- c(rep.int(0, m * (m + 1) / 2), rep.int(-W, m * (m - 1) / 2))
  u <- rep.int(0, nConsts)
  q <- rep.int(0, m^2)
  P <- diag(m^2)
  consts <- matrix(0, nrow = nConsts, ncol = m^2)
  k <- 0
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      k <- k + 1
      consts[k, (j - 1) * m + i] <- 1
      consts[k, (i - 1) * m + j] <- -1
    }
  }
  for (i in 1:m) {
    consts[k + i, ((i - 1) * m + 1):(i * m)] <- rep(1, m)
  }
  k <- k + m
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      k <- k + 1
      consts[k, (j - 1) * m + i] <- 1
    }
  }
  model <- osqp::osqp(P, q, consts, l, u, osqp::osqpSettings(verbose = FALSE))

  xMean <- colMeans(x)
  invVa <- solve(var(x) * (n - 1) / n)
  fit <- vector(mode = "list", length = n)
  residuals <- rep.int(0, n)
  wc <- t(apply(x, 1, function(xi) t(xi - xMean) %*% invVa)) # n by p
  totVa <- sum((scale(glVec, scale = FALSE))^2)
  if (nrow(wc) != n) wc <- t(wc) # for p=1
  if (substr(optns$metric, 1, 1) == "f") {
    for (i in 1:n) {
      w <- apply(wc, 1, function(wci) 1 + t(wci) %*% (x[i, ] - xMean))
      qNew <- apply(glVec, 2, weighted.mean, w) # m^2
      model$Update(q = -qNew)
      temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
      temp <- (temp + t(temp)) / 2 # symmetrize
      if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
      temp[temp > 0] <- 0 # off diagonal should be negative
      diag(temp) <- 0
      diag(temp) <- -colSums(temp)
      fit[[i]] <- temp
      residuals[i] <- sqrt(sum((gl[[i]] - temp)^2))
    }
    resVa <- sum(residuals^2)
    RSquare <- 1 - resVa / totVa
    AdjRSquare <- RSquare - (1 - RSquare) * p / (n - p - 1)
    if (nOut > 0) {
      predict <- vector(mode = "list", length = nOut)
      for (i in 1:nOut) {
        w <- apply(wc, 1, function(wci) 1 + t(wci) %*% (xOut[i, ] - xMean))
        qNew <- apply(glVec, 2, weighted.mean, w) # m^2
        model$Update(q = -qNew)
        temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
        temp <- (temp + t(temp)) / 2 # symmetrize
        if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
        temp[temp > 0] <- 0 # off diagonal should be negative
        diag(temp) <- 0
        diag(temp) <- -colSums(temp)
        predict[[i]] <- temp
      }
      res <- list(fit = fit, predict = predict, RSquare = RSquare, AdjRSquare = AdjRSquare, residuals = residuals, gl = gl, x = x, xOut = xOut, optns = optns)
    } else {
      res <- list(fit = fit, RSquare = RSquare, AdjRSquare = AdjRSquare, residuals = residuals, gl = gl, x = x, optns = optns)
    }
  } else if (substr(optns$metric, 1, 1) == "p") {
    # bAlpha <- matrix(colMeans(glAlphaVec), ncol = m)# m by m
    # eigenDecom <- eigen(bAlpha)
    # Lambda <- pmax(Re(eigenDecom$values), 0)# projection to M_m
    # U <- eigenDecom$vectors
    # qNew <- as.vector(U%*%diag(Lambda^(1/optns$alpha))%*%t(U))# inverse power
    # model$Update(q = -qNew)
    # omegaPlus <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
    # eigenDecom <- eigen(omegaPlus)
    # Lambda <- pmax(Re(eigenDecom$values), 0)
    # U <- eigenDecom$vectors
    # omegaPlusAlpha <- U%*%diag(Lambda^optns$alpha)%*%t(U)
    # totVa <- sum((t(glAlphaVec)-as.vector(omegaPlusAlpha))^2)# using Euclidean power metric
    for (i in 1:n) {
      w <- apply(wc, 1, function(wci) 1 + t(wci) %*% (x[i, ] - xMean))
      bAlpha <- matrix(apply(glAlphaVec, 2, weighted.mean, w), ncol = m) # m by m
      eigenDecom <- eigen(bAlpha)
      Lambda <- pmax(Re(eigenDecom$values), 0) # projection to M_m
      U <- eigenDecom$vectors
      qNew <- as.vector(U %*% diag(Lambda^(1 / optns$alpha)) %*% t(U)) # inverse power
      model$Update(q = -qNew)
      temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
      temp <- (temp + t(temp)) / 2 # symmetrize
      if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
      temp[temp > 0] <- 0 # off diagonal should be negative
      diag(temp) <- 0
      diag(temp) <- -colSums(temp)
      fit[[i]] <- temp
      residuals[i] <- sqrt(sum((gl[[i]] - temp)^2))
      # eigenDecom <- eigen(fit[[i]])
      # Lambda <- pmax(Re(eigenDecom$values), 0)
      # U <- eigenDecom$vectors
      # fitiAlpha <- U%*%diag(Lambda^optns$alpha)%*%t(U)
      # residuals[i] <- sqrt(sum((glAlpha[[i]]-fitiAlpha)^2))# using Euclidean power metric
    }
    resVa <- sum(residuals^2)
    RSquare <- 1 - resVa / totVa
    AdjRSquare <- RSquare - (1 - RSquare) * p / (n - p - 1)
    if (nOut > 0) {
      predict <- vector(mode = "list", length = nOut)
      for (i in 1:nOut) {
        w <- apply(wc, 1, function(wci) 1 + t(wci) %*% (xOut[i, ] - xMean))
        bAlpha <- matrix(apply(glAlphaVec, 2, weighted.mean, w), ncol = m) # m^2
        eigenDecom <- eigen(bAlpha)
        Lambda <- pmax(Re(eigenDecom$values), 0) # projection to M_m
        U <- eigenDecom$vectors
        qNew <- as.vector(U %*% diag(Lambda^(1 / optns$alpha)) %*% t(U)) # inverse power
        model$Update(q = -qNew)
        temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
        temp <- (temp + t(temp)) / 2 # symmetrize
        if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
        temp[temp > 0] <- 0 # off diagonal should be negative
        diag(temp) <- 0
        diag(temp) <- -colSums(temp)
        predict[[i]] <- temp
      }
      res <- list(fit = fit, predict = predict, RSquare = RSquare, AdjRSquare = AdjRSquare, residuals = residuals, gl = gl, x = x, xOut = xOut, optns = optns)
    } else {
      res <- list(fit = fit, RSquare = RSquare, AdjRSquare = AdjRSquare, residuals = residuals, gl = gl, x = x, optns = optns)
    }
  }
  class(res) <- "nr"
  res
}
