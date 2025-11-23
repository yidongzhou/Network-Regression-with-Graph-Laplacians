########################################################################
# Simulation for the paper "Network Regression with Graph Laplacians". #
########################################################################
source("src/gnr.R")
source("src/lnr.R")
source("src/kerFctn.R")
source("severn/functions_needed.R")
Q <- 1000
nVec <- c(50, 100, 200, 500, 1000)
N <- 9
m <- c(5, 5)
theta <- c(0.5, 0.2, 0.5)
d <- c(m[1]*(m[1] - 1)/2, m[1]*m[2], m[2]*(m[2] - 1)/2)
alpha <- 0.5
# initialization of OSQP solver
W <- 2^32 # bound on the weights of edges in the graph
nConsts <- sum(m)^2 # number of constraints
l <- c(rep.int(0, sum(m) * (sum(m) + 1) / 2), rep.int(-W, sum(m) * (sum(m) - 1) / 2))
u <- rep.int(0, nConsts)
q <- rep.int(0, sum(m)^2)
P <- diag(sum(m)^2)
consts <- matrix(0, nrow = nConsts, ncol = sum(m)^2)
k <- 0
for (i in 1:(sum(m) - 1)) {
  for (j in (i + 1):sum(m)) {
    k <- k + 1
    consts[k, (j - 1) * sum(m) + i] <- 1
    consts[k, (i - 1) * sum(m) + j] <- -1
  }
}
for (i in 1:sum(m)) {
  consts[k + i, ((i - 1) * sum(m) + 1):(i * sum(m))] <- rep(1, sum(m))
}
k <- k + sum(m)
for (i in 1:(sum(m) - 1)) {
  for (j in (i + 1):sum(m)) {
    k <- k + 1
    consts[k, (j - 1) * sum(m) + i] <- 1
  }
}
model <- osqp::osqp(P, q, consts, l, u, osqp::osqpSettings(verbose = FALSE))

# Helper function to create OSQP model (needed in parallel workers)
create_osqp_model <- function(m_size) {
  W <- 2^32
  nConsts <- m_size^2
  l <- c(rep.int(0, m_size * (m_size + 1) / 2), rep.int(-W, m_size * (m_size - 1) / 2))
  u <- rep.int(0, nConsts)
  q <- rep.int(0, m_size^2)
  P <- diag(m_size^2)
  consts <- matrix(0, nrow = nConsts, ncol = m_size^2)
  k <- 0
  for (i in 1:(m_size - 1)) {
    for (j in (i + 1):m_size) {
      k <- k + 1
      consts[k, (j - 1) * m_size + i] <- 1
      consts[k, (i - 1) * m_size + j] <- -1
    }
  }
  for (i in 1:m_size) {
    consts[k + i, ((i - 1) * m_size + 1):(i * m_size)] <- rep(1, m_size)
  }
  k <- k + m_size
  for (i in 1:(m_size - 1)) {
    for (j in (i + 1):m_size) {
      k <- k + 1
      consts[k, (j - 1) * m_size + i] <- 1
    }
  }
  osqp::osqp(P, q, consts, l, u, osqp::osqpSettings(verbose = FALSE))
}

# Parallel computing
library(doParallel)
NUM_CORES <- 36
cl <- makeCluster(NUM_CORES)
# Export required packages and source files to workers
clusterEvalQ(cl, {
  source("src/gnr.R")
  source("src/lnr.R")
  source("src/kerFctn.R")
  source("severn/functions_needed.R")
})
# Export the helper function and variables needed in parallel blocks
# Note: Variables like xSeq, LMeanRef4, etc. will be exported when they're created
clusterExport(cl, c("create_osqp_model", "m", "alpha", "d", "theta", "N", "Q", "nVec"))
registerDoParallel(cl)

#################################################################
# Part I: verification of the theoretical rates of convergence. #
#################################################################
# Simulation scenario I
set.seed(1)
ise1 <- foreach(n = nVec, .combine = 'cbind') %:%
  foreach(icount(Q), .combine = 'c') %dopar% {
    X <- runif(n, min = 0, max = 1)
    L <- list()
    LMean <- list()
    for(i in 1:n){
      LVec <- -rbeta(sum(d), shape1 = X[i], shape2 = 1-X[i])
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      L[[i]] <- temp

      LMeanVec <- -rep(X[i], sum(d))
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LMeanVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      LMean[[i]] <- temp
    }
    sim <- gnr(L, X)
    mean(sapply(1:n, function(i) sum((LMean[[i]]-sim$fit[[i]])^2)))
  }

mise1 <- colMeans(ise1)
fit1 <- lm(log(mise1)~log(nVec))
summary(fit1)
plot(log(mise1), log(nVec))

# Simulation scenario II
set.seed(1)
ise2 <- foreach(n = nVec, .combine = 'cbind') %:%
  foreach(icount(Q), .combine = 'c') %dopar% {
    X <- runif(n, min = 0, max = 1)
    L <- list()
    LMean <- list()
    for(i in 1:n){
      LVec <- -rbeta(sum(d), shape1 = X[i], shape2 = 1 - X[i])
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      eigenDecom <- eigen(temp)
      Lambda <- eigenDecom$values
      U <- eigenDecom$vectors
      L[[i]] <- U%*%diag(Lambda^(1/alpha))%*%t(U)
      
      LMeanVec <- -rep(X[i], sum(d))
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LMeanVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      eigenDecom <- eigen(temp)
      Lambda <- pmax(eigenDecom$values, 0)
      U <- eigenDecom$vectors
      temp <- U%*%diag(Lambda^(1/alpha))%*%t(U)
      model_local <- create_osqp_model(sum(m))
      model_local$Update(q = -temp)
      temp <- matrix(model_local$Solve()$x, ncol = sum(m))
      temp <- (temp + t(temp)) / 2 # symmetrize
      temp[temp > 0] <- 0 # off diagonal should be negative
      diag(temp) <- 0
      diag(temp) <- -colSums(temp)
      LMean[[i]] <- temp
    }
    sim <- gnr(L, X, optns = list(metric = "power", alpha = alpha))
    mean(sapply(1:n, function(i) sum((LMean[[i]]-sim$fit[[i]])^2)))
  }

mise2 <- colMeans(ise2)
fit2 <- lm(log(mise2)~log(nVec))
summary(fit2)
plot(log(mise2), log(nVec))

# Simulation scenario III
set.seed(1)
bw3 <- c(0.088, 0.075, 0.065, 0.054, 0.047)
ise3 <- foreach(n = nVec, .combine = 'cbind') %:%
  foreach(icount(Q), .combine = 'c') %dopar% {
    X <- runif(n, min = 0, max = 1)
    L <- list()
    LMean <- list()
    for(i in 1:n){
      LVec <- -rbeta(sum(d), shape1 = sin(pi * X[i]), shape2 = 1 - sin(pi * X[i]))
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      L[[i]] <- temp
      
      LMeanVec <- -rep(sin(pi * X[i]), sum(d))
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LMeanVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      LMean[[i]] <- temp
    }
    sim <- lnr(L, X, optns = list(bw = bw3[which(nVec == n)]))
    mean(sapply(1:n, function(i) sum((LMean[[i]]-sim$fit[[i]])^2)))
  }

mise3 <- colMeans(ise3)
fit3 <- lm(log(mise3)~log(nVec))
summary(fit3)
plot(log(mise3), log(nVec))

# Simulation scenario IV
set.seed(1)
xSeq <- seq(0, 1, length.out = 1001)
LMeanRef4 <- list()
for(i in 1:1001) {
  x <- xSeq[i]
  y <- matrix(0, nrow = sum(m), ncol = sum(m))
  for(j in 1:10000){
    LMeanVec <- -rbeta(sum(d), shape1 = sin(pi * x), shape2 = 1 - sin(pi * x))
    temp <- matrix(0, nrow = sum(m), ncol = sum(m))
    temp[lower.tri(temp)] <- LMeanVec
    temp <- temp + t(temp)
    diag(temp) <- -colSums(temp)
    eigenDecom <- eigen(temp)
    Lambda <- pmax(eigenDecom$values, 0)
    U <- eigenDecom$vectors
    temp <- U%*%diag(Lambda^alpha)%*%t(U)
    y <- y + temp/10000
  }
  eigenDecom <- eigen(y)
  Lambda <- pmax(eigenDecom$values, 0)
  U <- eigenDecom$vectors
  y <- U%*%diag(Lambda^(1/alpha))%*%t(U)

  model$Update(q = -y)
  y <- matrix(model$Solve()$x, ncol = sum(m))
  y <- (y + t(y)) / 2 # symmetrize
  y[y > 0] <- 0 # off diagonal should be negative
  diag(y) <- 0
  diag(y) <- -colSums(y)
  LMeanRef4[[i]] <- y
}
# Export LMeanRef4 and xSeq to parallel workers
clusterExport(cl, c("LMeanRef4", "xSeq"))

set.seed(1)
bw4 <- c(0.051, 0.044, 0.038, 0.032, 0.027)
ise4 <- foreach(n = nVec, .combine = 'cbind') %:%
  foreach(icount(Q), .combine = 'c') %dopar% {
    X <- runif(n, min = 0, max = 1)
    L <- list()
    LMean <- list()
    for(i in 1:n){
      LVec <- -rbeta(sum(d), shape1 = sin(pi * X[i]), shape2 = 1 - sin(pi * X[i]))
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      L[[i]] <- temp
      
      idx <- which.min(abs(X[i] - xSeq))
      LMean[[i]] <- LMeanRef4[[idx]]
    }
    sim <- lnr(L, X, optns = list(metric = 'power', alpha = alpha, bw = bw4[which(nVec == n)]))
    mean(sapply(1:n, function(i) sum((LMean[[i]] - sim$fit[[i]])^2)))
  }

mise4 <- colMeans(ise4)
fit4 <- lm(log(mise4)~log(nVec))
summary(fit4)
plot(log(mise4), log(nVec))

##################################################
# Part II: Networks with latent block structure. #
##################################################
# Simulation scenario I
set.seed(1)
isesbm1 <- foreach(n = nVec, .combine = cbind) %:%
  foreach(icount(Q), .combine = c) %dopar% {
    X <- runif(n, min = 0, max = 1)
    L <- list()
    LMean <- list()
    for(i in 1:n){
      bk1Vec <- -rbinom(d[1], 1, theta[1])*rbeta(d[1], shape1 = X[i], shape2 = 1-X[i])
      bk2Vec <- -rbinom(d[2], 1, theta[2])*rbeta(d[2], shape1 = X[i], shape2 = 1-X[i])
      bk3Vec <- -rbinom(d[3], 1, theta[3])*rbeta(d[3], shape1 = X[i], shape2 = 1-X[i])
      temp1 <- matrix(0, nrow = m[1], ncol = m[1])
      temp1[lower.tri(temp1)] <- bk1Vec
      temp1 <- temp1 + t(temp1)
      temp2 <- matrix(bk2Vec, nrow = m[1])
      temp3 <- matrix(0, nrow = m[2], ncol = m[2])
      temp3[lower.tri(temp3)] <- bk3Vec
      temp3 <- temp3 + t(temp3)
      temp <- rbind(cbind(temp1, temp2), cbind(t(temp2), temp3))
      diag(temp) <- -colSums(temp)
      L[[i]] <- temp

      temp1 <- matrix(0, nrow = m[1], ncol = m[1])
      temp1[lower.tri(temp1)] <- -rep(theta[1]*X[i], d[1])
      temp1 <- temp1 + t(temp1)
      temp2 <- matrix(-rep(theta[2]*X[i], d[2]), nrow = m[1])
      temp3 <- matrix(0, nrow = m[2], ncol = m[2])
      temp3[lower.tri(temp3)] <- -rep(theta[3]*X[i], d[3])
      temp3 <- temp3 + t(temp3)
      temp <- rbind(cbind(temp1, temp2), cbind(t(temp2), temp3))
      diag(temp) <- -colSums(temp)
      LMean[[i]] <- temp
    }
    intercept <- matrix(0, sum(m), sum(m))
    slope <- matrix(0, sum(m), sum(m))
    for(j in 1:sum(m)){
      for(k in j:sum(m)){
        coefGl <- coef(lm(sapply(L, function(Li) Li[j, k]) ~ X))
        intercept[j, k] <- coefGl[1]
        slope[j, k] <- coefGl[2]
        intercept[k, j] <- intercept[j, k]
        slope[k, j] <- slope[j, k]
      }
    }
    fit <- list()
    for(i in 1:n){
      fit[[i]] <- as.matrix(proj_sparse(intercept + X[i]*slope))
    }
    errsd <- mean(sapply(1:n, function(i) sum((LMean[[i]]-fit[[i]])^2)))

    sim <- gnr(L, X)
    err <- mean(sapply(1:n, function(i) sum((LMean[[i]]-sim$fit[[i]])^2)))
    c(errsd, err)
  }
isesbm1sd <- isesbm1[seq(1, 2*Q, by = 2), ]
isesbm1 <- isesbm1[seq(2, 2*Q, by = 2), ]

misesbm1 <- colMeans(isesbm1)
fitsbm1 <- lm(log(misesbm1)~log(nVec))
summary(fitsbm1)
plot(log(misesbm1), log(nVec))

# misesbm1sd <- colMeans(isesbm1sd)
# misesbm1sd - misesbm1
# 50 100 200 500 1000 
# 0.65399121 0.32821751 0.16554847 0.06726054 0.03354503 Dryden
# 0.65399022 0.32821639 0.16554715 0.06725911 0.03354389

# Simulation scenario II
set.seed(1)
isesbm2L2 <- matrix(0, nrow = 2 * Q, ncol = length(nVec))
for(n in nVec) {
  for(q in 1:Q) {
    X <- runif(n, min = 0, max = 1)
    L <- list()
    LMean <- list()
    for(i in 1:n){
      bk1Vec <- -rbinom(d[1], 1, theta[1])*rbeta(d[1], shape1 = X[i], shape2 = 1 - X[i])
      bk2Vec <- -rbinom(d[2], 1, theta[2])*rbeta(d[2], shape1 = X[i], shape2 = 1 - X[i])
      bk3Vec <- -rbinom(d[3], 1, theta[3])*rbeta(d[3], shape1 = X[i], shape2 = 1 - X[i])
      temp1 <- matrix(0, nrow = m[1], ncol = m[1])
      temp1[lower.tri(temp1)] <- bk1Vec
      temp1 <- temp1 + t(temp1)
      temp2 <- matrix(bk2Vec, nrow = m[1])
      temp3 <- matrix(0, nrow = m[2], ncol = m[2])
      temp3[lower.tri(temp3)] <- bk3Vec
      temp3 <- temp3 + t(temp3)
      temp <- rbind(cbind(temp1, temp2), cbind(t(temp2), temp3))
      diag(temp) <- -colSums(temp)
      eigenDecom <- eigen(temp)
      # Lambda <- pmax(eigenDecom$values, 0)
      Lambda <- eigenDecom$values
      U <- eigenDecom$vectors
      L[[i]] <- U%*%diag(Lambda^(1/alpha))%*%t(U)
      # L[[i]] <- temp

      temp1 <- matrix(0, nrow = m[1], ncol = m[1])
      temp1[lower.tri(temp1)] <- -rep(theta[1] * X[i], d[1])
      temp1 <- temp1 + t(temp1)
      temp2 <- matrix(-rep(theta[2] * X[i], d[2]), nrow = m[1])
      temp3 <- matrix(0, nrow = m[2], ncol = m[2])
      temp3[lower.tri(temp3)] <- -rep(theta[3] * X[i], d[3])
      temp3 <- temp3 + t(temp3)
      temp <- rbind(cbind(temp1, temp2), cbind(t(temp2), temp3))
      diag(temp) <- -colSums(temp)
      eigenDecom <- eigen(temp)
      Lambda <- pmax(eigenDecom$values, 0)
      U <- eigenDecom$vectors
      temp <- U%*%diag(Lambda^(1/alpha))%*%t(U)
      model$Update(q = -temp)
      temp <- matrix(model$Solve()$x, ncol = sum(m))
      temp <- (temp + t(temp)) / 2 # symmetrize
      temp[temp > 0] <- 0 # off diagonal should be negative
      diag(temp) <- 0
      diag(temp) <- -colSums(temp)
      LMean[[i]] <- temp
    }
    LRoot <- list()
    for(i in 1:n){
      LRoot[[i]] <- shapes::rootmat(L[[i]])
    }
    intercept <- matrix(0, sum(m), sum(m))
    slope <- matrix(0, sum(m), sum(m))
    for(j in 1:sum(m)){
      for(k in j:sum(m)){
        coefGl <- coef(lm(sapply(LRoot, function(Li) Li[j, k]) ~ X))
        intercept[j, k] <- coefGl[1]
        slope[j, k] <- coefGl[2]
        intercept[k, j] <- intercept[j, k]
        slope[k, j] <- slope[j, k]
      }
    }
    fit <- list()
    for(i in 1:n){
      fiti <- G_2(intercept + X[i]*slope)
      fit[[i]] <- as.matrix(proj_sparse(fiti %*% t(fiti)))
    }
    errsd <- mean(sapply(1:n, function(i) sum((LMean[[i]]-fit[[i]])^2)))

    sim <- gnr(L, X, optns = list(metric = "power", alpha = alpha))
    err <- mean(sapply(1:n, function(i) sum((LMean[[i]]-sim$fit[[i]])^2)))
    isesbm2L2[c(2 * q - 1, 2 * q), which(nVec == n)] <- c(errsd, err)
  }
}
isesbm2sdL2 <- isesbm2L2[seq(1, 2*Q, by = 2), ]
isesbm2L2 <- isesbm2L2[seq(2, 2*Q, by = 2), ]

# misesbm2L2 <- colMeans(isesbm2L2)
# fitsbm2L2 <- lm(log(misesbm2L2)~log(nVec))
# summary(fitsbm2L2)
# plot(log(misesbm2L2), log(nVec))
# 
# misesbm2sdL2 <- colMeans(isesbm2sdL2)
# misesbm2sdL2 - misesbm2L2
# 50 100 200 500 1000
# 16.3425910  8.0799012  4.2254665  1.6996689  0.8663227 Dryden
# 16.3425109  8.0798024  4.2254242  1.6996553  0.8663193

# Simulation scenario III
set.seed(1)
bw3 <- c(0.153, 0.133, 0.116, 0.096, 0.084)# n^-0.2
bw3sd <- c(0.153, 0.133, 0.116, 0.096, 0.084)
kh <- function(x, h){
  exp(-x^2/(2*h^2))/(sqrt(2*pi)*h)
}
isesbm3 <- foreach(n = nVec, .combine = 'cbind') %:%
  foreach(icount(Q), .combine = 'c') %dopar% {
    X <- runif(n, min = 0, max = 1)
    L <- list()
    LMean <- list()
    for(i in 1:n){
      bk1Vec <- -rbinom(d[1], 1, theta[1])*rbeta(d[1], shape1 = sin(pi*X[i]), shape2 = 1-sin(pi*X[i]))
      bk2Vec <- -rbinom(d[2], 1, theta[2])*rbeta(d[2], shape1 = sin(pi*X[i]), shape2 = 1-sin(pi*X[i]))
      bk3Vec <- -rbinom(d[3], 1, theta[3])*rbeta(d[3], shape1 = sin(pi*X[i]), shape2 = 1-sin(pi*X[i]))
      temp1 <- matrix(0, nrow = m[1], ncol = m[1])
      temp1[lower.tri(temp1)] <- bk1Vec
      temp1 <- temp1 + t(temp1)
      temp2 <- matrix(bk2Vec, nrow = m[1])
      temp3 <- matrix(0, nrow = m[2], ncol = m[2])
      temp3[lower.tri(temp3)] <- bk3Vec
      temp3 <- temp3 + t(temp3)
      temp <- rbind(cbind(temp1, temp2), cbind(t(temp2), temp3))
      diag(temp) <- -colSums(temp)
      L[[i]] <- temp

      temp1 <- matrix(0, nrow = m[1], ncol = m[1])
      temp1[lower.tri(temp1)] <- -rep(theta[1]*sin(pi*X[i]), d[1])
      temp1 <- temp1 + t(temp1)
      temp2 <- matrix(-rep(theta[2]*sin(pi*X[i]), d[2]), nrow = m[1])
      temp3 <- matrix(0, nrow = m[2], ncol = m[2])
      temp3[lower.tri(temp3)] <- -rep(theta[3]*sin(pi*X[i]), d[3])
      temp3 <- temp3 + t(temp3)
      temp <- rbind(cbind(temp1, temp2), cbind(t(temp2), temp3))
      diag(temp) <- -colSums(temp)
      LMean[[i]] <- temp
    }
    glVec <- matrix(unlist(L), ncol = sum(m)^2, byrow = TRUE) # n by m^2
    fit <- list()
    for(i in 1:n){
      wi <- sapply(X[i] - X, kh, bw3sd[which(nVec == n)])
      temp <- matrix(apply(glVec, 2, weighted.mean, wi), ncol = sum(m))
      fit[[i]] <- as.matrix(proj_sparse(temp))
    }
    errsd <- mean(sapply(1:n, function(i) sum((LMean[[i]]-fit[[i]])^2)))

    sim <- lnr(L, X, optns = list(bw = bw3[which(nVec == n)]))
    err <- mean(sapply(1:n, function(i) sum((LMean[[i]]-sim$fit[[i]])^2)))
    c(errsd, err)
  }
isesbm3sd <- isesbm3[seq(1, 2*Q, by = 2), ]
isesbm3 <- isesbm3[seq(2, 2*Q, by = 2), ]

misesbm3 <- colMeans(isesbm3)
fitsbm3 <- lm(log(misesbm3)~log(nVec))
summary(fitsbm3)
plot(log(misesbm3), log(nVec))

# misesbm3sd <- colMeans(isesbm3sd)
# misesbm3sd - misesbm3
# 50 100 200 500 1000
# 2.8483122 1.8821186 1.2421967 0.7043861 0.4657804 Dryden
# 1.4495197 0.8269499 0.4832614 0.2314146 0.1361364

# Simulation scenario IV
set.seed(1)
xSeq <- seq(0, 1, length.out = 1001)
LMeanRefsbm4 <- list()
for(i in 1:1001) {
  x <- xSeq[i]
  y <- matrix(0, nrow = sum(m), ncol = sum(m))
  for(j in 1:10000){
    bk1Vec <- -rbinom(d[1], 1, theta[1])*rbeta(d[1], shape1 = sin(pi*x), shape2 = 1-sin(pi*x))
    bk2Vec <- -rbinom(d[2], 1, theta[2])*rbeta(d[2], shape1 = sin(pi*x), shape2 = 1-sin(pi*x))
    bk3Vec <- -rbinom(d[3], 1, theta[3])*rbeta(d[3], shape1 = sin(pi*x), shape2 = 1-sin(pi*x))
    temp1 <- matrix(0, nrow = m[1], ncol = m[1])
    temp1[lower.tri(temp1)] <- bk1Vec
    temp1 <- temp1 + t(temp1)
    temp2 <- matrix(bk2Vec, nrow = m[1])
    temp3 <- matrix(0, nrow = m[2], ncol = m[2])
    temp3[lower.tri(temp3)] <- bk3Vec
    temp3 <- temp3 + t(temp3)
    temp <- rbind(cbind(temp1, temp2), cbind(t(temp2), temp3))
    diag(temp) <- -colSums(temp)
    eigenDecom <- eigen(temp)
    Lambda <- pmax(eigenDecom$values, 0)
    U <- eigenDecom$vectors
    temp <- U%*%diag(Lambda^alpha)%*%t(U)
    y <- y + temp/10000
  }
  eigenDecom <- eigen(y)
  Lambda <- pmax(eigenDecom$values, 0)
  U <- eigenDecom$vectors
  y <- U%*%diag(Lambda^(1/alpha))%*%t(U)

  model$Update(q = -y)
  y <- matrix(model$Solve()$x, ncol = sum(m))
  y <- (y + t(y)) / 2 # symmetrize
  y[y > 0] <- 0 # off diagonal should be negative
  diag(y) <- 0
  diag(y) <- -colSums(y)
  LMeanRefsbm4[[i]] <- y
}
# Export LMeanRefsbm4 and xSeq to parallel workers
clusterExport(cl, c("LMeanRefsbm4", "xSeq"))

set.seed(1)
bw4 <- c(0.115, 0.1, 0.087, 0.073, 0.063)# n^-0.2
bw4sd <- c(0.115, 0.1, 0.087, 0.073, 0.063)
kh <- function(x, h){
  exp(-x^2/(2*h^2))/(sqrt(2*pi)*h)
}
# Export kh function to parallel workers
clusterExport(cl, "kh")
isesbm4 <- foreach(n = nVec, .combine = 'cbind') %:%
  foreach(icount(Q), .combine = 'c') %dopar% {
    X <- runif(n, min = 0, max = 1)
    L <- list()
    LMean <- list()
    for(i in 1:n){
      bk1Vec <- -rbinom(d[1], 1, theta[1])*rbeta(d[1], shape1 = sin(pi*X[i]), shape2 = 1-sin(pi*X[i]))
      bk2Vec <- -rbinom(d[2], 1, theta[2])*rbeta(d[2], shape1 = sin(pi*X[i]), shape2 = 1-sin(pi*X[i]))
      bk3Vec <- -rbinom(d[3], 1, theta[3])*rbeta(d[3], shape1 = sin(pi*X[i]), shape2 = 1-sin(pi*X[i]))
      temp1 <- matrix(0, nrow = m[1], ncol = m[1])
      temp1[lower.tri(temp1)] <- bk1Vec
      temp1 <- temp1 + t(temp1)
      temp2 <- matrix(bk2Vec, nrow = m[1])
      temp3 <- matrix(0, nrow = m[2], ncol = m[2])
      temp3[lower.tri(temp3)] <- bk3Vec
      temp3 <- temp3 + t(temp3)
      temp <- rbind(cbind(temp1, temp2), cbind(t(temp2), temp3))
      diag(temp) <- -colSums(temp)
      L[[i]] <- temp

      idx <- which.min(abs(X[i]-xSeq))
      LMean[[i]] <- LMeanRefsbm4[[idx]]
    }
    LRoot <- list()
    for(i in 1:n){
      LRoot[[i]] <- shapes::rootmat(L[[i]])
    }
    glRootVec <- matrix(unlist(LRoot), ncol = sum(m)^2, byrow = TRUE) # n by m^2
    fit <- list()
    for(i in 1:n){
      wi <- sapply(X[i] - X, kh, bw4sd[which(nVec == n)])
      temp <- matrix(apply(glRootVec, 2, weighted.mean, wi), ncol = sum(m))
      fiti <- G_2(temp)
      fit[[i]] <- as.matrix(proj_sparse(fiti %*% t(fiti)))
    }
    errsd <- mean(sapply(1:n, function(i) sum((LMean[[i]]-fit[[i]])^2)))

    sim <- lnr(L, X, optns = list(metric = "power", alpha = alpha, bw = bw4[which(nVec == n)]))
    err <- mean(sapply(1:n, function(i) sum((LMean[[i]]-sim$fit[[i]])^2)))
    c(errsd, err)
  }
isesbm4sd <- isesbm4[seq(1, 2*Q, by = 2), ]
isesbm4 <- isesbm4[seq(2, 2*Q, by = 2), ]

misesbm4 <- colMeans(isesbm4)
fitsbm4 <- lm(log(misesbm4)~log(nVec))
summary(fitsbm4)
plot(log(misesbm4), log(nVec))

# misesbm4sd <- colMeans(isesbm4sd)
# misesbm4sd - misesbm4
# 50 100 200 500 1000
# 1.7505000 1.0556659 0.6183106 0.3146141 0.1883002 Dryden
# 1.4004693 0.8399016 0.4866702 0.2466061 0.1486243

###########################################################################
# Part III: Comparison between different types of regression and metrics. #
###########################################################################
# Comparison under simulation scenario I
set.seed(1)
bw3 <- c(0.013, 0.013, 0.010, 0.008, 0.007)
bw4 <- c(0.17, 0.15, 0.13, 0.11, 0.10)
ise <- foreach(n = nVec) %:%
  foreach(icount(Q), .combine = 'rbind') %dopar% {
    X <- runif(n, min = 0, max = 1)
    L <- list()
    LMean <- list()
    for(i in 1:n){
      LVec <- -rbeta(sum(d), shape1 = X[i], shape2 = 1-X[i])
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      L[[i]] <- temp
      
      LMeanVec <- -rep(X[i], sum(d))
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LMeanVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      LMean[[i]] <- temp
    }
    sim1 <- gnr(L, X)
    sim2 <- gnr(L, X, optns = list(metric = "power", alpha = alpha))
    sim3 <- lnr(L, X, optns = list(bw = bw3[which(nVec == n)]))
    sim4 <- lnr(L, X, optns = list(metric = 'power', alpha = alpha, bw = bw4[which(nVec == n)]))
    rowMeans(sapply(1:n, function(i) c(sum((LMean[[i]]-sim1$fit[[i]])^2), 
                                       sum((LMean[[i]]-sim2$fit[[i]])^2), 
                                       sum((LMean[[i]]-sim3$fit[[i]])^2), 
                                       sum((LMean[[i]]-sim4$fit[[i]])^2))))
  }
ise1 <- sapply(1:length(nVec), function(n) ise[[n]][, 1])
ise2 <- sapply(1:length(nVec), function(n) ise[[n]][, 2])
ise3 <- sapply(1:length(nVec), function(n) ise[[n]][, 3])
ise4 <- sapply(1:length(nVec), function(n) ise[[n]][, 4])

# Comparison under simulation scenario III
set.seed(1)
bw3 <- c(0.088, 0.075, 0.065, 0.054, 0.047)
bw4 <- c(0.070, 0.060, 0.052, 0.043, 0.038)
ise <- foreach(n = nVec) %:%
  foreach(icount(Q), .combine = 'rbind') %dopar% {
    X <- runif(n, min = 0, max = 1)
    L <- list()
    LMean <- list()
    for(i in 1:n){
      LVec <- -rbeta(sum(d), shape1 = sin(pi * X[i]), shape2 = 1 - sin(pi * X[i]))
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      L[[i]] <- temp
      
      LMeanVec <- -rep(sin(pi * X[i]), sum(d))
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LMeanVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      LMean[[i]] <- temp
    }
    sim1 <- gnr(L, X)
    sim2 <- gnr(L, X, optns = list(metric = "power", alpha = alpha))
    sim3 <- lnr(L, X, optns = list(bw = bw3[which(nVec == n)]))
    sim4 <- lnr(L, X, optns = list(metric = 'power', alpha = alpha, bw = bw4[which(nVec == n)]))
    rowMeans(sapply(1:n, function(i) c(sum((LMean[[i]]-sim1$fit[[i]])^2), 
                                       sum((LMean[[i]]-sim2$fit[[i]])^2), 
                                       sum((LMean[[i]]-sim3$fit[[i]])^2), 
                                       sum((LMean[[i]]-sim4$fit[[i]])^2))))
  }
ise1 <- sapply(1:length(nVec), function(n) ise[[n]][, 1])
ise2 <- sapply(1:length(nVec), function(n) ise[[n]][, 2])
ise3 <- sapply(1:length(nVec), function(n) ise[[n]][, 3])
ise4 <- sapply(1:length(nVec), function(n) ise[[n]][, 4])

####################################################################
# Part IV: Networks generated from Erdos-Renyi random graph model. #
####################################################################
# Simulation scenario I
set.seed(1)
iseer1 <- foreach(n = nVec, .combine = 'cbind') %:%
  foreach(icount(Q), .combine = 'c') %dopar% {
    X <- runif(n, min = 0, max = 1)
    L <- list()
    LMean <- list()
    for(i in 1:n){
      LVec <- -rbinom(sum(d), 1, N / sum(d))*rbeta(sum(d), shape1 = X[i], shape2 = 1 - X[i])
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      L[[i]] <- temp
      
      LMeanVec <- -rep(X[i] * N / sum(d), sum(d))
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LMeanVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      LMean[[i]] <- temp
    }
    sim <- gnr(L, X)
    mean(sapply(1:n, function(i) sum((LMean[[i]]-sim$fit[[i]])^2)))
  }

# miseer1 <- colMeans(iseer1)
# fiter1 <- lm(log(miseer1)~log(nVec))
# summary(fiter1)
# plot(log(miseer1), log(nVec))#################################

# Simulation scenario II
set.seed(1)
iseer2 <- foreach(n = nVec, .combine = 'cbind') %:%
  foreach(icount(Q), .combine = 'c') %dopar% {
    X <- runif(n, min = 0, max = 1)
    L <- list()
    LMean <- list()
    for(i in 1:n){
      LVec <- -rbinom(sum(d), 1, N / sum(d))*rbeta(sum(d), shape1 = X[i], shape2 = 1 - X[i])
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      eigenDecom <- eigen(temp)
      Lambda <- eigenDecom$values
      U <- eigenDecom$vectors
      L[[i]] <- U%*%diag(Lambda^(1/alpha))%*%t(U)

      LMeanVec <- -rep(X[i] * N / sum(d), sum(d))
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LMeanVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      eigenDecom <- eigen(temp)
      Lambda <- pmax(eigenDecom$values, 0)
      U <- eigenDecom$vectors
      temp <- U%*%diag(Lambda^(1/alpha))%*%t(U)
      # model$Update(q = -temp)
      # temp <- matrix(model$Solve()$x, ncol = sum(m))
      # temp <- (temp + t(temp)) / 2 # symmetrize
      # temp[temp > 0] <- 0 # off diagonal should be negative
      # diag(temp) <- 0
      # diag(temp) <- -colSums(temp)
      LMean[[i]] <- temp
    }
    sim <- gnr(L, X, optns = list(metric = "power", alpha = alpha))
    mean(sapply(1:n, function(i) sum((LMean[[i]]-sim$fit[[i]])^2)))
  }

# miseer2 <- colMeans(iseer2)
# fiter2 <- lm(log(miseer2)~log(nVec))
# summary(fiter2)
# plot(log(miseer2), log(nVec))

# Simulation scenario III
set.seed(1)
N <- 9
bw1 <- c(0.180, 0.157, 0.137, 0.114, 0.1)
iseer3 <- foreach(n = nVec, .combine = 'cbind') %:%
  foreach(icount(Q), .combine = 'c') %dopar% {
    X <- runif(n, min = 0, max = 1)
    L <- list()
    LMean <- list()
    for(i in 1:n){
      LVec <- -rbinom(sum(d), 1, N / sum(d))*rbeta(sum(d), shape1 = sin(pi * X[i]), shape2 = 1 - sin(pi * X[i]))
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      L[[i]] <- temp
      
      LMeanVec <- -rep(sin(pi * X[i]) * N / sum(d), sum(d))
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LMeanVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      LMean[[i]] <- temp
    }
    sim <- lnr(L, X, optns = list(bw = bw1[which(nVec == n)]))
    mean(sapply(1:n, function(i) sum((LMean[[i]]-sim$fit[[i]])^2)))
  }

# miseer3 <- colMeans(iseer3)
# fiter3 <- lm(log(miseer3)~log(nVec))
# summary(fiter3)
# plot(log(miseer3), log(nVec))

# Simulation scenario IV
set.seed(1)
xSeq <- seq(0, 1, length.out = 1001)
LMeanRefer4 <- list()
for(i in 1:1001) {
  x <- xSeq[i]
  y <- matrix(0, nrow = sum(m), ncol = sum(m))
  for(j in 1:10000){
    LMeanVec <- -rbinom(sum(d), 1, N / sum(d)) * rbeta(sum(d), shape1 = sin(pi * x), shape2 = 1 - sin(pi * x))
    temp <- matrix(0, nrow = sum(m), ncol = sum(m))
    temp[lower.tri(temp)] <- LMeanVec
    temp <- temp + t(temp)
    diag(temp) <- -colSums(temp)
    eigenDecom <- eigen(temp)
    Lambda <- pmax(eigenDecom$values, 0)
    U <- eigenDecom$vectors
    temp <- U%*%diag(Lambda^alpha)%*%t(U)
    y <- y + temp/10000
  }
  eigenDecom <- eigen(y)
  Lambda <- pmax(eigenDecom$values, 0)
  U <- eigenDecom$vectors
  y <- U%*%diag(Lambda^(1/alpha))%*%t(U)

  model$Update(q = -y)
  y <- matrix(model$Solve()$x, ncol = sum(m))
  y <- (y + t(y)) / 2 # symmetrize
  y[y > 0] <- 0 # off diagonal should be negative
  diag(y) <- 0
  diag(y) <- -colSums(y)
  LMeanRefer4[[i]] <- y
}
# Export LMeanRefer4 and xSeq to parallel workers
clusterExport(cl, c("LMeanRefer4", "xSeq"))

set.seed(1)
bw2 <- c(0.095, 0.082, 0.072, 0.060, 0.052)
iseer4 <- foreach(n = nVec, .combine = 'cbind') %:%
  foreach(icount(Q), .combine = 'c') %dopar% {
    X <- runif(n, min = 0, max = 1)
    L <- list()
    LMean <- list()
    for(i in 1:n){
      LVec <- -rbinom(sum(d), 1, N / sum(d)) * rbeta(sum(d), shape1 = sin(pi * X[i]), shape2 = 1 - sin(pi * X[i]))
      temp <- matrix(0, nrow = sum(m), ncol = sum(m))
      temp[lower.tri(temp)] <- LVec
      temp <- temp + t(temp)
      diag(temp) <- -colSums(temp)
      L[[i]] <- temp
      
      idx <- which.min(abs(X[i] - xSeq))
      LMean[[i]] <- LMeanRefer4[[idx]]
    }
    sim <- lnr(L, X, optns = list(metric = 'power', alpha = alpha, bw = bw2[which(nVec == n)]))
    mean(sapply(1:n, function(i) sum((LMean[[i]]-sim$fit[[i]])^2)))
  }

# miseer4 <- colMeans(iseer4)
# fiter4 <- lm(log(miseer4)~log(nVec))
# summary(fiter4)
# plot(log(miseer4), log(nVec))