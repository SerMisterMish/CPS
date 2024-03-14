library(rTensor)
library(tensr)
library(purrr)
library(Rssa)
library(Metrics)

tssa3 <- function(s, I = (length(s) + 2) %/% 3, L = I) {
  X <- tens3(s, I, L)
  result <- list()
  result$X <- X
  result$modes <- list(I = I, L = L, J = length(s) - I - L + 2)
  result$hosvd <- rTensor::hosvd(X)
  return(result)
}

# tens3.old <- function(s, I, L) {
#   require("rTensor")
#   v <- as.vector(s)
#   N <- length(v)
#   J <- N - I - L + 2
#   X <- array(NA, c(I, L, J))
#   for (i in 1:J) {
#     X[, , i] <- outer(1:I, 1:L, function(x, y) s[i + x + y - 2])
#   }
#   as.tensor(X)
# }

tens3 <- function(s, I, L) {
  require("rTensor")
  v <- as.vector(s)
  N <- length(v)
  J <- N - I - L + 2
  X <- sapply(1:J, function(j) Rssa::hankel(s[j:(j + I + L - 2)], I), simplify = "array")
  as.tensor(X)
}

# diagonal averaging i + j + k = const, const = 3:(I+L+J)
reconstruct.group3 <- function(X.tens) {
  X <- X.tens@data
  I <- length(X[, 1, 1])
  L <- length(X[1, , 1])
  J <- length(X[1, 1,])
  s <- vector(mode = "numeric", length = I + L + J - 2)
  for (C in 3:(I + L + J)) {
    sum <- 0
    count <- 0
    for (i in 1:(C - 2)) {
      for (l in 1:(C - 1 - i)) {
        if (i <= I && l <= L && C - i - l <= J) {
          sum <- sum + X[i, l, C - i - l]
          count <- count + 1
        }
      }
    }
    s[C - 2] <- sum / count
  }
  return(s)
}

make.group <- function(hosvd, group = 1:min(hosvd$Z@modes)) {
  ttl(hosvd$Z[group, group, group, drop = FALSE], list(
    as.matrix(hosvd$U[[1]][, group]),
    as.matrix(hosvd$U[[2]][, group]),
    as.matrix(hosvd$U[[3]][, group])),
      1:3)
}

make.groupHOOI <- function(X, group) {
  r <- length(group)
  hooi_x <- hooi(X@data, r = c(r, r, r), itermax = 1500, tol = 1.e-7)
  G <- hooi_x$G
  U <- hooi_x$U
  ## Reconstruct the hooi approximation.
  X_approx <- atrans(G, U)
  as.tensor(X_approx)
}

t3.reconstruct <- function(p, groups) {
  stopifnot(is.list(groups))
  lapply(lapply(groups, make.group, hosvd = p$hosvd), reconstruct.group3)
}

t3.reconstructHOOI <- function(p, groups) {
  require("tensr")
  stopifnot(is.list(groups))
  reconstruct.group3(make.groupHOOI(p$X, group = groups[[1]]))
}

rcnorm <- function(n, mean = 0, var = 1) {
  sqrt(var / 2) * rnorm(n, mean = Re(mean)) + 1i * sqrt(var / 2) * rnorm(n, mean = Im(mean))
}

snr.to.sd <- function(signal, SNR) {
  sqrt(mean(abs(signal)^2) / 10^(SNR / 10))
}

sd.to.snr <- function(signal, sd) {
  10 * log10(mean(abs(signal)^2) / sd^2)
}

fnorm <- function(tensr) { sqrt(sum(abs(tensr@data)^2)) }

set.seed(1)

N <- 24
# s <- exp((-0.01 + 2i * pi * 0.2) * 0:N) + exp((-0.02 + 2i * pi * 0.22) * 0:N)
s <- exp(-0.01 * 0:N) * cos(2 * pi * 0.2 * 0:N) + exp(-0.02 * 0:N) * cos(2 * pi * 0.22 * 0:N)
# s <- exp((2i * pi * 0.2) * 0:N) + exp((2i * pi * 0.22) * 0:N)
r <- 4
# r <- 2
SNR <- 30
sigma <- snr.to.sd(s, SNR)
# sigma <- 0.03

signest.sc <- function(I, L.t, L.m) {
  # s.n <- s + rcnorm(N + 1, var = sigma^2)
  s.n <- s + rnorm(N + 1, sd = sigma)

  estimates <- matrix(nrow = 2, ncol = length(s))

  if (!missing(L.m)) {
    hm <- hankel(s.n, L.m)
    hm.svd <- svd(hm, r, r)
    hat.hm <- hm.svd$u %*% diag(hm.svd$d[1:r]) %*% Conj(t(hm.svd$v))
    estimates[1,] <- hankel(hat.hm)
  }

  if (!missing(I) && !missing(L.t)) {
    ht <- tens3(s.n, I, L.t)
    # capture.output({ ht.hosvd <- rTensor::hosvd(ht, ranks = rep(r, 3)) })
    capture.output({ ht.hooi <- rTensor::tucker(ht, ranks = rep(r, 3)) })
    # ht.hat <- ttl(ht.hosvd$Z, ht.hosvd$U, 1:3)
    ht.hat <- ttl(ht.hooi$Z, ht.hooi$U, 1:3)
    estimates[2,] <- reconstruct.group3(ht.hat)
  }

  estimates
}

R <- 1000
signal.tens <- t(s %o% rep(1, 2)) %o% rep(1, R)
system.time({ se.comp.res <- replicate(R, signest.sc(I = 14, L.t = 8, L.m = 16)) })
rowMeans(sqrt(apply(abs(se.comp.res - signal.tens)^2, 1:2, mean)))