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
Q <- 12
c1 <- rnorm(Q)
c2 <- rnorm(Q)
# c1 <- rep(1, Q)
# c2 <- rep(1, Q)
s <- exp(-0.01 * 0:N) * cos(2 * pi * 0.2 * 0:N) %o% c1 + exp(-0.02 * 0:N) * cos(2 * pi * 0.22 * 0:N) %o% c2
r <- 4
r3 <- 2
# r <- 2
SNR <- 30
sigma <- mean(sapply(1:Q, function(i) snr.to.sd(s[, i], SNR)))
# sigma <- 0.02

mult.signest.sc <- function(L, Mat = TRUE, Tens = TRUE) {
  # s.n <- s + rcnorm(N + 1, var = sigma^2)
  s.n <- s + rnorm(N + 1, sd = sigma)

  estimates <- array(dim = c(nrow(s), Q, 2))

  if (Mat) {
    s.mssa <- ssa(s.n, L = L, kind = "mssa")
    estimates[, , 1] <- reconstruct(s.mssa, groups = list(1:r))[[1]]
  }

  if (Tens) {
    hm <- reduce(apply(s.n, 2, hankel, L = L, simplify = "list"), cbind)
    ht <- fold(hm, 1, 2:3, modes = c(L, N + 1 - L + 1, Q))
    # capture.output({ ht.hosvd <- rTensor::hosvd(ht, ranks = rep(r, 3)) })
    capture.output({ ht.hooi <- rTensor::tucker(ht, ranks = c(r, r, r3)) })
    # ht.hat <- ttl(ht.hosvd$Z, ht.hosvd$U, 1:3)
    ht.hat <- ttl(ht.hooi$Z, ht.hooi$U, 1:3)
    estimates[, , 2] <- apply(ht.hat@data, 3, hankel)
  }
  estimates
}

R <- 1000
signal.tens <- s %o% rep(1, R)
L.mat <- 22; L.tens <- 20
se.mult.comp.res <- list()
system.time({
  set.seed(1)
  se.mult.comp.res$mat <- replicate(R, mult.signest.sc(L = L.mat, Tens = FALSE), simplify = "array")[,,1,]
  set.seed(1)
  se.mult.comp.res$tens <- replicate(R, mult.signest.sc(L = L.tens, Mat = FALSE), simplify = "array")[,,2,]
})
mat.res <- se.mult.comp.res$mat
tens.res <- se.mult.comp.res$tens
mean(rowMeans(sqrt(apply(abs(mat.res - signal.tens)^2, 1:2, mean))))
mean(rowMeans(sqrt(apply(abs(tens.res - signal.tens)^2, 1:2, mean))))

# R <- 1000
# res2 <- list(mat = numeric(length(4:(N - r + 1))), tens = numeric(length(4:(N - r + 1))))
# for (L in 4:(N - r + 2)) {
#   set.seed(1)
#   system.time({ se.mult.comp.res <- replicate(R, mult.signest.sc(L = L), simplify = "array") })
#   mat.res <- se.mult.comp.res[, , 1,]
#   tens.res <- se.mult.comp.res[, , 2,]
#   res2$mat[L - 3] <- mean(rowMeans(sqrt(apply(abs(mat.res - signal.tens)^2, 1:2, mean))))
#   res2$tens[L - 3] <- mean(rowMeans(sqrt(apply(abs(tens.res - signal.tens)^2, 1:2, mean))))
# }
# which.min(res2$mat)
# which.min(res2$tens)