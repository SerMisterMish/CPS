library(rTensor)
library(tensr)
library(purrr)
library(Rssa)
library(Metrics)

# rcnorm <- function(n, mean = 0, var = 1) {
#   sqrt(var / 2) * rnorm(n, mean = Re(mean)) + 1i * sqrt(var / 2) * rnorm(n, mean = Im(mean))
# }

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
# w1 <- c(0.2, 0.18)
# w2 <- c(0.22, 0.24)
# l1 <- sample(w1, Q, replace = TRUE)
# l2 <- sample(w2, Q, replace = TRUE)
# s <- exp(-0.01 * 0:N) * cos(2 * pi * 0.2 * 0:N) %o% c1 + exp(-0.02 * 0:N) * cos(2 * pi * 0.22 * 0:N) %o% c2
s <- sapply(1:Q, function(i) exp(-0.01 * 0:N) * cos(2 * pi * 0.2 * 0:N + pi * i/ 6) * c1[i] +
                            exp(-0.02 * 0:N) * cos(2 * pi * 0.22 * 0:N + pi * i / 9) * c2[i])
# s <- sapply(1:Q, function(i) exp(-0.01 * 0:N) * cos(2 * pi * l1[i] * 0:N) * c1[i] +
#                             exp(-0.02 * 0:N) * cos(2 * pi * l2[i] * 0:N) * c2[i])
# r <- 8
r <- 4
# r3 <- 8
r3 <- 4
# r3 <- 2
SNR <- 30
# sigma <- mean(sapply(1:Q, function(i) snr.to.sd(s[, i], SNR)))
sigma <- 0.02

mult.signest.sc <- function(L, Mat = TRUE, Tens = TRUE) {
  # s.n <- s + rcnorm(N + 1, var = sigma^2)
  s.n <- s + rnorm((N + 1) * Q, sd = sigma)

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
# L.mat <- 22; L.tens <- 20 # - equal phases
L.mat <- 21; L.tens <- 21 # - linear phases
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
# signal.tens <- s %o% rep(1, R)
# res2 <- list(mat = numeric(length(4:(N - r + 2))), tens = numeric(length(4:(N - r + 2))))
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