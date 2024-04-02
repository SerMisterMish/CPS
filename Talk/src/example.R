library(Rssa)
library(lattice)
N <- 47
s1 <- exp(-0.01 * 0:(N - 1)) * cos(2 * pi * 0:(N - 1) / 3)
s2 <- exp(-0.02 * 0:(N - 1)) * cos(2 * pi * 0:(N - 1) / 4)
s <- s1 + s2

set.seed(1)
noise <- rnorm(N, sd = 0.4)

s.ssa <- ssa(ts(s + noise))
s.rec.signal <- reconstruct(s.ssa, groups = list(`Reconstructed signal` = 1:4))
s.rec.signal$Signal <- s
pdf("./Talk/src/img/decomp.pdf", width = 8, height = 4.5)
plot(s.rec.signal, add.residuals = FALSE, add.original = FALSE, plot.method = "xyplot",
     superpose = TRUE, auto.key = list(columns = 2))
dev.off()
s.rec.cos1 <-  reconstruct(s.ssa, groups = list(`Reconstructed` = 1:2))
s.rec.cos1$`Original` <- s1
s.rec.cos2 <-  reconstruct(s.ssa, groups = list(`Reconstructed` = 3:4))
s.rec.cos2$`Original` <- s2
pdf("./Talk/src/img/rec.pdf", width = 10, height = 4.5)
pl1 <- plot(s.rec.cos1, add.residuals = FALSE, add.original = FALSE, plot.method = "xyplot",
     superpose = TRUE, auto.key = list(columns = 2), main = "3-period cos")
pl2 <- plot(s.rec.cos2, add.residuals = FALSE, add.original = FALSE, plot.method = "xyplot",
     superpose = TRUE, auto.key = list(columns = 2), main = "4-period cos")
print(pl1, split = c(1, 1, 2, 1), more = TRUE)
print(pl2, split = c(2, 1, 2, 1), more = FALSE)
dev.off()
s.parest <- parestimate(ssa(ts(s + noise)), groups = list(cos1 = 1:4), method = "esprit")