#!/usr/bin/env Rscript
# Checks after model choice (manual 4.1, 4.3): prior-check PCA, model-choice LDA, OOB error vs trees.
# Prints fit numbers and writes 5 PNGs to results/figures/: one per axis pair (PC1-2, PC3-4, LD1-2, LD3-4), no titles,
# large labels, and the OOB curve. Needs only data.table (base graphics).
# Run from diyabc_rf/:  Rscript scripts/4_make_check_figs.R
library(data.table)
options(bitmapType="cairo")
out <- "results/figures"
cols <- c("#E41A1C","#377EB8","#4DAF4A","#000000","#FF7F00","#A65628","#F781BF","#999999","#984EA3","#66C2A5",
          "#FFD92F","#00CED1","#B3DE69")
# One scatter per axis pair: sims coloured by scenario, observed as a yellow diamond, legend outside on the right
scatter <- function(file, M, g, o, pr, xl, yl) {
  png(file, 1600, 1400, res=200); on.exit(dev.off())
  par(mar=c(5, 5.5, 1, 9), cex.lab=1.9, cex.axis=1.5, mgp=c(3.3, 1, 0))
  plot(M[,pr[1]], M[,pr[2]], pch=16, cex=.35, col=adjustcolor(cols[g],.35), xlab=xl, ylab=yl)
  points(o[pr[1]], o[pr[2]], pch=23, bg="yellow", cex=3, lwd=2.5)
  k <- length(labs)
  legend("topleft", inset=c(1.02, 0), xpd=TRUE, c(labs,"observed"), pch=c(rep(16,k),23), col=c(cols[1:k],"black"),
         pt.bg=c(rep(NA,k),"yellow"), pt.cex=c(rep(2,k),2.5), cex=1.5, bty="n")
}

## PCA (prior check)
pv <- scan("results/prior_check/pcaloc1_ACP.txt", nlines=1, quiet=TRUE)
pca <- fread("results/prior_check/pcaloc1_ACP.txt", skip=1, header=FALSE)
obs <- unlist(pca[V1==0, -1]); sim <- pca[V1!=0]
sc <- sim$V1; X <- as.matrix(sim[, -1])
labs <- paste0("S", seq_len(max(sc)))   # one label per scenario in header.txt (up to 13 colours)
cat(sprintf("PCA: variance PC1 %.3f PC2 %.3f\n", pv[3], pv[4]))
mahal_pct <- function(X, o, k) {         # percentile of obs Mahalanobis D2 among sims, first k axes
  X <- X[, 1:k, drop=FALSE]; mu <- colMeans(X); S <- cov(X)
  d <- mahalanobis(X, mu, S); do <- mahalanobis(o[1:k], mu, S); mean(d < do)
}
for (k in c(2, 5, 10, 20)) for (i in sort(unique(sc))) cat(sprintf("PCA %2d axes: obs Mahalanobis percentile within %s sims %.3f\n",
  k, labs[i], mahal_pct(X[sc==i,], obs, k)))
for (j in 1:2) cat(sprintf("PC%d: obs %.2f, sims range %.1f..%.1f, quantile in sims %.3f\n", j, obs[j], min(X[,j]), max(X[,j]), mean(X[,j] < obs[j])))
set.seed(1); idx <- sample(nrow(X), 20000)
for (pr in list(c(1,2), c(3,4)))
  scatter(file.path(out, sprintf("priorcheck_pca_PC%d_PC%d.png", pr[1], pr[2])), X[idx,], sc[idx], obs, pr,
          sprintf("PC%d (%.1f%%)",pr[1],100*pv[2+pr[1]]), sprintf("PC%d (%.1f%%)",pr[2],100*pv[2+pr[2]]))

## LDA (model choice)
lda <- fread("results/model_choice/modelchoice_all.lda", skip=1, header=FALSE)
n <- ncol(lda) - 1   # number of LDA axes: scenarios - 1
lo <- unlist(lda[1, 1:n]); L <- as.matrix(lda[-1, 1:n]); ls <- lda[-1][[n + 1]]
for (i in sort(unique(ls))) cat(sprintf("LDA %d axes: obs Mahalanobis percentile within %s sims %.3f\n",
  n, labs[i], mahal_pct(L[ls==i,], lo, n)))
set.seed(2); idx <- sample(nrow(L), 20000)
for (pr in if (n >= 4) list(c(1,2), c(3,4)) else list(c(1,2), c(2,3)))   # V5: 4 scenarios -> only 3 LDA axes
  scatter(file.path(out, sprintf("modelchoice_lda_LD%d_LD%d.png", pr[1], pr[2])), L[idx,], ls[idx], lo, pr,
          paste0("LD",pr[1]), paste0("LD",pr[2]))

## OOB error vs trees. Single-tree values are noisy: quote the error in the .confusion file.
e <- scan("results/model_choice/modelchoice_all.ooberror", quiet=TRUE)
png(file.path(out,"ooberror.png"), 1200, 800, res=150)
plot(e, type="l", ylim=c(0, max(e)), xlab="number of trees", ylab="OOB prior error rate", main="Error vs number of trees")
dev.off()
cat(sprintf("Mean error over the last 100 trees %.4f\n", mean(tail(e, 100))))
