library(openCyto)
library(flowClust)

set.seed(123)

source("test/get_tmix_filter.R")

gs <- load_gs("gs/LoveLab-Debug")
fs <- getData(gs, "CD3+")

## get the prior
prior <- prior_flowClust(fs, channels=c("CD4", "CD8"), K=2, nu0=2)

## no-prior gates
# gates_noprior <- flowClust.2d(fs, "CD4", "CD8", K=2, trans=1, plot=TRUE) ## doesn't work
gates_noprior_t0 <- fsApply(fs, function(ff) {
  flowClust.2d(ff, "CD4", "CD8", K=2, trans=0, target=c(3.5, 2), plot=TRUE, lambda=1-1E-5)
})

gates_noprior_t1 <- fsApply(fs, function(ff) {
  flowClust.2d(ff, "CD4", "CD8", K=2, trans=1, target=c(3.5, 2), plot=TRUE, lambda=1-1E-5)
})

## gates with prior
gates_prior_t0 <- fsApply(fs, function(ff) {
  flowClust.2d(ff, "CD4", "CD8", K=2, plot=TRUE, trans=0, target=c(3.5, 2), usePrior="yes", prior=prior, lambda=1-1E-5)
})

gates_prior_t1 <- fsApply(fs, function(ff) {
  flowClust.2d(ff, "CD4", "CD8", K=2, plot=TRUE, trans=1, target=c(3.5, 2), usePrior="yes", prior=prior, lambda=1 - 1E-5)
})

## noprior_t0
xyplot( CD8 ~ CD4, data=getData(gs, "CD3+"), filter=gates_noprior_t0, smooth=FALSE, xbin=32 )

## noprior_t1
xyplot( CD8 ~ CD4, data=getData(gs, "CD3+"), filter=gates_noprior_t1, smooth=FALSE, xbin=32 )

## prior_t0
xyplot( CD8 ~ CD4, data=getData(gs, "CD3+"), filter=gates_prior_t0, smooth=FALSE, xbin=32 )

## prior_t1
xyplot( CD8 ~ CD4, data=getData(gs, "CD3+"), filter=gates_prior_t1, smooth=FALSE, xbin=32 )

## get the tmix_filters, look for differences
## play with the prior first
fs <- getData(gs, "CD3+")
prior <- prior_flowClust(fs, channels=c("CD4", "CD8"), K=2)
fs <- fs[1]
tmix_prior_t1 <- fsApply(fs, function(ff) {
  list(
    get_tmix_result(ff, "CD4", "CD8", K=2, plot=TRUE, trans=1, target=c(3.5, 2), usePrior="yes", prior=prior, tol.init=1E-3),
    get_tmix_result(ff, "CD4", "CD8", K=2, plot=TRUE, trans=1, target=c(3.5, 2), usePrior="yes", prior=prior),
    get_tmix_result(ff, "CD4", "CD8", K=2, plot=TRUE, trans=1, target=c(3.5, 2), usePrior="yes", prior=prior, tol.init=1E-5),
    get_tmix_result(ff, "CD4", "CD8", K=2, plot=TRUE, trans=1, target=c(3.5, 2), usePrior="yes", prior=prior, tol.init=1E-5),
    get_tmix_result(ff, "CD4", "CD8", K=2, plot=TRUE, trans=1, target=c(3.5, 2), usePrior="yes", prior=prior, min.count=100, max.count=100),
    get_tmix_result(ff, "CD4", "CD8", K=2, plot=TRUE, trans=1, target=c(3.5, 2), usePrior="yes", prior=prior, level=0.8),
    get_tmix_result(ff, "CD4", "CD8", K=2, plot=TRUE, trans=1, target=c(3.5, 2), usePrior="yes", prior=prior, level=0.8)
  )
})

sapply(tmix_prior_t1, function(x) {
  lapply(x, function(x) x@lambda)
})

output <- sapply(tmix_prior_t1, function(x) {
  c(x@BIC, x@logLike)
})

## try to get a difference in the good, bad filters
debug(flowClust.2d)
gate <- flowClust.2d( getData(gs[[2]], "CD3+"), "CD4", "CD8", K=2, plot=TRUE, trans=1, target=c(3.5, 2), usePrior="yes", prior=prior)
xyplot( CD8 ~ CD4, getData(gs[[2]], "CD3+"), filter=gate )

## for the good version, we have lambda > 1

gate <- flowClust.2d( getData(gs[[4]], "CD3+"), "CD4", "CD8", K=2, plot=TRUE, trans=1, target=c(3.5, 2), usePrior="yes", prior=prior)
xyplot( CD8 ~ CD4, getData(gs[[4]], "CD3+"), filter=gate )

