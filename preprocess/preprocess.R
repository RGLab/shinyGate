library(flowWorkspace)
library(data.table.extras)
library(RColorBrewer)
library(dichromat)
library(hexbin)
library(flowWorkspace)
library(openCyto)
library(graph)
library(Kmisc)

file <- "data/LoveLab-FlowFrames-Logicle.rds"
dat <- readRDS(file)

## go back to matrices for now
dat <- lapply(dat, function(x) {
  as.data.frame( exprs(x) )
})
  
## pseudo-gates: SSC-A, FSC-A
dat <- lapply(dat, function(x) {
  x$`SSC-A` <- logicleTransform()(rpois(nrow(x), lambda=1E3))
  x$`FSC-A` <- logicleTransform()(rpois(nrow(x), lambda=1E3))
  return(x)
})

m <- lapply(dat, colnames)
common <- Reduce(intersect, m)

## turn into a flowSet
fs <- flowSet(lapply(dat, function(x) {
  flowFrame( as.matrix(x[ names(x) %in% common ]) )
}))

## convert factors to characters in the pData of each ff
fs <- fsApply(fs, function(ff) {
  pData( parameters(ff) ) <- 
    factor_to_char( pData( parameters(ff) ) )
  return(ff)
})

## and make an empty GatingSet
gs <- GatingSet(fs)
system( paste("rm -rf", normalizePath("gs/LoveLab")) )
save_gs(gs, path="gs/LoveLab")
