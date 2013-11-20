library(ggplot2)
library(flowWorkspace)
library(data.table.extras)
library(RColorBrewer)
library(dichromat)
library(hexbin)
library(flowWorkspace)
library(openCyto)
library(graph)
library(Kmisc)

source("../../Reports/LoveLabSingleCell/common_functions.R")

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

## there are some gates that aren't completely in common --
## we generate 'pseudo-gates' for these guys
m <- lapply(dat, colnames)
inclusionPlot(m)

allGates <- unique(unlist(m))
dat <- lapply(dat, function(x) {
  for (gate in allGates) {
    if (gate %nin% names(x)) {
      x[gate] <- 0
    }
  }
  return (x[allGates])
})

## turn into a flowSet
fs <- flowSet(lapply(dat, function(x) {
  flowFrame( as.matrix(x) )
}))

## convert factors to characters in the pData of each ff
fs <- fsApply(fs, function(ff) {
  pData( parameters(ff) ) <- 
    factor_to_char( pData( parameters(ff) ) )
  return(ff)
})

## make nicer names for fs
nm <- sampleNames(fs)

## they used an 'O' (letter-o) instead of a '0' (number-zero) in the names
nm <- gsub("TGO4", "TG04", nm)
plate <- sapply( strsplit(nm, "_"), function(x) paste(x[1:3], collapse="_") )
id <- gsub(".*_", "", nm)

sampleNames(fs) <- paste(
  paste0("P", as.integer(as.factor(plate))),
  id,
  sep=" : "
)

## and make an empty GatingSet
gs <- GatingSet(fs)
system( paste("rm -rf", normalizePath("gs/LoveLab")) )
save_gs(gs, path="gs/LoveLab")
