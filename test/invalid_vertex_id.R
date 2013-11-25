library(openCyto)
library(flowWorkspace)

gs <- load_gs("gs/LoveLab-InvalidNode")
lapply(gs, function(x) getNodes(x))
fs <- getData(gs, "CD14-")
fs <- flowSet(lapply(gs, function(x) {
  getData(x, "CD14-")
}))
