consistency <- function(gs) {
  allNodes <- lapply(gs, getNodes)
  intersection <- Reduce(intersect, allNodes)
  union <- Reduce(union, allNodes)
  node <- setdiff(union, intersection)
  if (length(node)) {
    parent <- getParent(gs[[1]], node, isPath=TRUE)
    samples <- names(allNodes)[unlist(lapply(allNodes, function(x) {
      !(node %in% x)
    }))]
  } else {
    samples <- parent <- character()
  }
  
  return (list(
    is_consistent=identical(intersection, union),
    node=node,
    parent=parent,
    samples=samples
  ))
}
