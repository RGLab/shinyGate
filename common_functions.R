consistency <- function(gs) {
  allNodes <- lapply(gs, getNodes)
  intersection <- Reduce(intersect, allNodes)
  union <- Reduce(union, allNodes)
  node <- setdiff(union, intersection)
  if (length(node)) {
    bad_samples <- names(allNodes)[unlist(lapply(allNodes, function(x) {
      !(node %in% x)
    }))]
    good_samples <- setdiff( names(allNodes), bad_samples )
    parent <- getParent(gs[[good_samples[1]]], node, isPath=TRUE)
    
  } else {
    bad_samples <- parent <- character()
  }
  
  return (list(
    is_consistent=identical(intersection, union),
    node=node,
    parent=parent,
    samples=bad_samples
  ))
}
