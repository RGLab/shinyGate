is_consistent <- function(gs) {
  allNodes <- lapply(gs, function(x) getNodes(x))
  if (!identical(
    Reduce(intersect, allNodes),
    Reduce(union, allNodes)
  )) {
    return (FALSE)
  } else {
    return (TRUE)
  }
}

## should only be one, if we are doing our job ...
inconsistent_node <- function(gs) {
  allNodes <- lapply(gs, function(x) getNodes(x))
  intersection <- Reduce(intersect, allNodes)
  union <- Reduce(union, allNodes)
  diff <- setdiff(union, intersection)
  return (diff)
}

is_missing_node <- function(gs, node) {
  unlist(lapply(gs, function(x) {
    !(node %in% getNodes(x))
  }))
}
