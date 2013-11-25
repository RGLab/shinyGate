get_tmix_filter <- function (fr, xChannel, yChannel, filterId = "", K = 2, usePrior = "no", 
  prior = list(NA), trans = 0, plot = FALSE, target = NULL, 
  transitional = FALSE, quantile = 0.9, translation = 0.25, 
  transitional_angle = NULL, min = NULL, max = NULL, ...) 
{
  if (!is.null(target)) {
    target <- as.numeric(target)
    if (length(target) != 2) {
      warning("The 'target' location must be a numeric vector of length 2.\n               Using largest cluster instead...")
      target <- NULL
    }
  }
  if (!(is.null(min) && is.null(max))) {
    fr <- .truncate_flowframe(fr, channels = c(xChannel, 
      yChannel), min = min, max = max)
  }
  if (usePrior == "yes" && identical(prior, list(NA))) {
    prior <- prior_flowClust(fr = fr, channels = c(xChannel, 
      yChannel), K = K)
  }
  tmix_filter <- tmixFilter(filterId, c(xChannel, yChannel), 
    K = K, trans = trans, usePrior = usePrior, prior = prior, 
    ...)
  return(tmix_filter)
}


get_tmix_result <- function (fr, xChannel, yChannel, filterId = "", K = 2, usePrior = "no", 
  prior = list(NA), trans = 0, plot = FALSE, target = NULL, 
  transitional = FALSE, quantile = 0.9, translation = 0.25, 
  transitional_angle = NULL, min = NULL, max = NULL, ...) 
{
  if (!is.null(target)) {
    target <- as.numeric(target)
    if (length(target) != 2) {
      warning("The 'target' location must be a numeric vector of length 2.\n               Using largest cluster instead...")
      target <- NULL
    }
  }
  if (!(is.null(min) && is.null(max))) {
    fr <- .truncate_flowframe(fr, channels = c(xChannel, 
      yChannel), min = min, max = max)
  }
  if (usePrior == "yes" && identical(prior, list(NA))) {
    prior <- prior_flowClust(fr = fr, channels = c(xChannel, 
      yChannel), K = K)
  }
  tmix_filter <- tmixFilter(filterId, c(xChannel, yChannel), 
    K = K, trans = trans, usePrior = usePrior, prior = prior, 
    ...)
  tmix_results <- try(filter(fr, tmix_filter), silent = TRUE)
  return(tmix_results)
}
