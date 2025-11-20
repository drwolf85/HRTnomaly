pif <- function(dta, nt = 100L, nss = NULL,
                max_depth = 12L, threshold = 0.95,
                proximity_type = c("single", "paired", "pivotal"),
                dist_fun = NULL) {
  if (is.null(dist_fun)) {
    dist_fun <- function(x, y) sum((x - y)^2)
  } else {
    if (!is.function(dist_fun))
      stop("The argument `dist_fun` must be `function(x, y)`.")
  }
  prx <- switch(proximity_type, paired = 2L,
                pivotal = 3L, 1L) # "single" (by default value) = 1L
  max_depth <- as.integer(max_depth)
  if (max_depth < 1)
    stop("The argument `max_depth` must be a positive integer number.")
  if (!is.list(dta))
    stop("The argument `dta` must be a `list` object.")
  nt <- as.integer(nt)
  if (nt < 1) stop("Provide a positive number of proximity isolation trees.")
  n <- length(dta)
  if (is.null(nss)) nss <- n * 0.25
  nss <- as.integer(nss)
  if (nss < (prx + 3L) || length(dta) < (prx + 3L))
    stop("Provide more data points or increase the size of the subsamples.")
  if (nss > length(dta))
    nss <- length(dta)
  rnv <- new.env()
  rnv <- parent.env(rnv)
  s <- .Call("pif", dta, prx, nt, nss, max_depth,
             quote(dist_fun(dta[[i]], dta[[j]])),
             rnv, pakcage = "HRTnomaly")
  attr(dta, "flag") <- s > quantile(s, prob = threshold)
  attr(dta, "scores") <- s
  return(dta)
}
