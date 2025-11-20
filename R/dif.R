dif <- function(dta, nt = 100L, nss = NULL, threshold = 0.95) {
	if (!is.data.frame(dta))
		stop("The argument `dta` must be a `data.frame` object.")
	nt <- as.integer(nt)
	if (nt < 1) stop("Provide a positive number of deep isolation trees")
	if (is.null(nss)) nss <- nrow(dta) * 0.25
	nss <- as.integer(nss)
	if (nss < 3 || nrow(dta) < 3) stop("Provide more data points or increase the size of the subsamples")
	if (nss > nrow(dta)) nss <- nrow(dta)
	whchr <- sapply(dta, is.character)
	dtchr <- as.data.frame(lapply(dta[, whchr], as.factor))
	dtchr <- if (prod(dim(dtchr)) > 0) model.matrix(~ . + 0, data = dtchr) else c()
	dtnum <- dta[, !whchr]
	dtnum <- if (prod(dim(dtnum)) > 0) model.matrix(~ . + 0, data = dtnum) else c()
	dtnum <- cbind(dtnum, dtchr)
	storage.mode(dtnum) <- "double"
	dimD <- dim(dtnum)
	s <- .C("dif", s = double(dimD[1]), dtnum, dimD, nt, nss, pakcage = "HRTnomaly")$s
	dta <- cbind.data.frame(dta, scores = s, flags = s > quantile(s, prob = threshold))
	return(dta)
}
