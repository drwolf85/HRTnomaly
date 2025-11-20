bayeswise <- function(a, prior = NULL, epochs = 1000L) {
  if (!is.numeric(epochs) || !is.finite(epochs))
    stop("The argument `epochs` must be a finite integer number")
  # Check if the dataset `a` contains values for the prior of each cell
  if (is.null(prior)) {
    if ("prior" %in% colnames(a)) {
      prior <- 1 - a$prior # Cell-level prior probability for regular cases!!!
    } else {
      # Non informative prior if `a` has no prior probabilities for each cell
      prior <- 0.5 # This is cell-level prior probability for regular cases!!!
    }
  } else {
    prior <- 1 - prior # Cell-level prior probability for regular cases!!!
  }
  a$prior <- prior

  storage.mode(epochs) <- "integer"
  wh <- c("strata", "unit_id", "master_varname")
  dtac <- a[, c(wh, "current_value_num")] %>%
    pivot_wider(names_from = any_of("master_varname"),
                values_from = matches("current_value_num"))
  ord <- order(dtac$strata)
  dtac <- dtac[ord, ]

  dtap <- a[, c(wh, "pred_value")] %>%
    pivot_wider(names_from = any_of("master_varname"),
                values_from = matches("pred_value"))
  ord <- order(dtap$strata)
  dtap <- dtap[ord, ]

  dtas <- a[, c(wh, "prior")] %>%
    pivot_wider(names_from = any_of("master_varname"),
                values_from = matches("prior"))
  ord <- order(dtas$strata)
  dtas <- dtas[ord, ]

  xc <- as.matrix(dtac[, -1L:-2L]) # Current
  xp <- as.matrix(dtap[, -1L:-2L]) # Past
  s <- as.matrix(dtas[, -1L:-2L]) # Prior Prob.
  g <- matrix(0, nrow(xc), ncol(xc))
  z <- matrix(0, nrow(xc), ncol(xc))
  h <- matrix(0, nrow(xc), ncol(xc))
  r <- matrix(0, nrow(xc), ncol(xc))
  t <- matrix(0, nrow(xc), ncol(xc))
  storage.mode(g) <- "integer"
  storage.mode(s) <- "double"
  scores <- .C("bayeswise", s = s, G = g, z = z, h = h, r = r, t = t, xc, xp,
               dim(xc), epochs, NAOK = TRUE, PACKAGE = "HRTnomaly")
  scores$s <- as.data.frame(scores$s)
  scores$G <- as.data.frame(scores$G)
  scores$z <- as.data.frame(scores$z)
  scores$h <- as.data.frame(scores$h)
  scores$r <- as.data.frame(scores$r)
  scores$t <- as.data.frame(scores$t)
  scores$s$unit_id <- dtac$unit_id
  scores$G$unit_id <- dtac$unit_id
  scores$z$unit_id <- dtac$unit_id
  scores$h$unit_id <- dtac$unit_id
  scores$r$unit_id <- dtac$unit_id
  scores$t$unit_id <- dtac$unit_id
  names(scores$s) <- c(colnames(xc), "unit_id")
  names(scores$G) <- c(colnames(xc), "unit_id")
  names(scores$z) <- c(colnames(xc), "unit_id")
  names(scores$h) <- c(colnames(xc), "unit_id")
  names(scores$r) <- c(colnames(xc), "unit_id")
  names(scores$t) <- c(colnames(xc), "unit_id")
  scores$s <- scores$s %>%
    pivot_longer(seq_len(ncol(xc)), values_drop_na = TRUE)
  scores$G <- scores$G %>%
    pivot_longer(seq_len(ncol(xc)), values_drop_na = TRUE)
  scores$z <- scores$z %>%
    pivot_longer(seq_len(ncol(xc)), values_drop_na = TRUE)
  scores$h <- scores$h %>%
    pivot_longer(seq_len(ncol(xc)), values_drop_na = TRUE)
  scores$r <- scores$r %>%
    pivot_longer(seq_len(ncol(xc)), values_drop_na = TRUE)
  scores$t <- scores$t %>%
    pivot_longer(seq_len(ncol(xc)), values_drop_na = TRUE)
  names(scores$s)[-1L] <- c("master_varname", "score")
  names(scores$G)[-1L] <- c("master_varname", "outlier")
  names(scores$z)[-1L] <- c("master_varname", "zScore")
  names(scores$h)[-1L] <- c("master_varname", "hScore")
  names(scores$r)[-1L] <- c("master_varname", "rScore")
  names(scores$t)[-1L] <- c("master_varname", "tScore")
  a <- a %>% inner_join(scores$z, by = c("master_varname", "unit_id"))
  a <- a %>% inner_join(scores$h, by = c("master_varname", "unit_id"))
  a <- a %>% inner_join(scores$r, by = c("master_varname", "unit_id"))
  a <- a %>% inner_join(scores$t, by = c("master_varname", "unit_id"))
  a <- a %>% inner_join(scores$s, by = c("master_varname", "unit_id"))
  a <- a %>% inner_join(scores$G, by = c("master_varname", "unit_id"))
  a$outlier <- as.logical(a$outlier)
  return(a)
}
