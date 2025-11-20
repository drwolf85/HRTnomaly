bayesHRT <- function(a, prior = NULL) {
  # Check if the dataset `a` contains values for the prior of each cell
  if (is.null(prior)) {
    if ("prior" %in% colnames(a)) {
        prior <- 1 - a$prior # This is cell-level prior probability for regular cases!!!
    } else {
      # Use non informative prior if the dataset `a` does not contain prior probabilities for each cell
      prior <- 0.5 # This is cell-level prior probability for regular cases!!!
    }
  } else {
    prior <- 1 - prior # This is cell-level prior probability for regular cases!!!
  }
  ## Historical residual and zero check
  hRes <- .C("history_res", double(nrow(a)), double(nrow(a)),
             as.double(a$current_value_num),
             as.double(a$pred_value), nrow(a),
             NAOK = TRUE, DUP = TRUE, PACKAGE = "HRTnomaly")[1L:2L]
  zScore <- hRes[[2L]] * prior # Used as priors for each cell
  hRes <- hRes[[1L]]
  dtah <- cbind.data.frame(a[, c("strata", "unit_id", "master_varname")], hRes, zScore)
  dtaz <- dtah[, -4L] %>% pivot_wider(names_from = any_of("master_varname"), values_from = matches("zScore"))
  dtah <- dtah[, -5L] %>% pivot_wider(names_from = any_of("master_varname"), values_from = matches("hRes"))
  dtah[, -1L:-2L][is.na(dtah[, -1L:-2L])] <- 0

  ## Tail-check
  dtac <- a[, c("strata", "unit_id", "master_varname", "current_value_num")] %>%
          pivot_wider(names_from = any_of("master_varname"), values_from = matches("current_value_num"))
  dtac[dtac <= 0] <- NA
  dtal <- as.matrix(log(dtac[, -1L:-2L]))

  gr <- factor(dtac$strata)

  tRes <- .C("tail_res", as.double(dtal), dim(dtal),
             gr, nlevels(gr), res = double(prod(dim(dtal))),
             NAOK = TRUE, PACKAGE = "HRTnomaly")$res
  tRes <- array(tRes, dim = dim(dtal))
  tRes[is.na(tRes)] <- 0

  ## Relational-check
  rRes <- 0
  dtae <- tRes
  rRes <- .C("relat_res", dtae = as.double(dtae),
             dim(dtal), PACKAGE = "HRTnomaly")$dtae
  rRes <- array(rRes, dim = dim(dtal))

  ## Putting things together using the highest posterior probability class
  pr_mat <- as.matrix(dtaz[, -1L:-2L])
  hRes <- as.matrix(dtah[, -1L:-2L])
  finals <- .C("post_results", as.double(pr_mat), integer(prod(dim(dtal))),
               dim(dtal), as.double(hRes), as.double(rRes), as.double(tRes),
               NAOK = TRUE, PACKAGE = "HRTnomaly")[1L:2L]
  pr_mat <- array(finals[[1L]], dim = dim(dtal)) # Outlier posterior probability!!!
  finals <- array(as.logical(finals[[2L]]), dim = dim(dtal)) # Outlier flag (if TRUE then outlier)

  ## Putting things together using a Fuzzy-Logic-Inspired procedure
  dtac[, -1L:-2L] <- array(pr_mat, dim = dim(dtal))
  # dtar <- dtac %>% pivot_longer(2 + seq_len(ncol(dtal)), values_drop_na = TRUE)
  dtar <- dtac %>% pivot_longer(cols = 3:dim(dtac)[2],
                                names_to = "master_varname",
                                values_to = "post_prob")
  a <- left_join(a, dtar)
  dtac[, -1L:-2L] <- finals
  dtar <- dtac %>% pivot_longer(cols = 3:dim(dtac)[2],
                                names_to = "master_varname",
                                values_to = "outlier")
  a <- left_join(a, dtar)
  return(a)
}
