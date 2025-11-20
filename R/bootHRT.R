bootHRT <- function(a, contamination = 0.08, boot_max_it = 1000L) {
  ## Historical and zero check
  hScore <- .C("history_check", double(nrow(a)), double(nrow(a)),
               as.double(a$current_value_num),
               as.double(a$pred_value), nrow(a),
               NAOK = TRUE, DUP = TRUE, PACKAGE = "HRTnomaly")[1L:2L]
  zScore <- hScore[[2L]]
  hScore <- hScore[[1L]]

  ## Tail-check
  dtac <- a[, c("strata", "unit_id", "master_varname", "current_value_num")] %>%
          pivot_wider(names_from = any_of("master_varname"), values_from = matches("current_value_num"))
  dtac[dtac <= 0] <- NA
  dtal <- as.matrix(log(dtac[, -1L:-2L]))
  
  gr <- factor(dtac$strata)
  
  tScore <- .C("tail_check", as.double(dtal), dim(dtal),
                gr, nlevels(gr), res = double(prod(dim(dtal))),
                NAOK = TRUE, PACKAGE = "HRTnomaly")$res
  
  ## Relational-check
  rScore <- 1
  dtae <- .C("normalize", as.double(dtal), dim(dtal),
              gr, nlevels(gr), res = double(prod(dim(dtal))),
              NAOK = TRUE, PACKAGE = "HRTnomaly")$res
  dtae[is.na(dtae)] <- 0
  rScore <- .C("relat_check", dtae = as.double(dtae),
                dim(dtal), PACKAGE = "HRTnomaly")$dtae
  rScore <- array(rScore, dim = dim(dtal))
  
  ## Putting things together using a Fuzzy-Logic-Inspired procedure
  dtac[, -1L:-2L] <- array(rScore * tScore, dim = dim(dtal))
  #dtar <- dtac %>% pivot_longer(2 + seq_len(ncol(dtal)), values_drop_na = TRUE)
  
  dtar <- dtac %>% pivot_longer(cols = 3:dim(dtac)[2], names_to = "master_varname", values_to = "rScore") #dtac %>% pivot_longer(5 + seq_len(ncol(dtal)), values_drop_na=TRUE)
  dtar <- left_join(a, dtar)

  ## Empirical Bayesian threshold distribution
  a$score <- zScore * hScore * dtar$rScore
  finScores <- na.omit(a$score[a$score!=0])
	storage.mode(contamination) <- "double"
	storage.mode(finScores) <- "double"
	storage.mode(boot_max_it) <- "integer"
	th_v <-.C("bayes_boot", th = double(boot_max_it), 
            boot_max_it, finScores, length(finScores), 
            contamination, PACKAGE = "HRTnomaly")$th
  mth <- mean(th_v)
  th_s <- c(quantile(th_v, c(0.25, 0.5, 0.75)), mth, mth + c(-1, 1) * sd(th_v))
  outly <- sapply(th_s, function(th) a$score < th & a$score != 0)
  colnames(outly) <- paste0("outlier_", c(1:3, 1:3), rep(c("qt", "mn"), each = 3))
  a <- cbind.data.frame(a, outly)
  attr(a, "thresholds") <- th_v

  return(a)
}
