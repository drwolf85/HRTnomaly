fuzzyHRT <- function(a, contamination = 0.08) {
  ## Historical and zero check
  hScore <- .C("history_check", double(nrow(a)), double(nrow(a)),
               as.double(a$current_value_num),
               as.double(a$pred_value), nrow(a),
               NAOK = TRUE, DUP = TRUE, PACKAGE = "HRTnomaly")[1L:2L]
  zScore <- hScore[[2L]]
  hScore <- hScore[[1L]]

  ## Tail-check
  dtac <- a[, c("strata", "unit_id", "master_varname", "current_value_num")] %>%
    pivot_wider(names_from = any_of("master_varname"),
                values_from = matches("current_value_num"))
  dtac[dtac <= 0] <- NA
  dtal <- as.matrix(log(dtac[, -1L:-2L]))

  gr <- factor(dtac$strata)
  # Smat <- double(prod(dim(dtal)))
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
  dtar <- dtac %>% pivot_longer(2 + seq_len(ncol(dtal)), values_drop_na = TRUE)

  dtar <- dtac %>% pivot_longer(cols = 3:dim(dtac)[2],
                                names_to = "master_varname",
                                values_to = "rScore")
  dtar <- left_join(a, dtar)

  a$score <- zScore * hScore * dtar$rScore
  th <- quantile(a$score[a$score != 0], contamination, na.rm = TRUE)
  a$outlier <- a$score < th & a$score != 0

  return(a)
}
