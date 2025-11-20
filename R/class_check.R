class_check <- function(pred, truth) {
  tb <- table(pred, truth)
  attr(tb, "overall") <- sum(diag(tb)) / sum(tb)
  attr(tb, "recall") <- diag(apply(tb, 1, function(xx) xx / sum(xx)))
  attr(tb, "precision") <- diag(apply(tb, 2, function(xx) xx / sum(xx)))
  attr(tb, "f1-score") <- 1 /
    (0.5 / attr(tb, "recall") + 0.5 / attr(tb, "precision"))
  tb <- unclass(tb)
  class(tb) <- "checkwise"
  return(tb)
}
