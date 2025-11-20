print.checkwise <- function(x, confusion = FALSE, ...) {
  cat("  Overall accuracy:\n", ...)
  cat(attr(x, "overall"), "\n", sep = "", ...)
  cat("  Recall:\n", ...)
  print(attr(x, "recall"))
  cat("  Precision:\n", ...)
  print(attr(x, "precision"))
  cat("  F1-score:\n", ...)
  print(attr(x, "f1-score"))
  if (confusion) {
    cat("  Confusion:\n", ...)
    print(as.table(x))
  }
}
