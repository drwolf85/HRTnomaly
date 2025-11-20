\name{print.checkwise}
\Rdversion{1.1}
\alias{print.checkwise}
\encoding{UTF-8}

\title{A Method to Print the Accuracy of Outlier Classification Results}

\description{
The function print the confusion matrix and accuracy results previously computed with the function \code{class_check}.
}

\usage{
  \method{print}{checkwise}(x, confusion = FALSE, ...)
}

\arguments{
  \item{x}{An S3 object of the class \code{checkwise}, typically computed with the function \code{class_check}.}
  \item{confusion}{A logical value, which is \code{FALSE} by default. If \code{TRUE}, the confusion matrix is printed after showing all accuracy metrics.}
  \item{...}{Additional arguments to pass to the function \code{cat}.}
}

\details{
  The function computes the confusion matrix using the function \code{table}. True positive and false negative are successively evaluated to compute overall accuracy, recall, precision, and F1-scores.
}

\value{An S3 class named \code{checkwise} with the confusion matrix, and other accuracy metrics appended as attribues.
\describe{
  \item{\code{attr(, "overall")}}{A numeric value between zero and one with the overall accuracy.}
  \item{\code{attr(, "recall")}}{A numeric vector of values between zero and one with the recall index for regular and outlier cells.}
  \item{\code{attr(, "precision")}}{A numeric vector of values between zero and one with the precision index for regular and outlier cells.}
  \item{\code{attr(, "f1-score")}}{A numeric vector of values between zero and one with the F1-scores for regular and outlier cells.}
}
}

\author{Luca Sartore \email{drwolf85@gmail.com}}

\examples{
# Load the package
library(HRTnomaly)
set.seed(2025L)
# Load the 'toy' data
data(toy)
# Detect cellwise outliers using Bayesian Analysis
res <- cellwise(toy[sample.int(100), ], 0.5, 10L)
class_check(res$outlier, res$anomaly_flag != "")
}

\keyword{outliers}
\keyword{distribution}
\keyword{probability}
