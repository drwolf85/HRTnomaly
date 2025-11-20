\name{class_check}
\Rdversion{1.1}
\alias{class_check}
\encoding{UTF-8}

\title{Evaluate the Accuracy of Outlier Classification Results}

\description{
The function computes the confusion matrix between the logical output of an outlier detection algorithm and a reference (ground-truth) logical vector. 
The function also calculates the overal accuracy of the results from the confusion matrix, including recall, precision, and F1-scores for the two classes (regular, versus outlier). 
}

\usage{
class_check(pred, truth)
}

\arguments{
  \item{pred}{A logical vector with the classification output from an anomaly detection algorithm.}
  \item{truth}{A logical vector with the observed classification as a reference (or ground truth).}
}

\details{
  The function computes the confusion matrix using the function \code{table}. True positive and false negative are successively evaluated to compute overall accuracy, recall, precision, and F1-scores.
}

\value{An S3 class named \code{checkwise} with the confusion matrix, and other accuracy metrices appended as attribues.
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
# Detect cellwise outliers using Cellwise Analysis
res <- cellwise(toy[sample.int(100), ], 0.05, 2L)
class_check(res$outlier, res$anomaly_flag != "")
}

\keyword{outliers}
\keyword{distribution}
\keyword{probability}
