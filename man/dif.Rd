\name{dif}
\Rdversion{1.1}
\alias{dif}
\encoding{UTF-8}

\title{Deep Isolation Forest}

\description{
The function builds a deep isolation forest that uses fuzzy logic to determine if a record is anomalous or not.
The function takes a wide-format \code{data.frame} object as input and returns it with two appended vectors.
The first vector contains the anomaly scores as numbers between zero and one, and the second vector provides 
a set of logical values indicating whether the records are outliers (\code{TRUE}) or not (\code{FALSE}).
}

\usage{
dif(dta, nt = 100L, nss = NULL, threshold = 0.95)
}

\arguments{
  \item{dta}{A wide-format \code{data.frame} object with records (stored by row).}
  \item{nt}{Number of deep isolation trees to build to form the forest. By default, it is set to \code{100}.}
  \item{nss}{Number of subsamples used to build a single deep isolation tree. 
  If set (by default) to \code{NULL}, the program will randomly select 25\% of the records provided to the \code{dta} argument.}
  \item{threshold}{A number between zero and one used as a threshold when identifying outliers from the anomaly scores. 
  By default, this argument is set to \code{0.95}, so that 5\% of the records is going to be classified as anomalous.}
}

\details{
The argument \code{dta} is proivded as an object of class \code{data.frame}. 
This object is considered as a wide-format \code{data.frame}. 
The use of the R-packages \code{dplyr}, \code{purrr}, and \code{tidyr} is highly recommended to simplify the conversion of datasets between long and wide formats.
}

\value{The wide-format \code{data.frame} is provided as input data and contains extra columns, i.e., for both anomaly scores and the outlier flags.}

\author{Luca Sartore \email{drwolf85@gmail.com}}

\examples{
# Load the package
library(HRTnomaly)
set.seed(2025L)
# Detect outliers in the `iris` dataset
res <- dif(iris)
}

\keyword{outliers}
\keyword{distribution}
\keyword{probability}
