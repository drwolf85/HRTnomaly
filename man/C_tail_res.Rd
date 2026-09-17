\name{C_tail_res}
\alias{C_tail_res}

\title{Native C Interface Wrapper for the Function \code{tail_res}}

\description{
  These compiled symbols represent the direct entry points to the
  underlying C source code, masked behind the \code{C_} namespace fix.

  The function \code{tail_res} computes the residuals for checking distribution tails of each data column by stratum.
}

%\usage{
%.C(C_tail_res, dta, dim, gr, ng, tRes)
%}
%
%\arguments{
%  \item{dta}{An input dataset.}
%  \item{dim}{The dimensions of the input dataset (i.e., a vector with the number of rows and columns).}
%  \item{gr}{An integer vector with stratum IDs.}
%  \item{ng}{The total number of strata.}
%  \item{tRes}{A vector where to store the residuals for detecting tail outliers.}
%}
%
\keyword{internal}

