\name{C_tail_check}
\alias{C_tail_check}

\title{Native C Interface Wrapper for the Function \code{tail_check}}

\description{
  These compiled symbols represent the direct entry points to the
  underlying C source code, masked behind the \code{C_} namespace fix.

  The function \code{tail_check} checks distribution tails of each data column by stratum.
}

%\usage{
%.C(C_tail_check, dta, dim, gr, ng, tScore)
%}
%
%\arguments{
%  \item{dta}{An input dataset.}
%  \item{dim}{The dimensions of the input dataset (i.e., a vector with the number of rows and columns).}
%  \item{gr}{An integer vector with stratum IDs.}
%  \item{ng}{The total number of strata.}
%  \item{tScore}{A vector where to store the scores for tail outliers.}
%}
%
\keyword{internal}

