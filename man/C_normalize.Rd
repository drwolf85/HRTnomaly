\name{C_normalize}
\alias{C_normalize}

\title{Native C Interface Wrapper for the Function \code{normalize}}

\description{
  These compiled symbols represent the direct entry points to the
  underlying C source code, masked behind the \code{C_} namespace fix.

  The function \code{normalize} performs the robust normalization by stratum of the data in each column.
}

%\usage{
%.C(C_normalize, dta, dim, gr, ng, res)
%}
%
%\arguments{
%  \item{dta}{An input dataset.}
%  \item{dim}{The dimensions of the input dataset (i.e., a vector with the number of rows and columns).}
%  \item{gr}{An integer vector with stratum IDs.}
%  \item{ng}{The total number of strata.}
%  \item{res}{A residual matrix assumed to be already initialized to zeros and to have the same dimensions of \code{dta}.}
%}
%
\keyword{internal}

