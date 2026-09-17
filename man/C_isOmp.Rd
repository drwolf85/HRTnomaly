\name{C_isOmp}
\alias{C_isOmp}

\title{Native C Interface Wrapper for the Function \code{isOmp}}

\description{
  These compiled symbols represent the direct entry points to the
  underlying C source code, masked behind the \code{C_} namespace fix.

  The function \code{isOmp} checks whether the package is configured to use OpenMP for parallel computing. It returns a logical value (\code{TRUE} if OpenMP is enabled, \code{FALSE} otherwise).
}

%\usage{
%.Call(C_isOmp)
%}

\keyword{internal}

