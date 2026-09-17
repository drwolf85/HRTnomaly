.onLoad <- function(lib, pkg) {
       .C(C_setNThreads, n = as.integer(1L))
}

.onAttach <- function(lib, pkg) {
  if (.Call(C_isOmp)) {
    packageStartupMessage(sprintf("Package compiled with openMP (version released in %.2f)",
      .Call(C_openMP_version)))
  }
  else {
    packageStartupMessage("Package compiled without openMP.")
  }
  nn <- setCores()
  if (nn > 0L)
    packageStartupMessage("Use the function setCores() to change the number of CPU cores.")
}

