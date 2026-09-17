#include "myomp.h"
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Lapack.h>
#include <R_ext/Boolean.h>

/**
 * @wrapper C_getNCores
 * @brief Get the max number of CPU cores
 * @param n Pointer to the number of CPU cores
 */
extern void getNCores(int *n) {
  #if __VOPENMP
    *n = omp_get_num_procs();
  #else
    *n = 1;
  #endif
}

/**
 * @wrapper C_getNThreads
 * @brief Get the number of threads to use
 * @param n Pointer to the number of threads
 */
extern void getNThreads(int *n) {
  #if __VOPENMP
    #pragma omp parallel default(shared)
    {
      #pragma omp masked
        *n = omp_get_num_threads();
    }
  #else
    *n = 0;
  #endif
}

/**
 * @wrapper C_setNThreads
 * @brief Set the number of threads to use
 * @param n Pointer to the number of threads
 */
extern void setNThreads(int *n) {
  #if __VOPENMP
    if (omp_get_num_procs() < *n) {
      *n = omp_get_num_procs();
      omp_set_num_threads(*n);
    }
    else if(*n > 0) {
      omp_set_num_threads(*n);
    }
    else {
      omp_set_num_threads(1);
      *n = 1;
    }
  #else
    *n = 1;
  #endif
}

/**
 * @wrapper C_isOmp
 * @brief Check if the package allows for parallel computing via OpenMP
 * @return SEXP
 */
extern SEXP isOmp(void) {
  SEXP ans;
  PROTECT(ans = allocVector(LGLSXP, 1));
  #if __VOPENMP
    LOGICAL(ans)[0] = TRUE;
  #else
    LOGICAL(ans)[0] = FALSE;
  #endif
  UNPROTECT(1);
  return ans;
}

/**
 * @wrapper C_openMP_version
 * @brief Retrieve the OpenMP version
 * @return SEXP
 */
extern SEXP openMP_version(void) {
  SEXP ans;
  PROTECT(ans = allocVector(REALSXP, 1));
  #if __VOPENMP
    REAL(ans)[0] = (double) _OPENMP * 1e-2;
  #else
    REAL(ans)[0] = NA_REAL;
  #endif
  UNPROTECT(1);
  return ans;
}
