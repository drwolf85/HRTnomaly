#ifndef HRTNOMALY_ALLOCATOR_H

    #define HRTNOMALY_ALLOCATOR_H

    #ifndef MI_CRAN_COMPLIANT
        #define MI_CRAN_COMPLIANT 0
    #endif

    #if MI_CRAN_COMPLIANT

        #include <mimalloc.h>

    #else

        #include <stdlib.h>

        #define mi_malloc malloc
        #define mi_calloc calloc
        #define mi_free free

    #endif

#endif
