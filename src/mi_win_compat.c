#include <stdlib.h>
#include <mimalloc.h>

void* mi_malloc(size_t size) mi_attr_noexcept {
  return malloc(size);
}

void* mi_calloc(size_t count, size_t size) mi_attr_noexcept {
  return calloc(count, size);
}

void mi_free(void* p) mi_attr_noexcept {
  free(p);
}