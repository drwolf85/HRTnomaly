library(HRTnomaly)

expected_exports <- c(
  "bayesHRT",
  "bayeswise",
  "bootHRT",
  "cellwise",
  "class_check",
  "dif",
  "fuzzyHRT",
  "gif",
  "pif",
  "setCores"
)

namespace <- asNamespace("HRTnomaly")
exports <- getNamespaceExports(namespace)

stopifnot(
  identical(sort(exports), sort(expected_exports)),
  !any(startsWith(exports, "C_")),
  !any(c(".onLoad", ".onAttach", "print.checkwise") %in% exports),
  is.function(getS3method("print", "checkwise", optional = TRUE))
)

registered <- getDLLRegisteredRoutines(getLoadedDLLs()[["HRTnomaly"]])
native_bindings <- paste0(
  "C_",
  c(names(registered[[".C"]]), names(registered[[".Call"]]))
)

stopifnot(
  length(registered[[".C"]]) == 18L,
  length(registered[[".Call"]]) == 3L,
  all(native_bindings %in% ls(namespace, all.names = TRUE)),
  !any(native_bindings %in% exports)
)
