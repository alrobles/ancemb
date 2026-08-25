#' Configure and inspect the optional Python backend
#'
#' The Python backend is deliberately optional. It provides data I/O and
#' specialised representation methods while `ancemb` retains ownership of
#' phylogenetic data, Stan inference, posterior summaries, and metrics.
#'
#' @param python Optional path to a Python executable.
#' @param virtualenv Optional virtualenv name or path.
#' @param condaenv Optional conda environment name.
#' @param required_modules Python modules to check.
#' @return Invisibly, a backend status list.
#' @export
#' @examples
#' \dontrun{
#' configure_python_backend(virtualenv = "ancemb")
#' }
configure_python_backend <- function(python = NULL, virtualenv = NULL,
                                      condaenv = NULL,
                                      required_modules = "numpy") {
  checkmate::assert_string(python, null.ok = TRUE)
  checkmate::assert_string(virtualenv, null.ok = TRUE)
  checkmate::assert_string(condaenv, null.ok = TRUE)
  checkmate::assert_character(required_modules, min.len = 1L,
                              any.missing = FALSE)
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Install the optional 'reticulate' package to configure Python.")
  }
  if (sum(!vapply(list(python, virtualenv, condaenv), is.null, logical(1))) > 1L) {
    stop("Specify only one of 'python', 'virtualenv', or 'condaenv'.")
  }
  if (!is.null(python)) {
    reticulate::use_python(python, required = TRUE)
  } else if (!is.null(virtualenv)) {
    reticulate::use_virtualenv(virtualenv, required = TRUE)
  } else if (!is.null(condaenv)) {
    reticulate::use_condaenv(condaenv, required = TRUE)
  }
  missing <- required_modules[
    !vapply(required_modules, reticulate::py_module_available, logical(1))
  ]
  status <- list(
    python = reticulate::py_config()$python,
    modules = required_modules,
    missing_modules = missing,
    available = !length(missing)
  )
  class(status) <- "ancemb_python_status"
  if (length(missing)) {
    warning("Python backend missing modules: ", paste(missing, collapse = ", "))
  }
  invisible(status)
}

#' @rdname configure_python_backend
#' @export
python_backend_status <- function(required_modules = "numpy") {
  configure_python_backend(required_modules = required_modules)
}

#' Load a NumPy array through the optional backend
#'
#' @param path Path to `.npy` or `.npz` file.
#' @param key Name of the array inside an `.npz` archive.
#' @return A numeric matrix or array.
#' @export
#' @examples
#' \dontrun{
#' load_embedding_array("embeddings.npy")
#' }
load_embedding_array <- function(path, key = NULL) {
  checkmate::assert_file_exists(path)
  checkmate::assert_string(key, null.ok = TRUE)
  if (grepl("\\.npy$", path, ignore.case = TRUE) &&
      requireNamespace("RcppCNPy", quietly = TRUE)) {
    return(RcppCNPy::npyLoad(path))
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Loading this file requires 'reticulate' or 'RcppCNPy'.")
  }
  np <- reticulate::import("numpy", delay_load = FALSE)
  loaded <- np$load(path, allow_pickle = FALSE)
  if (grepl("\\.npz$", path, ignore.case = TRUE)) {
    available <- as.character(reticulate::py_to_r(loaded$files))
    if (is.null(key)) {
      if (length(available) != 1L) {
        stop("An NPZ with multiple arrays requires 'key'.")
      }
      key <- available[[1L]]
    }
    if (!key %in% available) {
      stop("Key '", key, "' is not present in ", path, ".")
    }
    return(reticulate::py_to_r(loaded[[key]]))
  }
  reticulate::py_to_r(loaded)
}
