test_that("Python backend configuration validates mutually exclusive environments", {
  skip_if_not_installed("reticulate")
  expect_error(
    configure_python_backend(
      python = "python",
      virtualenv = "environment",
      required_modules = "numpy"
    ),
    "only one"
  )
})

test_that("NPY loading uses the R-native reader when available", {
  skip_if_not_installed("RcppCNPy")
  path <- tempfile(fileext = ".npy")
  on.exit(unlink(path), add = TRUE)
  expected <- matrix(as.numeric(seq_len(6)), nrow = 2)
  RcppCNPy::npySave(path, expected)
  expect_equal(load_embedding_array(path), expected)
})
