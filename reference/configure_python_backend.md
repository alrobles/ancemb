# Configure and inspect the optional Python backend

The Python backend is deliberately optional. It provides data I/O and
specialised representation methods while \`ancemb\` retains ownership of
phylogenetic data, Stan inference, posterior summaries, and metrics.

## Usage

``` r
configure_python_backend(
  python = NULL,
  virtualenv = NULL,
  condaenv = NULL,
  required_modules = "numpy"
)

python_backend_status(required_modules = "numpy")
```

## Arguments

- python:

  Optional path to a Python executable.

- virtualenv:

  Optional virtualenv name or path.

- condaenv:

  Optional conda environment name.

- required_modules:

  Python modules to check.

## Value

Invisibly, a backend status list.

## Examples

``` r
if (FALSE) { # \dontrun{
configure_python_backend(virtualenv = "ancemb")
} # }
```
