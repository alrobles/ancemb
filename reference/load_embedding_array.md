# Load a NumPy array through the optional backend

Load a NumPy array through the optional backend

## Usage

``` r
load_embedding_array(path, key = NULL)
```

## Arguments

- path:

  Path to \`.npy\` or \`.npz\` file.

- key:

  Name of the array inside an \`.npz\` archive.

## Value

A numeric matrix or array.

## Examples

``` r
if (FALSE) { # \dontrun{
load_embedding_array("embeddings.npy")
} # }
```
