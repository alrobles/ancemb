# Build the phylogenetic shared-height matrix for all nodes

Returns the matrix \`M\` where \`M\[i,j\]\` is the height of the MRCA of
nodes \`i\` and \`j\`. This is the generalisation of \`ape::vcv\` to all
nodes.

## Usage

``` r
build_all_node_covariance_r(tree)
```

## Arguments

- tree:

  A \`phylo\` object.

## Value

A symmetric matrix with \`n_tips + n_nodes\` rows/columns.
