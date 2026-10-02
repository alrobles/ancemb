# Reconstruct ancestral embeddings with \`phytools::fastAnc\`

Convenience wrapper around \`phytools::fastAnc\` that returns the full
matrix of ancestral node estimates (one column per embedding dimension).

## Usage

``` r
reconstruct_ancestral_ml_r(tree, tip_embeddings)
```

## Arguments

- tree:

  A \`phylo\` object.

- tip_embeddings:

  Numeric matrix (\`n_tips x D\`).

## Value

Numeric matrix (\`n_nodes x D\`). Row 1 is the root estimate.
