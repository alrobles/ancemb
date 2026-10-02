# Simulate BM evolution of an embedding along a phylogeny

Simulate BM evolution of an embedding along a phylogeny

## Usage

``` r
simulate_bm_embedding_r(tree, z_root, rates)
```

## Arguments

- tree:

  A \`phylo\` object.

- z_root:

  Numeric vector, root state (length \`D\`).

- rates:

  Numeric vector, BM rate per dimension (length \`D\`).

## Value

List with \`tip_values\` (\`n_tips\` x \`D\`), \`node_values\`
(\`n_nodes\` x \`D\`), \`all_values\` (\`n_total\` x \`D\`) and
\`root_value\`.
