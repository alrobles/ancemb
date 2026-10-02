# Assert that a phylo tree has usable branch lengths

Brownian Motion and Ornstein-Uhlenbeck models are defined in terms of
branch lengths; several \`ape\`/\`phytools\` calls fail later with less
actionable messages when \`edge.length\` is missing or malformed.

## Usage

``` r
.assert_branch_lengths(tree)
```

## Arguments

- tree:

  A \`phylo\` object.

## Value

\`tree\`, invisibly.
