# Internal helper that compiles and samples a Stan model

Internal helper that compiles and samples a Stan model

## Usage

``` r
.fit_embedding_model(
  data,
  family,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  ...
)
```

## Arguments

- data:

  Stan data list.

- family:

  Either \`"bm"\` or \`"ou"\`.

- chains, iter_warmup, iter_sampling:

  Passed to \`\$sample()\`.

- ...:

  Additional arguments passed to \`\$sample()\`.
