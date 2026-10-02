# Detect model family from a Stan fit

Uses the presence of \`rate\[1\]\` (BM) or \`sigma2\[1\]\` (OU) in the
posterior draws.

## Usage

``` r
.detect_family(fit)
```

## Arguments

- fit:

  A \`CmdStanMCMC\` object.

## Value

Character \`"bm"\` or \`"ou"\`.
