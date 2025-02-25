#oncoexpr

Hi, this is a project for Oncology RNAseq Analysis.

## Installation

You can install OncoExpr directly from GitHub using the following R command:

```R
remotes::install_github('atgenomix/oncoexpr')
```

```R
library(oncoexpr)
```

Note: A spark connection is necessary for 'oncoExprAppSpark' shiny app.
```R
oncoExprAppSpark()
```
