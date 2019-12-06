## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE,
  fig.width = 7
)

## ---- echo=FALSE---------------------------------------------------------
isINLA <- requireNamespace('INLA', quietly = TRUE)
