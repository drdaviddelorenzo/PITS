## Resubmission

This is a resubmission. The initial submission was flagged by the incoming
pre-test for "Overall checktime 12 min > 10 min". The vignette
`cdss-cfr-example.Rmd` ran large Monte Carlo simulations (including a full
parameter grid) during rebuilding. The number of simulation replications in the
vignette has been reduced to illustrative sizes (n_sim = 100-200, smaller
parameter grid), with a note advising users to use n_sim >= 1000 for real
analyses. Vignette re-building now completes in under 20 seconds and the full
`R CMD check --as-cran` runs in about one minute locally.

## Submission summary

PITS provides simulation-based statistical power analysis for Interrupted Time
Series (ITS) study designs, and is the R companion to a methods manuscript
currently in preparation.

## Test environments

* local: macOS, R 4.6.0 (R CMD check --as-cran)
* GitHub Actions (r-lib/actions): macOS-latest (release), windows-latest
  (release), ubuntu-latest (devel, release, oldrel-1)
* win-builder: R-devel and R-release
* macOS builder (R-release)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'David de Lorenzo <drdaviddelorenzo@gmail.com>'
  New submission

  This NOTE is expected for a first submission.

## Notes for the reviewer

* The package uses UK (British) English in its documentation and messages
  (e.g. "optimisation", "behaviour", "centred"). Any spell-check flags refer to
  these intentional British spellings.

* The methods manuscript that this package accompanies is not yet published. A
  citation (currently "manuscript in preparation", see inst/CITATION) and a
  reference with a DOI in the Description field will be added in a patch release
  as soon as the paper is accepted.
