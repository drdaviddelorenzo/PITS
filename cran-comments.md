## Resubmission

This is a resubmission addressing a CRAN reviewer comment (Benjamin Altmann):
functions must not write to the user's home filespace (including the working
directory) by default.

* `export_results()`: the `dir` argument no longer has a default; the user must
  supply an output directory. The function errors with an informative message if
  `dir` is not given.
* `run_its_power()`: `output_dir` now defaults to `NULL`; when
  `save_output = TRUE` the user must supply a directory, otherwise the function
  errors. With the default `save_output = FALSE`, nothing is written.

No function now writes to the working directory or home filespace by default.
Examples, tests and vignettes that write files use `tempdir()`.

### Earlier resubmission (check time)

A previous resubmission reduced the Monte Carlo replications in the
`cdss-cfr-example.Rmd` vignette (n_sim = 100-200, smaller parameter grid) to fix
an "Overall checktime > 10 min" pre-test NOTE. Vignette re-building now completes
in under 20 seconds; the full `R CMD check --as-cran` runs in about one minute
locally.

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
