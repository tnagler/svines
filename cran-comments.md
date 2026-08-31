This release adds support for discrete variables and compatibility with
rvinecopulib 1.0.0. It also updates the package to require C++17.

## Test environments

Local checks on Linux Mint 21.3 with R 4.6.1:

* `R CMD check --no-manual` with the CRAN release of rvinecopulib
  (0.7.3.1.0)
* `R CMD check --no-manual` with the development snapshot for the upcoming
  rvinecopulib 1.0.0 release (`final-sync`, commit 98c43e0; package version
  1.0.0.1.0)
* the complete testthat suite with Eigen assertions enabled against the same
  development snapshot (89 tests passed)

GitHub Actions checks using CRAN dependencies:

* Ubuntu: R release, R devel, and R 4.3
* Windows: R release
* macOS: R release

All five GitHub Actions `R CMD check` jobs completed successfully. The
separate test-coverage and pkgdown workflows also completed successfully.

## R CMD check results

0 errors | 0 warnings | 0 notes
