## Submission

This is a new release: iPEB 0.1.0.

## Test environments

* local: macOS 26.5.2, aarch64-apple-darwin20, R 4.3.1
* win-builder: R release and R Under development (unstable) (x86_64-w64-mingw32)
* macOS builder: R release (aarch64-apple-darwin20)
* GitHub Actions (r-lib/actions, R CMD check):
  - macOS-latest, R release
  - windows-latest, R release
  - ubuntu-latest, R release / R devel / R oldrel-1

## R CMD check results

0 errors | 0 warnings | 1 note.

The single NOTE is the standard "New submission" flag. Where it also lists
"autocorrelated" as possibly misspelled in the DESCRIPTION, this is a
correctly-spelled technical term (autocorrelated residuals).

A locally-produced NOTE, "unable to verify current time", appears only on the
maintainer's machine when the time server is unreachable; it is unrelated to the
package.

## Notes

* A manuscript describing the method is in preparation (matching the DESCRIPTION);
  a reference with DOI will be added on acceptance.
* The example dataset (`ipeb_example`) is entirely synthetic. The real cohort
  used in the accompanying paper (PLCO) is controlled-access via the NCI Cancer
  Data Access System and is not distributed with the package.
* Examples that fit models are wrapped in `\donttest{}` to keep total example
  time small; they run in well under 5 seconds with `--run-donttest`.
* `run_app()` launches an interactive Shiny application and is therefore wrapped
  in `\dontrun{}`; `shiny` is a suggested dependency and its absence is handled.
