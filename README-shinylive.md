# Hosting the iPEB explorer with shinylive

The workflow in `.github/workflows/shinylive.yaml` compiles `inst/shiny-app/` to
a static [shinylive](https://posit-dev.github.io/r-shinylive/) site and publishes
it to GitHub Pages at <https://bitansa.github.io/iPEB/>.

R runs in the visitor's browser through WebAssembly (webR). Nothing is sent to a
server, so a clinician can load a CSV of patient biomarker values without that
file leaving their computer. There is no server to pay for or maintain, and no
usage limit.

## One-time repository setup

1. On GitHub, open **Settings -> Pages**.
2. Under **Build and deployment -> Source**, choose **GitHub Actions**
   (not "Deploy from a branch").
3. Push to `main`, or run the workflow manually from the **Actions** tab
   (`shinylive` -> **Run workflow**).

The first successful run publishes the site; later pushes that touch
`inst/shiny-app/`, `R/`, or `DESCRIPTION` refresh it automatically.

## Why the app is assembled instead of installed

`shinylive::export()` can only install packages that exist as WebAssembly
binaries at <https://repo.r-wasm.org>, which mirrors CRAN on its own schedule.
iPEB reached CRAN on 2026-09-17, and its WebAssembly binary did not appear in
that repository in the days afterwards, so `library(iPEB)` failed at runtime in
the browser with "there is no package called 'iPEB'" even though the export and
deploy steps both succeeded.

The workflow therefore assembles a **self-contained** app rather than relying on
that repository:

- the package's `R/*.R` sources are copied to `_app/ipeb-R/` and sourced
  directly by `app.R` (replacing `library(iPEB)`);
- `data/ipeb_example.rda` is copied next to the app and loaded with `load()`
  (replacing the `data()` call);
- `R/run_app.R` is skipped, since it exists only to launch this same app.

This works because iPEB is pure R (`NeedsCompilation: no`) and every call it
makes into nlme is namespace-qualified (`nlme::lme`, and so on), so no attach
step is needed and nlme ships with R itself. The only package shinylive has to
fetch is `shiny`.

The sources are copied at build time from `R/`, so there is no duplicated code
in the repository: `R/` remains the single source of truth, and the deployed app
updates automatically whenever `R/`, `data/`, or the app changes.

If iPEB does eventually appear in the webR repository, this setup keeps working
and needs no change.

## Testing the export locally

Run the "Assemble self-contained app" step from
`.github/workflows/shinylive.yaml` in an R session at the package root, then:

```r
install.packages("shinylive")
shinylive::export(appdir = "_app", destdir = "_site")
httpuv::runStaticServer("_site")   # opens the exported app in a browser
```

Delete `_app/` and `_site/` afterwards; both are build output and neither is
tracked.
