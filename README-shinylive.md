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

## Dependency availability (read this if the export step fails)

`shinylive::export()` installs the app's R dependencies as WebAssembly binaries
from <https://repo.r-wasm.org>, which mirrors CRAN. The app needs `shiny`,
`iPEB`, and `nlme`; `shiny` and `nlme` have been there for a long time.

`iPEB` reached CRAN on 2026-09-17, and the webR repository rebuilds from CRAN on
its own schedule, so the binary may not exist for a few days after publication.
If the export step fails with a message about `iPEB` not being available, that
is the reason: wait a few days and re-run the workflow from the Actions tab.
Nothing needs to change in the code. iPEB is a pure R package
(`NeedsCompilation: no`), so it should build for WebAssembly without trouble.

To check availability before re-running, search for `iPEB` in
<https://repo.r-wasm.org/bin/emscripten/contrib/4.4/PACKAGES>.

## Testing the export locally

```r
install.packages("shinylive")
shinylive::export(appdir = "inst/shiny-app", destdir = "_site")
httpuv::runStaticServer("_site")   # opens the exported app in a browser
```

Delete `_site/` afterwards; it is build output and is not tracked.
