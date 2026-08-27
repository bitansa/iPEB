#' Launch the iPEB Shiny explorer
#'
#' Opens the companion Shiny application bundled with the package. The app lets
#' you load the built-in synthetic example (or upload a long-format CSV), choose
#' the objective, operating specificity, windowing, layer options and feature
#' selection, then fit \code{\link{ipeb}} and view the marker weights and the
#' held-out evaluation interactively.
#'
#' The app requires the \pkg{shiny} package (a suggested dependency); install it
#' with \code{install.packages("shiny")} if needed. A hosted, browser-based
#' version for users without R will be linked from the package README once deployed.
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}} (for
#'   example \code{launch.browser} or \code{port}).
#'
#' @return Invisibly \code{NULL}; called for the side effect of launching the app.
#'
#' @examples
#' \dontrun{
#' run_app()
#' }
#'
#' @seealso \code{\link{ipeb}}
#' @export
run_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The 'shiny' package is required to run the app. ",
         "Install it with install.packages(\"shiny\").", call. = FALSE)
  }
  app_dir <- system.file("shiny-app", package = "iPEB")
  if (app_dir == "" || !file.exists(file.path(app_dir, "app.R"))) {
    stop("Could not locate the bundled Shiny app; try reinstalling iPEB.",
         call. = FALSE)
  }
  shiny::runApp(app_dir, ...)
  invisible(NULL)
}
