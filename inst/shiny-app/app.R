# Companion Shiny app for the iPEB package.
# Launched via iPEB::run_app(), or directly with
#   shiny::runApp(system.file("shiny-app", package = "iPEB"))
# Requires the 'shiny' package (a Suggested dependency of iPEB).

library(shiny)
library(iPEB)

ui <- fluidPage(
  titlePanel("iPEB explorer"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("source", "Data",
        c("Bundled example (ipeb_example)" = "example",
          "Upload a CSV" = "upload")),
      conditionalPanel("input.source == 'upload'",
        fileInput("file", "CSV file", accept = ".csv"),
        helpText("Long format: one row per subject-visit.")),
      uiOutput("col_ui"),
      uiOutput("marker_ui"),
      uiOutput("cov_ui"),
      selectInput("objective", "Objective",
        c("Sensitivity" = "sensitivity", "Lead time" = "leadtime",
          "Combined" = "combined")),
      sliderInput("alpha", "Operating specificity (alpha, for fitting)",
        0.5, 0.99, 0.95, 0.01),
      numericInput("window_train", "Training window (months; blank = whole trajectory)",
        value = NA, min = 1),
      selectInput("slope", "Random slope", c("auto", "on", "off")),
      selectInput("innovation", "Innovations", c("auto", "ar1", "iid")),
      selectInput("select", "Feature selection", c("none", "backward")),
      conditionalPanel("input.select != 'none'",
        numericInput("n_markers", "Target panel size (blank = prune by gate)",
          value = NA, min = 1)),
      textInput("specs", "Evaluation specificities (comma-separated)",
        "0.90, 0.95, 0.99"),
      numericInput("window_eval", "Evaluation window (months; blank = same as training)",
        value = NA, min = 1),
      actionButton("run", "Fit and evaluate", class = "btn-primary")
    ),
    mainPanel(
      h4("Fitted model"), verbatimTextOutput("summary"),
      h4("Marker weights"), plotOutput("weights", height = "300px"),
      h4("Evaluation on the test subjects"), tableOutput("eval")
    )
  )
)

server <- function(input, output, session) {

  raw <- reactive({
    if (input$source == "example") {
      data("ipeb_example", package = "iPEB", envir = environment())
      ipeb_example
    } else {
      req(input$file)
      utils::read.csv(input$file$datapath, stringsAsFactors = FALSE)
    }
  })

  output$col_ui <- renderUI({
    cols <- names(raw()); req(cols)
    pick <- function(id, label, guess)
      selectInput(id, label, cols, selected = if (guess %in% cols) guess else cols[1])
    tagList(
      pick("id", "Subject id", "id"),
      pick("case", "Case indicator (1/0)", "case"),
      pick("time", "Visit time", "time"),
      pick("ttd", "Days to diagnosis", "time_to_dx")
    )
  })

  output$marker_ui <- renderUI({
    cols <- names(raw()); req(cols)
    used <- c(input$id, input$case, input$time, input$ttd, "split")
    guess <- setdiff(cols, used)
    numeric_guess <- guess[vapply(guess, function(cc) is.numeric(raw()[[cc]]), logical(1))]
    checkboxGroupInput("markers", "Markers", choices = numeric_guess,
                       selected = head(numeric_guess, 6))
  })

  output$cov_ui <- renderUI({
    cols <- names(raw()); req(cols)
    used <- c(input$id, input$case, input$time, input$ttd, "split", input$markers)
    checkboxGroupInput("covariates", "Covariates (optional)",
                       choices = setdiff(cols, used), selected = character(0))
  })

  fit_res <- eventReactive(input$run, {
    df <- raw()
    validate(need(length(input$markers) >= 1, "Select at least one marker."))
    if ("split" %in% names(df)) {
      tr <- df[df$split == "train", , drop = FALSE]
      te <- df[df$split == "test", , drop = FALSE]
    } else {
      ids <- unique(df[[input$id]])
      tr_ids <- sample(ids, round(0.7 * length(ids)))
      tr <- df[df[[input$id]] %in% tr_ids, , drop = FALSE]
      te <- df[!df[[input$id]] %in% tr_ids, , drop = FALSE]
    }
    specs <- as.numeric(strsplit(input$specs, ",")[[1]])
    specs <- specs[is.finite(specs) & specs > 0 & specs < 1]
    wtr <- if (is.null(input$window_train) || is.na(input$window_train)) Inf else input$window_train
    wev <- if (is.null(input$window_eval) || is.na(input$window_eval)) NULL else input$window_eval
    nm  <- if (is.null(input$n_markers) || is.na(input$n_markers)) NULL else input$n_markers
    fit <- ipeb(tr, markers = input$markers, id = input$id, case = input$case,
                time = input$time, time_to_dx = input$ttd,
                covariates = if (is.null(input$covariates)) character(0) else input$covariates,
                objective = input$objective, alpha = input$alpha, window = wtr,
                slope = input$slope, innovation = input$innovation,
                select = input$select, n_markers = nm)
    list(fit = fit, eval = evaluate(fit, te, specificities = specs, window = wev))
  })

  output$summary <- renderPrint(print(fit_res()$fit))
  output$weights <- renderPlot(plot(fit_res()$fit))
  output$eval <- renderTable(fit_res()$eval, digits = 3)
}

shinyApp(ui, server)
