#' merge_data UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_merge_data_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
        h1("Merge data"),
      ),
      fluidRow(
        # Box for uploading and merging the data
        shinydashboardPlus::box(
          id = ns("box_upload"),
          width = 5,
          solidHeader = TRUE,
          status = "purple",
          title = div(
            id = ns("box_upload_header"),
            "Upload data"
          ),
          fileInput(
            ns("skyline_file"),
            "Upload Skyline CSV file",
            accept = ".csv"
          ),
          fileInput(
            ns("glycounter_files"),
            HTML("Upload GlyCounter <i>OxoSignal</i> text files"),
            accept = ".txt",
            multiple = TRUE
          ),
          actionButton(ns("merge_data"), "Merge data")
        ),
        # Box for displaying and downloading the data
        shinydashboardPlus::box(
          id = ns("box_download"),
          width = 7,
          solidHeader = TRUE,
          status = "purple",
          title = div(
            id = ns("box_download_header"),
            "View and export data"
          ),
          downloadButton(
            ns("download_data"),
            HTML("Download <i>xlsx</i> file")
          ),
          br(),
          DT::dataTableOutput(ns("merged_data"))
        )
      )
    )
  )
}


#' merge_data Server Functions
#'
#' @noRd
mod_merge_data_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # TODO: Checks on validity of the uploaded data.

    skyline_data <- reactive({
      req(input$skyline_file$datapath)
      FALSE  # Call function here
    }) |> bindEvent("merge_data")

    glycounter_data <- reactive({
      req(input$glycounter_files$datapath)
      FALSE  # Call function here
    }) |> bindEvent("merge_data")

    merged_data <- reactive({
      req(skyline_data(), glycounter_data())
      FALSE  # Call function here
    }) |> bindEvent("merge_data")

    # Control status of "Merge data" button
    shinyjs::toggleState("merge_data", is_truthy(merged_data()))

    # Display merged data in table
    output$merged_data <- DT::renderDT({
      req(merged_data())
      DT::datatable(
        data = merged_data(),
        filter = "top",
        options = list(
          scrollX = TRUE,
          pageLength = 6,
          columnDefs = list(list(className = "dt-center", targets = "_all"))
        )
      )
    })
  })
}

## To be copied in the UI
# mod_merge_data_ui("merge_data_1")

## To be copied in the server
# mod_merge_data_server("merge_data_1")
