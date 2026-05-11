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
          numericInput(
            ns("mz_tolerance_ppm"),
            HTML(
              "Tolerance around theoretical <i>m/z</i> values (ppm)
              in Skyline data"
            ),
            value = 10, min = 1, step = 1
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
      load_skyline_data(input$skyline_file$datapath)
    }) |> bindEvent(input$merge_data)

    glycounter_data <- reactive({
      req(input$glycounter_files$datapath)
      # Extract filenames of OxoSignal files (original and in memory)
      original_names <- input$glycounter_files$name
      oxosignal_indices <- grepl("_OxoSignal\\.txt$", original_names)
      oxosignal_files <- input$glycounter_files$datapath[oxosignal_indices]
      oxosignal_names <- original_names[oxosignal_indices]
      # Process the OxoSignal files
      # TODO: Check against zero OxoSignal files.
      # TODO: Allow for uploading zip file.
      load_glycounter_data(
        setNames(as.list(oxosignal_files), oxosignal_names)
      )
    }) |> bindEvent(input$merge_data)

    merged_data <- reactive({
      req(skyline_data(), glycounter_data())
      fragment_cols <- extract_fragment_cols(glycounter_data())
      skyline_prepped <- prepare_skyline_data(
        skyline_data(), input$mz_tolerance_ppm
      )
      glycounter_candidates <- extract_glycounter_candidates(
        skyline_prepped, glycounter_data()
      )
      glycounter_summary <- summarize_glycounter_data(
        glycounter_candidates, fragment_cols
      )
      merge_data(skyline_prepped, glycounter_summary)
    })


    # Display merged data in table
    output$merged_data <- DT::renderDT({
      req(merged_data())
      DT::datatable(
        data = merged_data(),
        filter = "top",
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          columnDefs = list(list(className = "dt-center", targets = "_all"))
        )
      )
    })

    # Download data in xlsx format.
    output$download_data <- downloadHandler(
      filename = function() {
        paste0(
          format(Sys.Date(), "%Y%m%d"), "_", format(Sys.time(), "%H%M"),
          "_merged_data.xlsx"
        )
      },
      content = function(file) {
        writexl::write_xlsx(merged_data(), path = file)
      }
    )


    # Control status of buttons.
    observe({
      # Button for merging data
      shinyjs::toggleState(
        id = "merge_data",
        condition = all(
          !is.null(input$skyline_file$datapath),
          !is.null(input$glycounter_files$datapath)
        )
      )
      # Button for downloading data
      shinyjs::toggleState("download_data", is_truthy(merged_data()))
    })


    return(list(
      merged_data = merged_data,
      mz_tolerance_ppm = reactive(input$mz_tolerance_ppm)
    ))
  })
}

## To be copied in the UI
# mod_merge_data_ui("merge_data_1")

## To be copied in the server
# mod_merge_data_server("merge_data_1")
