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

  )
}


#' merge_data Server Functions
#'
#' @noRd
mod_merge_data_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

  })
}

## To be copied in the UI
# mod_merge_data_ui("merge_data_1")

## To be copied in the server
# mod_merge_data_server("merge_data_1")
