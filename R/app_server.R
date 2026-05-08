#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # Increase max. uploaded file size from 5MB to 200MB
  options(shiny.maxRequestSize=200*1024^2)  # Change as necessary

  # Extract results from tabs
  results_merge_data <- mod_merge_data_server("merge_data_1")
}
