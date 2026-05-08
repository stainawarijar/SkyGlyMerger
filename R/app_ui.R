#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    shinydashboardPlus::dashboardPage(
      title = "SkyGlyMerger", skin = "purple",
      header = shinydashboard::dashboardHeader(
        title = "SkyGlyMerger"
      ),
      sidebar = shinydashboardPlus::dashboardSidebar(
        shinydashboard::sidebarMenu(
          id = "tabs",
          shinydashboard::menuItem(
            text = HTML("&nbspMerge data"),
            tabName = "merge_data", icon = icon("file-import")
          )
        )
      ),
      body = shinydashboard::dashboardBody(
        # Code to keep title and icons visible when collapsing the sidebar
        tags$style(
          '
          @media (min-width: 768px){
            .sidebar-mini.sidebar-collapse .main-header .logo {
                width: 230px;
            }
            .sidebar-mini.sidebar-collapse .main-header .navbar {
                margin-left: 230px;
            }
          }
          '
        ),
        shinydashboard::tabItems(
          shinydashboard::tabItem(
            "merge_data", mod_merge_data_ui("merge_data_1")
          )
        )
      )
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "SkyGlyMerger"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
