#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(shiny)
library(shinydashboard)
library(SeqNet)
library(igraph)
library(readr)


ui <- dashboardPage(
  dashboardHeader(
    title = div(
      img(src = "logo.png", height = "50px"),
      "GRN VTOOLS",
      style = "display: flex; align-items: center; justify-content: center;"
    ),
    titleWidth = 300,
    tags$li(
      style = "position: absolute; top = 10px; right = 10px;", class = "dropdown",
      actionButton("logout", "Log Out", icon = icon("sign-out"))
    )
  ),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Data", tabName = "data", icon = icon("database"),
               menuSubItem("Upload Data", tabName = "upload_data"),
               menuSubItem("Normalize Data", tabName = "normalize"),
               menuSubItem("Exploratory Expression", tabName = "eda"),
               menuSubItem("Expression based clustering", tabName = "ebc")
      ),
      menuItem("Tools", tabName = "tools", icon = icon("wrench"),
               uiOutput("tool_input")
      ),
      menuItem("Generate Network", tabName = "generate_network", icon = icon("sitemap")),
      menuItem("Download", tabName = "download", icon = icon("download")),
      menuItem("Log in", tabName = "log_in", icon = icon("sign-in")),
      menuItem("Register", tabName = "register", icon = icon("user-plus"))
    )
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "home",
              h2("Welcome to GRN VTOOLS"),
              img(src = "home.jpg", height = "500px", style = "display: block; margin: 0 auto;"),
              p("Join us on this exciting journey as we provide you with the tools you need to unravel the mysteries of gene regulatory networks. Experience the power and versatility of GRN VTOOLS and unlock a deeper understanding of gene expression and regulation.", style = "color:#800080;font-weight: bold; text-align: center;")
      ),
      
      tabItem(tabName = "upload_data",
              h2("Upload Data"),
              p("Upload your data files here:"),
              radioButtons("num_files", label = "Number of Files:",
                           choices = list("One File" = "one", "Two Files" = "two"), 
                           selected = "one"),
              conditionalPanel(
                condition = "input.num_files == 'one'",
                fileInput("data_file1", "Upload Data File", accept = c(".csv"))
              ),
              conditionalPanel(
                condition = "input.num_files == 'two'",
                fileInput("data_file1", "Upload Data File 1", accept = c(".csv")),
                fileInput("data_file2", "Upload Data File 2", accept = c(".csv"))
              )
      ),
      
      tabItem(tabName = "normalize",
              h2("Normalize Data"),
              actionButton("norm_button", label = "Normalize", 
                           icon = icon("sitemap"), 
                           style = "background-color: purple; color: white")
              
              ,
              box(
                
                title = "normalized data",
                status = "primary",
                solidHeader = TRUE,
                width = 20,
                
                fluidRow(
                  column(width = 6,  # Specify the width for the first cell (e.g., 50%)
                         plotOutput("norm_plot")  # Plot output for the first cell
                  ),
                  column(width = 6,  # Specify the width for the second cell (e.g., 50%)
                         tableOutput("normalized_table")  # Table output for the second cell
                  )
                )
                
              )
              
      ),
      
      tabItem(tabName = "eda",
              h2("Exploratory Expression"),
              conditionalPanel(
                condition = "input.num_files == 'one'",
                fileInput("data_file3", "Upload Data File", 
                          accept = c(".csv", ".tsv", ".txt"))
              ),
              
              
              
              
              
              actionButton("de_button", label = "differential Expression", 
                           icon = icon("sitemap"), 
                           style = "background-color: purple; color: white")
              
              ,
              box(
                
                title = "differential Expression",
                status = "primary",
                solidHeader = TRUE,
                width = 20,
                
                fluidRow(
                  splitLayout(cellWidths = c("100%"),plotOutput("de_plot")
                              
                              
                  )
                  
                )
              )
      ),
      tabItem(tabName = "ebc",
              h2("Expression based clustering"),
              
              conditionalPanel(
                condition = "input.num_files == 'one'",
                fileInput("data_file3", "Upload Data File", 
                          accept = c(".csv", ".tsv", ".txt"))
              ),
              actionButton("ebc_button", label = "Expression based clustering", 
                           icon = icon("sitemap"), 
                           style = "background-color: purple; color: white")
              
              ,
              box(
                
                title = "Expression based clustering",
                status = "primary",
                solidHeader = TRUE,
                width = 20,
                
                fluidRow(
                  splitLayout(cellWidths = c("100%"),plotOutput("ebc_plot")
                              
                              
                  )
                  
                )
              )
      ),
      
      tabItem(tabName = "tools",
              h2("Tools"),
              p("Choose a network generation tool:"),
              uiOutput("network_tool_boxes")
      ),
      
      tabItem(tabName = "generate_network",
              h2("Generate Network"),
            
              fluidRow(
                column(6, uiOutput("diane_network_box")),
                column(6, uiOutput("seqnet_network_box"))
              ),
              actionButton("generate_button", "Generate Network"),
              
              conditionalPanel(
                "input.tool_choice.includes('seqnet')",
                actionButton("get_sigma_button", "Get Sigma"),  # Add the Get Sigma button
                actionButton("change_color", "Change Color")  # Add the Change Color button here
              )
              
              
              
              
      ),
      
      tabItem(tabName = "log_in",
              h2("Log in"),
              div(
                id = "login_button",
                textInput("username", "Username"),
                passwordInput("password", "Password"),
                br(),
                actionButton("login", "Log In")
              )
      ),
      
      tabItem(tabName = "register",
              h2("Register an Account"),
              div(
                id = "register_button",
                tags$p("Don't have an account? Register below."), 
                textInput("new_username", "New Username"),
                passwordInput("new_password", "New Password"),
                passwordInput("confirm_password", "Confirm Password"),
                actionButton("register", "Register")
              )
      )
    )
  )
)

# ... (Previous code)

