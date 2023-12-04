#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#





ui <- dashboardPage(
  dashboardHeader(
    title = div(
      img(src = "logo.png", height = "50px"),
      "GRN VTOOLS",
      style = "display: flex; align-items: center; justify-content: center;"
    ),
    titleWidth = 300,
    # Add the Log Out button to the header
    tags$li(
      style = "position: absolute; top: 10px; right: 10px;", class = "dropdown",
      hidden(
        div(id = "logout_button",
            actionButton("logout", "Log Out", icon = icon("sign-out"))
        )
      )
    )
  ),
  
  
  dashboardSidebar(
    useShinyjs(),
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      
      div(id = "register_menu", 
          menuItem("Register", tabName = "register", icon = icon("user-plus"))),
      br(),
      
      div(id = "login_menu",
          menuItem("Log in", tabName = "log_in", icon = icon("sign-in"))),
      br(),
      # Data menu item
      
      div(id = "data_menu",
          hidden(
            div(id = "data_item",
                menuItem("Data", tabName = "data", icon = icon("database"),
                         menuSubItem("Upload Data", tabName = "upload_data"),
                         menuSubItem("Normalize Data", tabName = "normalize"),
                         menuSubItem("Differential Expression", tabName = "eda"),
                         menuSubItem("Expression based clustering", tabName = "ebc")
                )
            )
          )
      ),
      br(),
      # Tools menu item
      div(id = "tools_menu",
          hidden(
            div(id = "tools_item",
                menuItem("Tools", tabName = "tools", icon = icon("wrench"))
            )
          ),
          uiOutput("tool_input")
      ), #ADDED
      br(),
      # Generate Network menu item
      div(id = "network_menu",
          hidden(
            div(id = "network_item",
                menuItem("Generate Network", tabName = "generate_network", icon = icon("sitemap"))
            )
          )
      )
      ,
      br(),
      div(id = "download_menu",
          hidden(
            div(id = "download_item",
                menuItem("Download", tabName = "download", icon = icon("download"))
            )
          )
      )
      
      
      
      
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
      )
      ,
      
      tabItem(tabName = "normalize",
              h2("Normalize Data"),
              actionButton("norm_button", label = "Normalize", 
                           icon = icon("sitemap"), 
                           style = "background-color: purple; color: white"),
              box(
                title = "Normalized Data",
                status = "primary",
                solidHeader = TRUE,
                width = 12,
                fluidRow(
                  column(width = 6,
                         plotOutput("norm_plot")
                  ),
                  column(width = 6,
                         tableOutput("normalized_table")
                  )
                )
              )
      ),
      
      tabItem(tabName = "eda",
              h2("Exploratory Expression"),
              actionButton("de_button", label = "Differential Expression", 
                           icon = icon("sitemap"), 
                           style = "background-color: purple; color: white"),
              box(
                title = "Differential Expression",
                status = "primary",
                solidHeader = TRUE,
                width = 12,
                plotOutput("de_plot")
              )
      ),
      tabItem(tabName = "ebc",
              h2("Expression based clustering"),
              
              
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
                actionButton("change_color", "Change Color"),  # Add the Change Color button here
                actionButton("revert_color", "Revert Color", class = "btn-primary", style = "display: none;")
              )
              
              
              
              
      ),
      
      tabItem(tabName = "download",
              h2("Download Data"),
              downloadButton("download_data", "Download Abiotic Stresses Data")
      ),
      
      
      tabItem(tabName = "log_in",  ##changed name here
              h2("Log in"),
              
              
              div(id = "login_button", 
                  
                  textInput("username", "Username"),
                  passwordInput("password", "Password"),
                  br(),
                  actionButton("login", "Log In") )
              
      ),
      tabItem(tabName = "register",
              h2("Register an Account"),
              
              div(id = "register_button",
                  tags$p("Don't have an account? Register below."), 
                  tags$p("Make sure your username is 8 characters long, your password is 8 characters long at least, has 2 special characters at least, and 1 capital letter at least.", style="color:#FF0000"),
                  
                  textInput("new_username", "New Username"),
                  passwordInput("new_password", "New Password"),
                  passwordInput("confirm_password", "Confirm Password"),
                  actionButton("register", "Register")
              )  
              
              
      ))
  )
)

# ... (Previous code)







