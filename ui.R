ui <- dashboardPage(
  
  # Define the dashboard header
  dashboardHeader(title = div(
    ##img(src = normalizePath("C:/Users/user/Desktop/GRNN/GRN/GRNNEW/logo.jpeg"), height = "30px"),
    img(src="logo.png", height = "50px"),
    "GRN VTOOLS",
    style = "display: flex; align-items: center; justify-content: center;"
  ), 
  titleWidth = 300,
  # Add the Log Out button to the header
  tags$li(
    style = "position: absolute; top: 10px; right: 10px;", class = "dropdown",
    actionButton("logout", "Log Out", icon = icon("sign-out"))
  )
  ),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      
      menuItem("Data", tabName = "data", icon = icon("database"),
               menuSubItem("Upload Data", tabName = "upload_data"),
               menuSubItem("Normalize Data", tabName = "normalize"),
               menuSubItem("Exploratory Expression", tabName = "eda")
      ),
      
      menuItem("Tools", tabName = "tools", icon = icon("wrench")),
      menuItem("Generate Network", tabName = "generate_network", icon = icon("sitemap")),
      menuItem("Download", tabName = "download", icon = icon("download")),
      menuItem("Log in", tabName = "log_in", icon = icon("sign-in")),
      menuItem("Register", tabName = "register", icon = icon("user-plus"))
    )
  ),
  
  dashboardBody(
    # Use shinyjs to hide/show UI elements dynamically
    useShinyjs(),
    tabItems(
      tabItem(tabName = "home",
              style = "display: flex; flex-direction: column; align-items: center; justify-content: center;",
              h2(style = "text-align: center;","Welcome to GRN VTOOLS"),
              img(src="Home.png", height = "500px", style = "display: block; margin: 0 auto;"),
              p("Join us on this exciting journey as we provide you with the tools you need to unravel the mysteries of gene regulatory networks. Experience the power and versatility of GRN VTOOLS and unlock a deeper understanding of gene expression and regulation.", style="color:#800080;font-weight: bold; text-align: center;"),
              
              # Add more content here
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
              p("Data normalization content goes here.")
              # Add content for data normalization
      ),
      
      tabItem(tabName = "eda",
              h2("Exploratory Expression"),
              p("Exploratory data analysis content goes here.")
              # Add content for exploratory data analysis
      ),
      
      tabItem(tabName = "tools",
              h2("Tools"),
              p("Choose a network generation tool:"),
              selectInput("tool_choice", label = "Choose Network Generation Tool:", 
                          choices = list("SeqNet"), 
                          selected = "SeqNet")),
      
      tabItem(tabName = "generate_network",
              h2("Generate Network"),
              p("Network generation content goes here."),
              actionButton("network_button", label = "Generate Network", 
                           icon = icon("sitemap"), 
                           style = "background-color: purple; color: white"),
              box(
                title = "Gene Network",
                status = "primary",
                solidHeader = TRUE,
                width = 12,
                plotOutput("network_plot")
              )
      ),
      
      tabItem(tabName = "log_in",
              h2("Log in"),
              conditionalPanel(
                condition = "!Cookies.get('loggedIn')",
                div(id = "login_button",
                    textInput("username", "Username"),
                    passwordInput("password", "Password"),
                    br(),
                    actionButton("login", "Log In")
                )
              ),
              conditionalPanel(
                condition = "Cookies.get('loggedIn')",
                div(id = "logout_button1",
                    actionButton("logout", "Log Out")
                )  
                
              )
      ),
      
      tabItem(tabName = "register",
              h2("Register an Account"),
              conditionalPanel(
                condition = "!Cookies.get('loggedIn')",
                div(id = "register_button",
                    tags$p("Don't have an account? Register below."), 
                    textInput("new_username", "New Username"),
                    passwordInput("new_password", "New Password"),
                    passwordInput("confirm_password", "Confirm Password"),
                    actionButton("register", "Register")
                )
              ),
              conditionalPanel(
                condition = "Cookies.get('loggedIn')",
                div(id = "logout_button2",
                    actionButton("logout", "Log Out")
                )  
                
              )
      )
    )
  )
)
