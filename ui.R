

dashboardPage(
  
  # Define the dashboard header
  dashboardHeader(title = div(
    
    img(src="logo.png", height = "30px"),
    "GRN VTOOLS",
    style = "display: flex; align-items: center; justify-content: center;"
  ), 
  titleWidth = 300,
  
  # Add the Log Out button to the header
  tags$li(
    style = "position: absolute; top: 10px; right: 10px;", class = "dropdown",
    actionButton("logout", "Log Out" , icon = icon("sign-out"))
  )
  ),
  
  
  # Define the sidebar
  dashboardSidebar(
    
    # Use shinyjs to hide/show UI elements dynamically
    useShinyjs(),
    # Sidebar menu
    sidebarMenu(
      # Home menu item
      menuItem("Home", tabName = "home", icon = icon("home")),
      
      div(id = "register_menu", 
      menuItem("Register", tabName = "Register", icon = icon("user-plus"))),
      br(),
      
      div(id = "login_menu",
      menuItem("Log in", tabName = "Log_in", icon = icon("sign-in"))),
      br(),
      # Data menu item
      div(id = "data_menu",
      menuItem("Data", tabName = "data", icon = icon("database"),
               # Normalize Data submenu item
               menuSubItem("Upload Data", tabName = "upload_data"),
               menuSubItem("Normalize Data", tabName = "normalize"),
               menuSubItem("Differential Expression", tabName = "de"),
               menuSubItem("Expression based clustering", tabName = "ebc")
      )),
      br(),
      # Tools menu item
      div(id = "tools_menu",
      menuItem("Tools", tabName = "tools", icon = icon("wrench")),  
      uiOutput("tool_input")), #ADDED
      br(),
      # Generate Network menu item
      div(id = "network_menu",
      menuItem("Generate Network", tabName = "generate_network", icon = icon("sitemap"))),
      br(),
      div(id = "download_menu",
      menuItem("Download", tabName = "download", icon = icon("download")))
    
      
      
      
    )
  ),
  
  # Define the main body
  dashboardBody(
    # Tab items
    tabItems(
      # Home tab
      tabItem(tabName = "home",
              h2("Welcome to GRN VTOOLS"),
              img(src="Home.png", height = "500px", style = "display: block; margin: 0 auto;"),
             
              p("Join us on this exciting journey as we provide you with the tools you need to unravel the mysteries of gene regulatory networks. Experience the power and versatility of GRN VTOOLS and unlock a deeper understanding of gene expression and regulation.", style="color:#800080;font-weight: bold; text-align: center;"),
              
             
      ),
      
      # Data tab
      tabItem(tabName = "upload_data",
              h2("Upload Data"),
              p("Upload your data files here:"),
              radioButtons("num_files", label = "Number of Files:",
                           choices = list("One File" = "one",
                                          "Two Files" = "two"), 
                           selected = "one"),
              conditionalPanel(
                condition = "input.num_files == 'one'",
                fileInput("data_file1", "Upload Data File", 
                          accept = c(".csv"))
              ),
              conditionalPanel(
                condition = "input.num_files == 'two'",
                fileInput("data_file1", "Upload Data File 1", 
                          accept = c(".csv")),
                fileInput("data_file2", "Upload Data File 2", 
                          accept = c(".csv"))
              )
      ),
      
      # Normalize tab
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
                  splitLayout(cellWidths = c("100%"),plotOutput("norm_plot") 
                              
                  )
                  
                )
              )
      ),
      
      # EDA tab
      # EDA tab
      tabItem(tabName = "eda",
              h2("differential Expression"),
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
                  splitLayout(cellWidths = c("100%"),tableOutput("de_plot")
                              
                  )
                  
                )
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
      
      # Tools tab
      tabItem(tabName = "tools",
              h2("Tools"),
              p("Choose a network generation tool:"),
              selectInput("tool_choice", label = "Choose Network Generation Tool:", 
                          choices = list("SeqNet", "DIANE"), 
                          selected = "SeqNet")),
      
      # Generate Network tab
      tabItem(tabName = "generate_network",
              h2("Generate Network"),
              
              
              actionButton("network_button", label = "Generate Network", 
                           icon = icon("sitemap"), 
                           style = "background-color: purple; color: white")
              
              ,
              box(
                
                title = "Gene Network",
                status = "primary",
                solidHeader = TRUE,
                width = 20,
                
                fluidRow(
                  splitLayout(cellWidths = c("50%", "50%"),uiOutput("network_plot"),
                              uiOutput("network_plot_2")
                              
                  )
                  
                )
              )
      ),
      
    
     tabItem(tabName = "Log_in",  ##changed name here
            h2("Log in"),
              
        
            div(id = "login_button", 
                
                textInput("username", "Username"),
                passwordInput("password", "Password"),
                br(),
                actionButton("login", "Log In") )
             
     ),

         
         
          
     
     
      
      tabItem(tabName = "Register",
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
      
              
              
            
              
              
             
              
              
              
              
      ))
 


