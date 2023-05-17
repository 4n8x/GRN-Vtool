#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#




# Define server logic required to draw a histogram





#---








#---
isSeqNet<-function(V){
  
  if(v=="SeqNet")  
    
    return(TRUE)
}
function(input, output, session) {
  
  
  observeEvent(input$network_button, {
    
    data1<-read_csv(input$data_file1$datapath ,show_col_types = FALSE)
    
    data1<-as.data.frame(data1)
    sapply(data1,class)
    
    
    
    
    
    
    
    p <- 100 # Number of nodes (genes) in the network.
    n <- 100 # Number of samples to generate.
    network <- random_network(p) # Create a random network of size p.
    network <- gen_partial_correlations(network) # Generate weights.
    df_ref <- SeqNet::reference$rnaseq  # The reference dataset
    
    g <-network
    
    
    output$network_plot<-renderPlot({plot(g)})
    
    output$network_plot_2<-renderPlot({plot(g)})
    
    
    
    
  })##end network butten
  
  
  
  ## LOGIN,REGISTER CODES
  
  
  # Create a new user database connection
  user_con <- dbConnect(RSQLite::SQLite(), "./data/users.db")
  
  # Create a new table if it doesn't exist
  dbExecute(user_con, "CREATE TABLE IF NOT EXISTS users (id INTEGER PRIMARY KEY, username TEXT UNIQUE, password_hash TEXT)")
  
  # Define a reactive value to store the user ID and username
  user_info <- reactiveValues(id = NULL, username = NULL)
  
  ## Hide the logout button initially
  ## shinyjs::hideElement("#logout_button")
  hide("logout_button")
  
  # Define an observer to handle login attempts
  observeEvent(input$login, {
    # Check if the input values are not empty
    if (input$username == "" || input$password == "") {
      # If either input value is empty, show an error message
      showModal(
        modalDialog(
          title = "Error",
          "Please enter your username and password.",
          easyClose = TRUE
        )
      )
      
    } else {
      # Check if the credentials match a user in the database
      query <- sprintf("SELECT * FROM users WHERE username = '%s'", input$username)
      result <- dbGetQuery(user_con, query)
      
      if (nrow(result) == 1 && !is.na(result$password_hash) && digest(input$password, algo = "sha256") == result$password_hash) {
        # If the credentials are valid, store the user ID and username in the reactive values
        user_info$id <- result$id
        user_info$username <- input$username
        
        # Show a welcome message in a modal dialog
        showModal(
          modalDialog(
            title = "Welcome",
            paste0("Welcome, ", user_info$username, "!"),
            easyClose = TRUE
          )
        )
        
      } else {
        # If the credentials are invalid, show an error message
        showModal(
          modalDialog(
            title = "Error",
            "Invalid username or password.",
            easyClose = TRUE
          )
        )
      }
      
      # Reset the login form
      updateTextInput(session, "username", value = "")
      updateTextInput(session, "password", value = "")
    }
  })
  
  # Define an observer to handle registration attempts
  observeEvent(input$register, {
    # Check if the input values are not empty
    if (input$new_username == "" || input$new_password == "" || input$confirm_password == "") {
      # If any input values are empty, show an error message
      showModal(
        modalDialog(
          title = "Error",
          "Please enter a username and password, and confirm your password.",
          easyClose = TRUE
        )
      )
      
    } else {
      # Check if the new username is available
      query <- sprintf("SELECT COUNT(*) FROM users WHERE username = '%s'", input$new_username)
      result <- dbGetQuery(user_con, query)
      
      if (result[[1]] != 0) {
        # If the new username is already taken, show an error message
        showModal(
          modalDialog(
            title = "Error",
            "Username is already taken.",
            easyClose = TRUE
          )
        )
        
      } else if (input$new_password != input$confirm_password) {
        # If the passwords do not match, show an error message
        showModal(
          modalDialog(
            title = "Error",
            "Passwords do not match.",
            easyClose = TRUE
          )
        )
        
      } else {
        # If the new username is available and passwords match, add the user to the database
        password_hash <- digest(input$new_password, algo = "sha256")
        query <- sprintf("INSERT INTO users (username, password_hash) VALUES ('%s', '%s')", input$new_username, password_hash)
        dbExecute(user_con, query)
        
        # Log in the user automatically after registration
        user_info$id <- dbGetQuery(user_con, paste0("SELECT id FROM users WHERE username = '", input$new_username, "'"))$id
        user_info$username <- input$new_username
        
        # Show a welcome message in a modal dialog
        showModal(
          modalDialog(
            title = "Welcome",
            paste0("Welcome, ", user_info$username, "!"),
            easyClose = TRUE
          )
        )
        
        # Reset the registration form
        updateTextInput(session, "new_username", value = "")
        updateTextInput(session, "new_password", value = "")
        updateTextInput(session, "confirm_password", value = "")
        hide("register_button") ##hiding here
      }
    }
  })
  
  # Define an observer to handlelogout attempts
  observeEvent(input$logout, {
    # Clear the user info from the reactive values
    user_infoid <- NULL
    user_info$username <- NULL
    
    ##adding message
    showModal(
      modalDialog(
        title = "Thank you and come back!",
        paste0("Thank you ", user_info$username, "!"),
        easyClose = TRUE
      )
    )
    # Show login and registration buttons
    
    ##shinyjs::showElement("#login_buttons")
    show("login_button")
    show("register_button")
    ##shinyjs::hideElement("#logout_button")
    hide("logout_button")
    # Reset the login and registration forms
    updateTextInput(session, "username", value = "")
    updateTextInput(session, "password", value = "")
    updateTextInput(session, "new_username", value = "")
    updateTextInput(session, "new_password", value = "")
    updateTextInput(session, "confirm_password", value = "")
  })
  
  ##  Define a reactive expression to check if the user is logged in
  user_logged_in <- reactive({
    !is.null(user_info$id)
  })
  
  ## Define an observer to show/hide UI elements depending on whether the user is logged in
  observe({
    if (user_logged_in()) {
      # Hide the login and registration buttons, show the logout button
      ##shinyjs::hideElement("#login_buttons")
      hide("login_button")
      hide("register_button")
      ##shinyjs::showElement("#logout_button")
      show("logout_button")
      
      
      
    } else {
      # Show the login and registration buttons, hide the logout button
      ##shinyjs::showElement("#login_buttons")
      show("login_button")
      show("register_button")
      ## shinyjs::hideElement("#logout_button")
      hide("logout_button")
    }
  })
  
  
  
  ## Disconnect from the user database when the app is closed
  session$onSessionEnded(function() {
    dbDisconnect(user_con)
  })  ##end of login and register codes
  
}##end function
