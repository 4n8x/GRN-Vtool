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
  user_con <- dbConnect(RSQLite::SQLite(), "./data/users.db")  ##"./data/users.db" when we deploy it and rsconnect::deployApp(appDir = ".", appFiles = c("data")) when we deploy the app again in order for the users database to be deployed 
  
  # Create a new table if it doesn't exist
  dbExecute(user_con, "CREATE TABLE IF NOT EXISTS users (id INTEGER PRIMARY KEY, username TEXT UNIQUE, password_hash TEXT)")
  
  # Define an observer to handle login attempts
  observeEvent(input$login, {
    # Check if the credentials match a user in the database
    query <- sprintf("SELECT * FROM users WHERE username = '%s'", input$username)
    result <- dbGetQuery(user_con, query)
    
    if (nrow(result) == 1 && !is.na(result$password_hash) && digest(input$password) == result$password_hash) {
      # If the credentials are valid, store the user ID in the session
      session$user_id <- result$id
      ## session$doRedirect("/main")  ##rawan page but no need for redirection at the moment
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
  })
  
  # Define an observer to handle registration attempts
  observeEvent(input$register, {
    # Hash the password before inserting it into the database
    password_hash <- digest(input$new_password)
    
    # Check if the username already exists in the database
    query <- sprintf("SELECT * FROM users WHERE username = '%s'", input$new_username)
    result <- dbSendStatement(user_con, query)
    user_exists <- nrow(dbFetch(result)) > 0
    dbClearResult(result)
    
    if (user_exists) {
      # If the username already exists, show an error message
      showModal(
        modalDialog(
          title = "Error",
          "Username already exists.",
          easyClose = TRUE
        )
      )
    } else {
      # Insert a new user into the database
      query <- sprintf("INSERT INTO users (username, password_hash) VALUES ('%s', '%s')", input$new_username, password_hash)
      dbExecute(user_con, query)
      
      # Show a success message
      showModal(
        modalDialog(
          title = "Success",
          "Registration successful.",
          easyClose = TRUE
        )
      )
    }
  })
  
  # Close the user database connection when the app is stopped
  onStop(function() {
    dbDisconnect(user_con)
  })  ##end of login and register codes
  
}##end function








