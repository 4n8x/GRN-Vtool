#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


server <- function(input, output, session) {
  
  # Initialize an empty igraph graph object
  g <- reactiveVal(NULL)
  data("abiotic_stresses")
  
  fileUploaded <- reactive({
    uploaded = FALSE
    if (!is.null(input$data_file1)) {
      file1 <- try(read.csv(input$data_file1$datapath,row.names='target_id',fileEncoding="UTF-8-BOM"))
      if (!inherits(file1, "try-error") && nrow(file1) > 0) {
        uploaded = TRUE
      }
    }
    if (!is.null(input$data_file2)) {
      file2 <- try(read.csv(input$data_file2$datapath), silent = TRUE)
      if (!inherits(file2, "try-error") && nrow(file2) > 0) {
        uploaded = TRUE
      }
    }
    uploaded
  })
  
  observeEvent(input$data_file1, {
    if (!is.null(input$data_file1)) {
      file1 <- try(read.csv(input$data_file1$datapath,row.names='target_id',fileEncoding="UTF-8-BOM"))
      if (!inherits(file1, "try-error") && nrow(file1) > 0) {
        abiotic_data <<- file1  # Store data in global variable
        tcc_object <<- DIANE::normalize(abiotic_data, abiotic_stresses$conditions, 
                                        iteration = FALSE)
      } else {
        showModal(modalDialog(
          title = "Warning",
          "The first uploaded file is empty or invalid. Please upload a non-empty CSV file."
        ))
      }
    }
  })
  
  # Observe event for the second file input
  observeEvent(input$data_file2, {
    if (!is.null(input$data_file2)) {
      file2 <- try(read.csv(input$data_file2$datapath), silent = TRUE,row.names='target_id')
      if (!inherits(file2, "try-error") && nrow(file2) > 0) {
        abiotic_data <<- file2  # Store data in global variable
        tcc_object <<- DIANE::normalize(abiotic_data, abiotic_stresses$conditions, 
                                        iteration = FALSE)
      } else {
        showModal(modalDialog(
          title = "Warning",
          "The second uploaded file is empty or invalid. Please upload a non-empty CSV file."
        ))
      }
    }
  })
  
  
  # Create an output to enable conditional panels
  output$fileUploaded <- reactive({
    fileUploaded()
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)
  
  
  output$download_data <- downloadHandler(
    filename = function() {
      paste("abiotic_stresses", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Choose the part of abiotic_stresses you want to download. For example:
      # If you want to download 'raw_counts'
      tcc_object <<- DIANE::filter_low_counts(tcc_object, 10*length(abiotic_stresses$conditions))
      
      normalized_counts <<- TCC::getNormalizedData(tcc_object)
      data_to_download <- normalized_counts
      
      # Or, if you want to download 'normalized_counts'
      # data_to_download <- abiotic_stresses$normalized_counts
      
      # Ensure the data is in a consistent format (data frame)
      if (is.data.frame(data_to_download)) {
        write.csv(data_to_download, file, row.names = FALSE)
      } else {
        stop("Selected data is not a data frame")
      }
    }
  )
  # Generate tool selection input based on user choices
  output$tool_input <- renderUI({
    if (input$num_files == "one" | input$num_files == "zero") {
      checkboxGroupInput("tool_choice", label = "Choose Network Generation Tool:",
                         choices = list("DIANE" = "diane", "SeqNet" = "seqnet")
      )
    } else {
      radioButtons("tool_choice", label = "Choose Network Generation Tool:",
                   choices = list("DIANE" = "diane", "SeqNet" = "seqnet")
      )
    }
  })
  
  
  # Generate network boxes based on tool choices
  output$diane_network_box <- renderUI({
    if ("diane" %in% input$tool_choice) {
      if (input$num_files == "one") {
        box(
          title = "Gene Network (DIANE)",
          status = "primary",
          solidHeader = TRUE,
          width = 20, # Increase the width to 20 columns for bigger box
          plotOutput("diane_network_plot"),
          align = "center" # Center align the box
        )
      } else if (input$num_files == "two") {
        fluidRow(
          column(6, box(
            title = "Gene Network (DIANE - File 1)",
            status = "primary",
            solidHeader = TRUE,
            width = 20, # Increase the width to 20 columns for bigger box
            plotOutput("diane_network_plot_1"),
            align = "center" # Center align the box
          )),
          column(6, box(
            title = "Gene Network (DIANE - File 2)",
            status = "primary",
            solidHeader = TRUE,
            width = 20, # Increase the width to 20 columns for bigger box
            plotOutput("diane_network_plot_2"),
            align = "center" # Center align the box
          ))
        )
      }
    }
  })
  
  output$seqnet_network_box <- renderUI({
    
    if ("seqnet" %in% input$tool_choice) {
      if (input$num_files == "one") {
        box(
          title = "Gene Network (SeqNet)",
          status = "primary",
          solidHeader = TRUE,
          width = 20, # Increase the width to 20 columns for bigger box
          visNetworkOutput("seqnet_network_plot"),
          align = "center" # Center align the box
        )
      } else if (input$num_files == "two") {
        fluidRow(
          column(6, box(
            title = "Gene Network (SeqNet - File 1)",
            status = "primary",
            solidHeader = TRUE,
            width = 20, # Increase the width to 20 columns for bigger box
            visNetworkOutput("seqnet_network_plot_1"),
            align = "center" # Center align the box
          )),
          column(6, box(
            title = "Gene Network (SeqNet - File 2)",
            status = "primary",
            solidHeader = TRUE,
            width = 20, # Increase the width to 20 columns for bigger box
            visNetworkOutput("seqnet_network_plot_2"),
            align = "center" # Center align the box
          ))
        )
      }
    }
  })
  
  nodes_df_reactive <- reactiveVal()
  edges_df_reactive <- reactiveVal()
  
  
  
  #DIANE code
  
  
  observeEvent(input$norm_button, {
    
    # Load the data directly from the SeqNet package
    if (fileUploaded()) {
      if ("seqnet" %in% input$tool_choice){
        print("norm")
        
        
        # Extract the normalized counts from the tcc_object
        normalized_counts <<- TCC::getNormalizedData(tcc_object)
        
        # Output the normalized data to a table in the UI
        output$norm_plot<<-renderPlot({plot(tcc_object)})}
      
      
    }
    else{
      showModal(modalDialog(
        title = "Warning",
        "Please upload a file before normalizing data."
      ))
    }
  })
  
  observeEvent(input$de_button, {
    if (fileUploaded()) {
      if (input$num_files == "one"| input$num_files == "two") {
        
        if ("seqnet" %in% input$tool_choice){
          print("differential started")
          # Load the data directly from the SeqNet package
          #@ data("abiotic_stresses")
          
          tcc_object <<- DIANE::filter_low_counts(tcc_object, 10*length(abiotic_stresses$conditions))
          fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = abiotic_stresses$conditions)
          print("fit done")
          
          tags <- DIANE::estimateDEGs(fit, reference = "C", perturbation = "H", p.value = 1)
          print("tags done")
          
          output$de_plot<<-renderPlot({ DIANE::draw_DEGs(tags, fdr = 0.01, lfc = 1)})}
        else if ("diane" %in% input$tool_choice) {
          
          
          
        }
        else{
          showModal(modalDialog(
            title = "Warning",
            "Please chose a method before generating Differential Expression."
          ))
        }
      }
      
    }
    else{
      showModal(modalDialog(
        title = "Warning",
        "Please upload a file before generating Differential Expression."
      ))
    }
  })
  observeEvent(input$ebc_button, {
    if (fileUploaded()) {
      print("clustering started")
      if (input$num_files == "one" | input$num_files == "two") {
        
        if ("seqnet" %in% input$tool_choice){
          # data("abiotic_stresses")
          genes <<- abiotic_stresses$heat_DEGs
          clustering <<- run_coseq(conds = unique(abiotic_stresses$conditions), 
                                   data = abiotic_stresses$normalized_counts, genes = genes, K = 6:9)
          
          #> ****************************************
          #> coseq analysis: Normal approach & arcsin transformation
          #> K = 6 to 9 
          #> Use seed argument in coseq for reproducible results.
          #> ****************************************
          
          print("clustering done")
          output$ebc_plot<<-renderPlot({DIANE::draw_profiles(data =abiotic_stresses$normalized_counts, clustering$membership, conds = unique(abiotic_stresses$conditions), k = 3) })
        }
        else  if ("diane" %in% input$tool_choice) {
          
          
          
          
          
          
          
        }
        else{
          showModal(modalDialog(
            title = "Warning",
            "Please chose a method before generating Differential Expression."
          ))
        }
      }
      
    }
    else{
      showModal(modalDialog(
        title = "Warning",
        "Please upload a file before clustering data."
      ))
    }
    
  })
  #> idhar hai network
  
  observeEvent(input$generate_button, {
    if (is.null(input$tool_choice)) {
      showModal(modalDialog(
        title = "Tool Selection Required",
        "Please select a network generation tool from the Tools menu."
      ))
    }
  })
  
  # For some specific action in your app
  observeEvent(input$some_action_button, {
    if (!is.null(input$some_required_input)) {
      # Proceed with the action
    } else {
      showModal(modalDialog(
        title = "Missing Information",
        "Please provide all the required information before proceeding."
      ))
    }
  })
  
  
  observeEvent(input$generate_button, {
    if ("seqnet" %in% input$tool_choice) {
      if (fileUploaded()) {
        
        
        data("abiotic_stresses")
        edges <- abiotic_stresses$heat_DEGs_regulatory_links 
        # Get the row names of the matrix 'edges'
        edges <- edges[-(19:25), ]
        
        row_names <- rownames(edges)
        
        # Print the row names
        print(row_names)
        column_names <- colnames(edges)
        
        # Print the column names
        print(column_names)
        common_genes <- intersect(row_names, column_names)
        
        edges <- edges[common_genes, common_genes]
        row_names2 <- rownames(edges)
        
        # Print the row names
        print(row_names2)
        column_names2 <- colnames(edges)
        
        # Print the column names
        print(column_names2)
        symmetric_edges <- pmax(edges, t(edges))
        
        threshold <- 0  # setting threshold
        binary_edges <- ifelse(symmetric_edges > threshold, 1, 0)
        
        # create the network with the adjacency matrix
        nw <- create_network_from_adjacency_matrix(binary_edges)
        
        
        nodes_df <- data.frame(
          id = nw$node_names,  # row names from matrix
          label = nw$node_names,
          name = nw$node_names,
          color = "orange",  # Assigns a default color to each node
          size = 30,  # Assigns a default size to each node
          value = runif(nrow(binary_edges))  # Assigns a random value to each node for now
        )
        
        nodes_df_reactive(nodes_df)
        
        # Create edges_df
        # create a data frame from the matrix where each '1' represents an edge
        edges_data <- nw$modules[[1]]$edges
        
        # If the edges are in a matrix form where each row represents an edge and the first two columns
        # are the indices of the connected nodes,converting it to a data frame:
        edges_df <- data.frame(
          from = edges_data[, 1],
          to = edges_data[, 2]
        )
        
        # if there are corresponding node names, it can map them
        edges_df$from <- nw$node_names[edges_df$from]
        edges_df$to <- nw$node_names[edges_df$to]
        
        # Add a color attribute to edges
        edges_df$color <- "grey"
        get_connected_nodes <- function(clicked_node, edges_df) {
          # Find rows in edges_df where the clicked node is either 'from' or 'to'
          connected <- subset(edges_df, from == clicked_node | to == clicked_node)
          # Return the unique node names connected to the clicked node
          unique(c(connected$from, connected$to))
        }
        
        edges_df_reactive(edges_df)
        observe({
          # Assuming nodes_df has an 'id' column with the node ids
          updateSelectInput(session, "nodeSelector",
                            choices = nodes_df$id)
        })
        observeEvent(input$nodeSelector, {
          print("here")
          print(input$nodeSelector)
          connected_nodes_node <- get_connected_nodes(input$nodeSelector, edges_df)
          
          showModal(modalDialog(
            title = "Node Information",
            paste("Clicked Node:", input$nodeSelector, 
                  "\nConnected Nodes:", paste(connected_nodes_node, collapse = ", "))
          ))
        })
        observeEvent(input$interactive_graph_nodes, {
          if(!identical(input$interactive_graph_nodes$nodes, list())) {
            clicked_node <- input$interactive_graph_nodes$nodes[1]
            
            connected_nodes_node <- get_connected_nodes(clicked_node, edges_df)
            
            showModal(modalDialog(
              title = "Node Information",
              paste("Clicked Node:", clicked_node, 
                    "\nConnected Nodes:", paste(connected_nodes_node, collapse = ", "))
            ))
          }
        }, ignoreNULL = TRUE)
        output$seqnet_network_plot <- renderVisNetwork({
          visNetwork(nodes_df, edges_df, width = "100%", height = "8000px") %>%
            visIgraphLayout() %>%
            visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
            visInteraction(dragNodes = TRUE) %>%
            visEvents(click = "function(nodes){
    Shiny.setInputValue('interactive_graph_nodes', nodes, {priority: 'event'});
  }")
        })
        shinyjs::hide("generate_button")
      }
      else{
        showModal(modalDialog(
          title = "Warning",
          "Please upload a file before generating network."
        ))
      }
    }
    else if ("seqnet" %in% input$tool_choice) {
      print("SeqNet network bein generated here")
      if (fileUploaded()) {
        # data("abiotic_stresses")
        edges <- abiotic_stresses$heat_DEGs_regulatory_links 
        # Get the row names of the matrix 'edges'
        edges <- edges[-(19:25), ]
        
        row_names <- rownames(edges)
        
        # Print the row names
        print(row_names)
        column_names <- colnames(edges)
        
        # Print the column names
        print(column_names)
        common_genes <- intersect(row_names, column_names)
        
        edges <- edges[common_genes, common_genes]
        row_names2 <- rownames(edges)
        
        # Print the row names
        print(row_names2)
        column_names2 <- colnames(edges)
        
        # Print the column names
        print(column_names2)
        symmetric_edges <- pmax(edges, t(edges))
        
        threshold <- 0  # setting threshold
        binary_edges <- ifelse(symmetric_edges > threshold, 1, 0)
        
        # create the network with the adjacency matrix
        nw <- create_network_from_adjacency_matrix(binary_edges)
        
        
        nodes_df <- data.frame(
          id = nw$node_names,  # row names from matrix
          label = nw$node_names,
          name = nw$node_names,
          color = "orange",  # Assigns a default color to each node
          size = 30,  # Assigns a default size to each node
          value = runif(nrow(binary_edges))  # Assigns a random value to each node for now
        )
        
        nodes_df_reactive(nodes_df)
        
        # Create edges_df
        # create a data frame from the matrix where each '1' represents an edge
        edges_data <- nw$modules[[1]]$edges
        
        # If the edges are in a matrix form where each row represents an edge and the first two columns
        # are the indices of the connected nodes,converting it to a data frame:
        edges_df <- data.frame(
          from = edges_data[, 1],
          to = edges_data[, 2]
        )
        
        # if there are corresponding node names, it can map them
        edges_df$from <- nw$node_names[edges_df$from]
        edges_df$to <- nw$node_names[edges_df$to]
        
        # Add a color attribute to edges
        edges_df$color <- "grey"
        get_connected_nodes <- function(clicked_node, edges_df) {
          # Find rows in edges_df where the clicked node is either 'from' or 'to'
          connected <- subset(edges_df, from == clicked_node | to == clicked_node)
          # Return the unique node names connected to the clicked node
          unique(c(connected$from, connected$to))
        }
        
        edges_df_reactive(edges_df)
        
        observeEvent(input$interactive_graph_nodes, {
          if(!identical(input$interactive_graph_nodes$nodes, list())) {
            clicked_node <- input$interactive_graph_nodes$nodes[1]
            
            connected_nodes_node <- get_connected_nodes(clicked_node, edges_df)
            
            showModal(modalDialog(
              title = "Node Information",
              paste("Clicked Node:", clicked_node, 
                    "\nConnected Nodes:", paste(connected_nodes_node, collapse = ", "))
            ))
          }
        }, ignoreNULL = TRUE)
        output$seqnet_network_plot <- renderVisNetwork({
          visNetwork(nodes_df, edges_df, width = "100%", height = "8000px") %>%
            visIgraphLayout() %>%
            visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
            visInteraction(dragNodes = TRUE) %>%
            visEvents(click = "function(nodes){
    Shiny.setInputValue('interactive_graph_nodes', nodes, {priority: 'event'});
  }")
        })
        shinyjs::hide("generate_button")
      }
      else{
        showModal(modalDialog(
          title = "Warning",
          "Please upload a file before generating network."
        ))
      }
    }
  })
  
  
  #end of diane 
  
  
  
  observeEvent(input$change_color, {
    if ("seqnet" %in% input$tool_choice) {
      if (!is.null(nodes_df_reactive())) {
        print("hiding the button")
        shinyjs::hide("change_color")
        shinyjs::show("revert_color")
        updated_nodes_df <- nodes_df_reactive()
        updated_nodes_df$color <- "red"
        nodes_df_reactive(updated_nodes_df)
        get_connected_nodes <- function(clicked_node, edges_df_closure) {
          # Call the closure to get the dataframe
          edges_df <- edges_df_closure()
          
          # Ensure edges_df is a dataframe
          if (!is.data.frame(edges_df)) {
            stop("The object returned by edges_df_closure is not a dataframe")
          }
          
          # Find rows in edges_df where the clicked node is either 'from' or 'to'
          connected <- subset(edges_df, edges_df$from == clicked_node | edges_df$to == clicked_node)
          
          # Print edges_df for debugging
          # print(edges_df)
          
          # Define your data array
          data <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
          
          # Return the data array
          unique(c(connected$from, connected$to))
          
          
        }
        
        
        
        
        observeEvent(input$interactive_graph_nodes, {
          if(!identical(input$interactive_graph_nodes$nodes, list())) {
            clicked_node <- input$interactive_graph_nodes$nodes[1]
            
            connected_nodes_node <- get_connected_nodes(clicked_node, edges_df_reactive)
            
            showModal(modalDialog(
              title = "Node Information",
              paste("Clicked Node:", clicked_node, 
                    "\nConnected Nodes:", paste(connected_nodes_node, collapse = ", "))
            ))
          }
        }, ignoreNULL = TRUE)
        
        output$seqnet_network_plot <- renderVisNetwork({
          visNetwork(nodes_df_reactive(), edges_df_reactive(), width = "100%", height = "8000px") %>%
            visIgraphLayout() %>%
            visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
            visInteraction(dragNodes = TRUE) %>%
            visEvents(click = "function(nodes){
    Shiny.setInputValue('interactive_graph_nodes', nodes, {priority: 'event'});
  }")
          
        })
      }
      
    }
  })
  
  observeEvent(input$revert_color, {
    if ("seqnet" %in% input$tool_choice) {
      if (!is.null(nodes_df_reactive())) {
        updated_nodes_df <- nodes_df_reactive()
        updated_nodes_df$color <- "orange"  # Revert color to orange
        nodes_df_reactive(updated_nodes_df)
        shinyjs::show("change_color")
        shinyjs::hide("revert_color")
        
        # Update the network plot
        output$seqnet_network_plot <- renderVisNetwork({
          visNetwork(nodes_df_reactive(), edges_df_reactive(), width = "100%", height = "8000px") %>%
            visIgraphLayout() %>%
            visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
            visInteraction(dragNodes = TRUE) %>%
            visEvents(click = "function(nodes){
            Shiny.setInputValue('interactive_graph_nodes', nodes, {priority: 'event'});
          }")
        })
      }
    }
  })
  
  
  
  observeEvent(input$get_sigma_button, {
    if ("seqnet" %in% input$tool_choice) {
      g_copy <- g()
      if (!is.null(g_copy)) {
        # Calculate partial correlations for the network
        nw <- gen_partial_correlations(g_copy)
        sigma_network <- get_sigma(nw)
        
        # Display the covariance matrix in a modal dialog
        showModal(
          modalDialog(
            title = "Covariance Matrix for Network",
            tableOutput("sigma_network_table")
          )
        )
        
        output$sigma_network_table <- renderTable({
          sigma_network
        })
      }
    }
  })
  ##Database connection,login,register
  
  user_con <- dbConnect(RMySQL::MySQL(),
                        dbname = "login",
                        host = "login-grnvtools.cs2w1zfwm4d4.eu-north-1.rds.amazonaws.com",
                        port = 3306,
                        user = "admin",
                        password = "rS268043!")
  
  
  # Create a new table if it doesn't exist
  dbExecute(user_con, "CREATE TABLE IF NOT EXISTS users (id INTEGER PRIMARY KEY, username TINYTEXT , password_hash TEXT)")
  
  user_info <- reactiveValues(id = NULL, username = NULL)
  
  user_logged_in <- reactive({
    !is.null(user_info$id)
  })
  
  
  
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
      input_password <- charToRaw(input$password)
      input_username <- charToRaw(input$username)
      database_username <- charToRaw(result$username)
      if (nrow(result) == 1 && !is.na(result$password_hash) && digest(input$password, algo = "sha256") == result$password_hash && (database_username == input_username)) {
        # If the credentials are valid, store the user ID and username in the reactive values
        user_info$id <- result$id
        user_info$username <- input$username
        ##edited
        
        
        # Show a welcome message in a modal dialog
        showModal(
          modalDialog(
            title = "Welcome",
            paste0("Welcome, ", user_info$username, "!, proceed to visualize GRN starting from Data tab "),
            easyClose = TRUE
          )
        )
        shinyjs::show("logout")
        shinyjs::show("data_menu")
        shinyjs::show("tools_menu")
        shinyjs::show("network_menu")
        shinyjs::show("download_menu")
        shinyjs::hide("login_button")
        shinyjs::hide("register_button")
        shinyjs::hide("register_menu")
        shinyjs::hide("login_menu")
        
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
  })##observe login attempts
  
  
  # Define an observer to handle registration attempts
  observeEvent(input$register, {
    length_username <- nchar(input$new_username)
    length_password <- nchar(input$new_password)
    special_char <- input$new_password
    capital_letter <- input$new_password
    special_char_count <- nchar(gsub("[[:alnum:]\\s]", "", special_char))
    
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
      
    }
    
    ##Adding password and username being 8 characters long
    
    else if (length_username <8 || length_password < 8) {
      # If any input values are less than 8 characters show error message
      showModal(
        modalDialog(
          title = "Error",
          "Please enter a username and password that are at least 8 characters long",
          easyClose = TRUE
        )
      )
      
    }
    
    ##Adding passwords having special characters (at least 2)
    
    else if (special_char_count < 2) {
      
      showModal(
        modalDialog(
          title = "Error",
          "Please enter a password that has at least 2 special character.",
          easyClose = TRUE
        )
      )
      
    } 
    
    else if(!grepl("[A-Z]", capital_letter)){
      
      showModal(
        modalDialog(
          title = "Error",
          "Please enter a password that has at least 1 capital letter.",
          easyClose = TRUE
        )
      )
      
    }
    
    
    
    else {
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
        
        
        
        
        showModal(
          modalDialog(
            title = "Welcome",
            paste0("Welcome, ", user_info$username, "!, proceed to visualize GRN starting from Data tab."),
            easyClose = TRUE
          )
        )
        
        
        ##Hide/show test later
        
        
        # Reset the registration form
        updateTextInput(session, "new_username", value = "")
        updateTextInput(session, "new_password", value = "")
        updateTextInput(session, "confirm_password", value = "")
        
        
      }
    }
  }) #Observe register attempts
  
  
  # Clear the user info from the reactive values
  observeEvent(input$logout, {
    
    
    user_info$id <- NULL
    user_info$username <- NULL
    
    ##adding message
    showModal(
      modalDialog(
        title = "Thank you for using GRNVtool!",
        paste0("Thank you", user_info$username, "!"),
        easyClose = TRUE
      )
    )
    
    
    shinyjs::runjs("setTimeout(function() { window.location.href = '/' }, 2000);")
    
    
    
    # Reset the login and registration form (comment for now)
    updateTextInput(session, "username", value = "")
    updateTextInput(session, "password", value = "")
    updateTextInput(session, "new_username", value = "")
    updateTextInput(session, "new_password", value = "")
    updateTextInput(session, "confirm_password", value = "")
    
  })
  
  ##redirect test 
  
  
  
  observe({
    if (user_logged_in()) {
      
      shinyjs::hide("login_button")
      shinyjs::hide("register_button")
      shinyjs::hide("login_menu")
      shinyjs::hide("register_menu")
      shinyjs::show("logout")
      shinyjs::show("data_item")
      shinyjs::show("data_menu")
      shinyjs::show("tools_item")
      shinyjs::show("network_item")
      shinyjs::show("download_item")
      shinyjs::show("tools_menu")
      shinyjs::show("logout_button")
      
      shinyjs::show("network_menu")
      shinyjs::show("download_menu")
    } else{
      
      
      shinyjs::show("login_button")
      shinyjs::show("register_button")
      
      shinyjs::hide("logout")
      shinyjs::hide("data_menu")
      shinyjs::hide("tools_menu")
      shinyjs::hide("network_menu")
      shinyjs::hide("download_menu")
      
      
      
    }
    
  })
  
  ## Disconnect from the user database when the app is closed
  session$onSessionEnded(function() {
    dbDisconnect(user_con)
  })  ##end of login and register codes
  
  
  
  
}
