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
  
  
  # Generate tool selection input based on user choices
  output$tool_input <- renderUI({
    if (input$num_files == "one") {
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
  
  
  
  
  observeEvent(input$generate_button, {
    if ("seqnet" %in% input$tool_choice) {
      # Generate a random network
      p <- 100
      n <- 100
      network <- random_network(p)
      network <- gen_partial_correlations(network)
      g(network)  # Assign the generated network to the 'g' object
      
      # Create nodes_df with a full_label column for concatenated labels
      nodes_df <- data.frame(id = 1:length(network$node_names), 
                             label = paste( 1:length(network$node_names)),
                             name = paste0(network$node_names), 
                             color = "orange", 
                             size = 30, 
                             value = runif(length(network$node_names)))
      
      # Create edges_df
      edges_list <- lapply(network$modules, function(module) {
        edges_df <- as.data.frame(module$edges[, 1:2])
        colnames(edges_df) <- c("from", "to")
        edges_df$color <- "grey"
        return(edges_df)
      })
      edges_df <- do.call(rbind, edges_list)
      
      # Filter nodes with no edges
      nodes_with_edges <- unique(c(edges_df$from, edges_df$to))
      nodes_df <- nodes_df[nodes_df$id %in% nodes_with_edges, ]
      
      
      
      observeEvent(input$interactive_graph_nodes, {
        print(input$interactive_graph_nodes$nodes)
        if(!identical(input$interactive_graph_nodes$nodes, list())) {
          clicked_node <- input$interactive_graph_nodes$nodes[1]
          list_length <- length(input$interactive_graph_nodes$edges)
          showModal(modalDialog(
            title = "Node Information",
            paste("Clicked Node:", clicked_node, 
                  "\nConnected Nodes:", list_length)
          ))
        }
      }, ignoreNULL = TRUE)
      
      # Render visNetwork
      print(nodes_df$label)
      
      
      if (input$num_files == "one") {
        output$seqnet_network_plot <- renderVisNetwork({
          if (!is.null(g())) {
            visNetwork(nodes_df, edges_df, width = "100%", height = "8000px") %>%
              visIgraphLayout() %>%
              visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
              visInteraction(dragNodes = TRUE) %>%
              visEvents(click = "function(nodes){
        Shiny.setInputValue('interactive_graph_nodes', nodes, {priority: 'event'});
      }")
          }
        })
      } else if (input$num_files == "two") {
        output$seqnet_network_plot_1 <- renderVisNetwork({
          if (!is.null(g())) {
            visNetwork(nodes_df, edges_df, width = "100%", height = "8000px") %>%
              visIgraphLayout() %>%
              visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
              visInteraction(dragNodes = TRUE) %>%
              visEvents(click = "function(nodes){
        Shiny.setInputValue('interactive_graph_nodes', nodes, {priority: 'event'});
      }")
          }
        })
        output$seqnet_network_plot_2 <- renderVisNetwork({
          if (!is.null(g())) {
            visNetwork(nodes_df, edges_df, width = "100%", height = "8000px") %>%
              visIgraphLayout() %>%
              visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
              visInteraction(dragNodes = TRUE) %>%
              visEvents(click = "function(nodes){
        Shiny.setInputValue('interactive_graph_nodes', nodes, {priority: 'event'});
      }")
          }
        })
      }
    }
  })
  
  #DIANE code
  
  
  observeEvent(input$norm_button, {
    
    
    
    data1<<-read_csv(input$data_file1$datapath ,show_col_types = FALSE)
    
    data1<<-as.data.frame(data1)
    sapply(data1,class)
    
    
    data("abiotic_stresses")
    data("gene_annotations")
    data("regulators_per_organism")
    
    
    DIANE::draw_distributions(abiotic_stresses$raw_counts, boxplot = TRUE)
    #normalization using the defult paramaters
    tcc_object <<- DIANE::normalize(abiotic_stresses$raw_counts, norm_method = 'tmm', iteration = FALSE)
    
    #Low counts removal (allow only genes with more than 10 counts per sample in average)
    
    threshold = 10*length(abiotic_stresses$conditions)
    tcc_object <<- DIANE::filter_low_counts(tcc_object, threshold)
    normalized_counts <<- TCC::getNormalizedData(tcc_object)
    
    
    output$norm_plot<<-renderPlot({plot(tcc_object)})
    #--------------------------------------------------
    data <- reactive({
      req(input$data_file1)
      raw_data <- read.csv(input$data_file1$datapath, row.names = 1, header = TRUE, sep = ",")
      return(raw_data)
    })
    
    normalized_data <- eventReactive(input$norm_button, {
      # Load data
      counts <- data()
      
      # Perform normalization using edgeR's TMM (Trimmed Mean of M-values) method
      normalized_counts <- calcNormFactors(counts)
      
      return(normalized_counts)
    })
    output$normalized_table <- renderTable({
      # Display the normalized data in a table format
      # Use 'normalized_data()' to access the normalized data
      normalized_data()
    })
    
    
    
    
  })
  observeEvent(input$de_button, {
    if ("diane" %in% input$tool_choice) {
      
      data("abiotic_stresses")
      data("gene_annotations")
      data("regulators_per_organism")
      
      #Differential expression analysis
      fit <<- DIANE::estimateDispersion(tcc = tcc_object, conditions = abiotic_stresses$conditions)
      #> Warning in edgeR::DGEList(counts = tcc$count, lib.size = tcc$norm.factors, :
      #> norm factors don't multiply to 1
      
      
      
      topTags <<- DIANE::estimateDEGs(fit, reference = "C", perturbation = "H", p.value = 0.01, lfc = 2)
      
      # adding annotations
      DEgenes <<- topTags$table
      DEgenes[,c("name", "description")] <- gene_annotations$`Arabidopsis thaliana`[
        match(get_locus(DEgenes$genes, unique = FALSE), rownames(gene_annotations$`Arabidopsis thaliana`)),
        c("label", "description")]
      
      knitr::kable(head(DEgenes, n = 10))
      
      #Network inference
      #Regulators for your organism
      aggregated_data <<- aggregate_splice_variants(data = normalized_counts)
      
      #Grouping highly correlated regulators
      
      genes <<- get_locus(topTags$table$genes)
      regressors <<- intersect(genes, regulators_per_organism[["Arabidopsis thaliana"]])
      
      # use normalized counts if you did not aggregate splice variants
      grouping <- DIANE::group_regressors(aggregated_data, genes, regressors)
      #> [1] "adding tf  AT1G74890  to group  2 because correlation of -0.912173913043478 to mean"
      #> [1] 19
      #> [1] "adding tf  AT5G59570  to group  4 because correlation of -0.919130434782609 to mean"
      #> [1] 18
      
      grouped_counts <<- grouping$counts
      
      grouped_targets <<- grouping$grouped_genes
      
      grouped_regressors <<- grouping$grouped_regressors
      
      
      genes <<- topTags$table$genes
      
      tags <- DIANE::estimateDEGs(fit, reference = "C", perturbation = "H", p.value = 1)
      
      
      output$de_plot<<-renderPlot({ DIANE::draw_DEGs(tags)})
      
      
      
      
    }
  })
  
  observeEvent(input$ebc_button, {
    if ("diane" %in% input$tool_choice) {
      
      
      #Expression based clustering
      genes <- topTags$table$genes
      
      clustering <- DIANE::run_coseq(conds = unique(abiotic_stresses$conditions), data = normalized_counts, 
                                     genes = genes, K = 6:9, transfo = "arcsin", model = "Normal", seed = 123)
      #> ****************************************
      #> coseq analysis: Normal approach & arcsin transformation
      #> K = 6 to 9 
      #> Use seed argument in coseq for reproducible results.
      #> ****************************************
      
      
      
      
      
      
      
      
      output$ebc_plot<<-renderPlot({DIANE::draw_profiles(data = normalized_counts, clustering$membership, conds = unique(abiotic_stresses$conditions), k = 3) })
      
      
      
      
      
      
      
      
      
    }
  })
  
  observeEvent(input$generate_button, {
    if ("diane" %in% input$tool_choice) {
      
      
      set.seed(123)
      mat <<- DIANE::network_inference(grouped_counts, conds = abiotic_stresses$conditions, 
                                       targets = grouped_targets, regressors = grouped_regressors, 
                                       nCores = 1, verbose = FALSE)
      
      network <<- DIANE::network_thresholding(mat, n_edges = length(genes))
      
      data <- network_data(network, regulators_per_organism[["Arabidopsis thaliana"]], 
                           gene_annotations$`Arabidopsis thaliana`)
      
      knitr::kable(head(data$nodes))
      genes <<- topTags$table$genes
      
      output$eda_plot<<-renderUI( {DIANE::draw_heatmap(normalized_counts, subset = genes, 
                                                       title = "Log expression for DE genes under heat stress")})
      
      #viualize
      #(the interactive graph network will be seen here)
      DIANE::draw_network(data$nodes, data$edges)
      
      nGenes = length(grouped_targets)
      nRegulators = length(grouped_regressors)
      
      res <<- data.frame(density = seq(0.001, 0.1, length.out = 20),
                         nEdges = sapply(seq(0.001, 0.1, length.out = 20),
                                         get_nEdges, nGenes, nRegulators))
      
      
      
      
      output$diane_network_box<<-renderUI({DIANE::draw_network(data$nodes, data$edges)})
      
      
      
      
      
    }
  })
  
  
  #end of diane 
  
  
  
  observeEvent(input$change_color, {
    if ("seqnet" %in% input$tool_choice) {
      g_copy <- g()
      if (!is.null(g_copy)) {
        # Modify the plot using plot_modules
        modified_network <- plot_modules(
          g_copy,
          node_scale = 4,
          edge_scale = 1,
          node_color = adjustcolor("orange", 0.5),
          group_color = palette.colors(9, "Set 1"),
          generate_layout = igraph::nicely,
          include_vertex_labels = TRUE,
          show_legend = FALSE,
          legend_position = "topright",
          legend_horizontal = FALSE,
          display_plot = FALSE  # Don't display the plot immediately
        )
        
        if (input$num_files == "one") {
          output$seqnet_network_plot <- renderPlot({
            print(modified_network)  # Display the modified network plot
          })
        } else if (input$num_files == "two") {
          output$seqnet_network_plot_1 <- renderPlot({
            print(modified_network)  # Display the modified network plot
          })
          output$seqnet_network_plot_2 <- renderPlot({
            print(modified_network)  # Display the modified network plot
          })
        }
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
