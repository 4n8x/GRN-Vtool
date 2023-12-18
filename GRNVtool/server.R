#
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
  first_click <- reactiveVal(FALSE)
  can_download1 <<- reactiveVal(FALSE)
  can_download2 <<- reactiveVal(FALSE)
  can_download3 <<- reactiveVal(FALSE)
  can_download4 <<- reactiveVal(FALSE)
  
  user_selected_node <<- reactiveVal(FALSE)
  
  observe({
    if (input$num_files == "two") {
      shinyjs::show("data_file2")
    } else {
      shinyjs::hide("data_file2")
    }
  })
  fileUploaded <- reactive({
    uploaded = FALSE
    if (!is.null(input$data_file1)) {
      print(input$data_file1$datapath)
      file1 <<- try(read.csv(input$data_file1$datapath,row.names='target_id',fileEncoding="UTF-8-BOM"))
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
        print(class(tcc_object))
        expression_data <- file1
        
        
        regulators <- targets <- row_ids <- rownames(file1)
        print(abiotic_stresses$conditions)
        network_matrix <- network_inference(expression_data, 
                                            conds =abiotic_stresses$conditions,
                                            targets = targets, 
                                            regressors = regulators,
                                            nTrees = 10, 
                                            nCores = 4)
        
        # Step 4: Threshold the Network
        # Thresholding to keep significant edges
        network <- network_thresholding(network_matrix, n_edges = 100)  # Adjust n_edges as needed
        
        # Step 5: Prepare Data for Visualization
        network_data_for_vis <- network_data(network, regulators)
        
        # Step 6: Visualize the Network
        nodes <- network_data_for_vis$nodes
        edges2 <<- network_data_for_vis$edges
      } else {
        showModal(modalDialog(
          title = "Warning",
          "The first uploaded file is invalid, Upload a valid CSV file."
        ))
      }
    }
  })
  
  observeEvent(input$data_file2, {
    if (!is.null(input$data_file2)) {
      file2 <<- try(read.csv(input$data_file2$datapath,row.names='target_id',fileEncoding="UTF-8-BOM"))
      if (!inherits(file2, "try-error") && nrow(file2) > 0) {
        abiotic_data2 <<- file2  # Store data in global variable
        tcc_object2 <<- DIANE::normalize(abiotic_data2, abiotic_stresses$conditions, 
                                         iteration = FALSE)
      } else {
        showModal(modalDialog(
          title = "Warning",
          "The second uploaded file is invalid, Upload a valid CSV file."
        ))
      }
    }
  })
  
  
  # Create an output to enable conditional panels
  output$fileUploaded <- reactive({
    fileUploaded()
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)
  
  
  observeEvent(input$download_data4, {
    print("here")
    # all_conditions_met <- can_download1() && can_download2() && can_download3() && can_download4()
    #   print(all_conditions_met)
    data_to_download <- edges2
    file_path <- file.path(getwd(), "nodes.csv")
    write.csv(data_to_download, file_path, row.names = FALSE)
    
  })
  
  observeEvent(input$download_data, {
    print("here")
    #   print(all_conditions_met)
    data_to_download <- normalized_counts
    file_path <- file.path(getwd(), "norm.csv")
    write.csv(data_to_download, file_path, row.names = FALSE)
    
  })
  observeEvent(input$download_data2, {
    print("here")
    # all_conditions_met <- can_download1() && can_download2() && can_download3() && can_download4()
    #   print(all_conditions_met)
    fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = abiotic_stresses$conditions)
    
    
    tags <- DIANE::estimateDEGs(fit, reference = "C", perturbation = "H", p.value = 1)
    
    file_path <- file.path(getwd(), "diffrential.csv")
    write.csv(tags$table, file_path, row.names = FALSE)
    
    
  })
  observeEvent(input$download_data3, {
    print("here")
    # all_conditions_met <- can_download1() && can_download2() && can_download3() && can_download4()
    #   print(all_conditions_met)
    genes <<- abiotic_stresses$heat_DEGs
    
    clustering <<- run_coseq(conds = unique(abiotic_stresses$conditions), 
                             data = abiotic_stresses$normalized_counts, genes = genes, K = 6:9)
    
    print(clustering$model@allResults)
    
    data_to_download <- clustering$model@allResults
    file_path <- file.path(getwd(), "clustering.csv")
    write.csv(data_to_download, file_path, row.names = FALSE)
    
  })
  observeEvent(input$download_data2, {
    print("here")
    # all_conditions_met <- can_download1() && can_download2() && can_download3() && can_download4()
    #   print(all_conditions_met)
    data_to_download <- file1
    file_path <- file.path(getwd(), "data.csv")
    write.csv(data_to_download, file_path, row.names = FALSE)
    
  })
  observeEvent(input$download_data, {
    print("here")
    # all_conditions_met <- can_download1() && can_download2() && can_download3() && can_download4()
    #   print(all_conditions_met)
    data_to_download <- file1
    file_path <- file.path(getwd(), "data.csv")
    write.csv(data_to_download, file_path, row.names = FALSE)
    
  })
  
  
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
          visNetworkOutput("diane_network_plot"),
          align = "center" # Center align the box
        )
      } else if (input$num_files == "two") {
        fluidRow(
          column(6, box(
            title = "Gene Network (DIANE - File 1)",
            status = "primary",
            solidHeader = TRUE,
            width = 20, # Increase the width to 20 columns for bigger box
            visNetworkOutput("diane_network_plot_1"),
            align = "center" # Center align the box
          )),
          column(6, box(
            title = "Gene Network (DIANE - File 2)",
            status = "primary",
            solidHeader = TRUE,
            width = 20, # Increase the width to 20 columns for bigger box
            visNetworkOutput("diane_network_plot_2"),
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
          column(12, box(
            title = "Gene Network (SeqNet - File 1)",
            status = "primary",
            solidHeader = TRUE,
            width = 20, # Increase the width to 20 columns for bigger box
            visNetworkOutput("seqnet_network_plot_1"),
            align = "center" # Center align the box
          )),
          column(12, box(
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
      
      if(input$num_files=="one"){
        
        can_download2(TRUE)
        
        print("norm")
        
        
        # Extract the normalized counts from the tcc_object
        normalized_counts <<- TCC::getNormalizedData(tcc_object)
        
        # Output the normalized data to a table in the UI
        output$norm_plot<<-renderPlot({plot(tcc_object)})
        
      }
      else{
        can_download2(TRUE)
        normalized_counts <<- TCC::getNormalizedData(tcc_object)
        
        # Output the normalized data to a table in the UI
        output$norm_plot<<-renderPlot({plot(tcc_object)})
        
        normalized_counts2 <<- TCC::getNormalizedData(tcc_object2)
        
        # Output the normalized data to a table in the UI
        output$norm_plot2<<-renderPlot({plot(tcc_object2)})
      }
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
      can_download3(TRUE)
      if (input$num_files == "one") {
        can_download3(TRUE)
        
        
        print("differential started")
        # Load the data directly from the SeqNet package
        #@ data("abiotic_stresses")
        
        tcc_object <<- DIANE::filter_low_counts(tcc_object, 10*length(abiotic_stresses$conditions))
        fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = abiotic_stresses$conditions)
        print("fit done")
        
        tags <- DIANE::estimateDEGs(fit, reference = "C", perturbation = "H", p.value = 1)
        print("tags done")
        
        output$de_plot<<-renderPlot({ DIANE::draw_DEGs(tags, fdr = 0.01, lfc = 1)})
        
      }
      else{
        can_download3(TRUE)
        tcc_object <<- DIANE::filter_low_counts(tcc_object, 10*length(abiotic_stresses$conditions))
        fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = abiotic_stresses$conditions)
        print("fit done")
        
        tags <- DIANE::estimateDEGs(fit, reference = "C", perturbation = "H", p.value = 1)
        print("tags done")
        
        output$de_plot<<-renderPlot({ DIANE::draw_DEGs(tags, fdr = 0.01, lfc = 1)})
        
        tcc_object2 <<- DIANE::filter_low_counts(tcc_object2, 10*length(abiotic_stresses$conditions))
        fit2 <- DIANE::estimateDispersion(tcc = tcc_object2, conditions = abiotic_stresses$conditions)
        print("fit2 done")
        
        tags2 <- DIANE::estimateDEGs(fit2, reference = "C", perturbation = "H", p.value = 1)
        print("tags2 done")
        
        output$de_plot2<<-renderPlot({ DIANE::draw_DEGs(tags2, fdr = 0.01, lfc = 1)})
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
      if (input$num_files == "one" ) {
        
        can_download1(TRUE)
        
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
      else{
        can_download1(TRUE)
        
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
        
        genes <<- abiotic_stresses$heat_DEGs
        clustering <<- run_coseq(conds = unique(abiotic_stresses$conditions), 
                                 data = abiotic_stresses$normalized_counts, genes = genes, K = 6:9)
        
        #> ****************************************
        #> coseq analysis: Normal approach & arcsin transformation
        #> K = 6 to 9 
        #> Use seed argument in coseq for reproducible results.
        #> ****************************************
        
        print("clustering done")
        output$ebc_plot2<<-renderPlot({DIANE::draw_profiles(data =abiotic_stresses$normalized_counts, clustering$membership, conds = unique(abiotic_stresses$conditions), k = 3) })
        
        
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
    if (is.null(input$tool_choice)) {
      showModal(modalDialog(
        title = "Tool Selection Required",
        "Please select a network generation tool from the Tools menu."
      ))
    }
    else if  ("seqnet" %in% input$tool_choice){
      if (fileUploaded()) {
        can_download4(TRUE)
        if (input$num_files == "one" ){
          can_download4(TRUE)
          # Reset the flags at the beginning
          first_click(FALSE)
          user_selected_node(FALSE)
          
          data("abiotic_stresses")
          edges <- abiotic_stresses$heat_DEGs_regulatory_links
          edges <- edges[-(19:25), ]
          
          row_names <- rownames(edges)
          column_names <- colnames(edges)
          common_genes <- intersect(row_names, column_names)
          
          edges <- edges[common_genes, common_genes]
          symmetric_edges <- pmax(edges, t(edges))
          
          threshold <- 0  # setting threshold
          binary_edges <- ifelse(symmetric_edges > threshold, 1, 0)
          
          # create the network with the adjacency matrix
          nw <- create_network_from_adjacency_matrix(binary_edges)
          
          nodes_df <- data.frame(
            id = nw$node_names,
            label = nw$node_names,
            name = nw$node_names,
            color = "orange",
            size = 30,
            value = runif(nrow(binary_edges))
          )
          
          nodes_df_reactive(nodes_df)
          
          edges_data <- nw$modules[[1]]$edges
          edges_df <- data.frame(
            from = edges_data[, 1],
            to = edges_data[, 2]
          )
          
          # if there are corresponding node names, it can map them
          edges_df$from <- nw$node_names[edges_df$from]
          edges_df$to <- nw$node_names[edges_df$to]
          
          # Add a color attribute to edges
          edges_df$color <- "grey"
          
          
          edges_df_reactive(edges_df)
          observe({
            # Assuming nodes_df has an 'id' column with the node ids
            updateSelectInput(session, "nodeSelector",
                              choices = nodes_df$id)
          })
          observeEvent(input$nodeSelector, {
            local <- first_click()
            if(local){
              print("here")
              print(input$nodeSelector)
              if(input$nodeSelector !="AT1G03840" ){
                connected_nodes_node <- get_connected_nodes(input$nodeSelector, edges_df)
                
                showModal(modalDialog(
                  title = "Node Information",
                  paste("Clicked Node:", input$nodeSelector, 
                        "\nConnected Nodes:", paste(connected_nodes_node, collapse = ", "))
                ))
              }
            }
            else(
              first_click(TRUE)
            )
          })
          get_connected_nodes <- function(clicked_node, edges_df) {
            # Find rows in edges_df where the clicked node is either 'from' or 'to'
            connected <- subset(edges_df, from == clicked_node | to == clicked_node)
            # Return the unique node names connected to the clicked node
            unique(c(connected$from, connected$to))
          }
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
              visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
              
              visInteraction(dragNodes = TRUE) %>%
              visEvents(click = "function(nodes){
    Shiny.setInputValue('interactive_graph_nodes', nodes, {priority: 'event'});
  }")
          })
          shinyjs::hide("generate_button")
        }
        else if (input$num_files == "two"){
          can_download4(TRUE)
          # Reset the flags at the beginning
          first_click(FALSE)
          user_selected_node(FALSE)
          
          data("abiotic_stresses")
          edges <- abiotic_stresses$heat_DEGs_regulatory_links
          edges <- edges[-(19:25), ]
          
          row_names <- rownames(edges)
          column_names <- colnames(edges)
          common_genes <- intersect(row_names, column_names)
          
          edges <- edges[common_genes, common_genes]
          symmetric_edges <- pmax(edges, t(edges))
          
          threshold <- 0  # setting threshold
          binary_edges <- ifelse(symmetric_edges > threshold, 1, 0)
          
          # create the network with the adjacency matrix
          nw <- create_network_from_adjacency_matrix(binary_edges)
          
          nodes_df <- data.frame(
            id = nw$node_names,
            label = nw$node_names,
            name = nw$node_names,
            color = "orange",
            size = 30,
            value = runif(nrow(binary_edges))
          )
          
          nodes_df_reactive(nodes_df)
          
          edges_data <- nw$modules[[1]]$edges
          edges_df <- data.frame(
            from = edges_data[, 1],
            to = edges_data[, 2]
          )
          
          # if there are corresponding node names, it can map them
          edges_df$from <- nw$node_names[edges_df$from]
          edges_df$to <- nw$node_names[edges_df$to]
          
          # Add a color attribute to edges
          edges_df$color <- "grey"
          
          
          edges_df_reactive(edges_df)
          observe({
            # Assuming nodes_df has an 'id' column with the node ids
            updateSelectInput(session, "nodeSelector",
                              choices = nodes_df$id)
          })
          observeEvent(input$nodeSelector, {
            local <- first_click()
            if(local){
              print("here")
              print(input$nodeSelector)
              if(input$nodeSelector !="AT1G03840" ){
                connected_nodes_node <- get_connected_nodes(input$nodeSelector, edges_df)
                
                showModal(modalDialog(
                  title = "Node Information",
                  paste("Clicked Node:", input$nodeSelector, 
                        "\nConnected Nodes:", paste(connected_nodes_node, collapse = ", "))
                ))
              }
            }
            else(
              first_click(TRUE)
            )
          })
          get_connected_nodes <- function(clicked_node, edges_df) {
            # Find rows in edges_df where the clicked node is either 'from' or 'to'
            connected <- subset(edges_df, from == clicked_node | to == clicked_node)
            # Return the unique node names connected to the clicked node
            unique(c(connected$from, connected$to))
          }
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
          output$seqnet_network_plot_1 <- renderVisNetwork({
            visNetwork(nodes_df, edges_df, width = "100%", height = "8000px") %>%
              visIgraphLayout() %>%
              visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
              visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
              
              visInteraction(dragNodes = TRUE) %>%
              visEvents(click = "function(nodes){
    Shiny.setInputValue('interactive_graph_nodes', nodes, {priority: 'event'});
  }")
            
          })
          output$seqnet_network_plot_2 <- renderVisNetwork({
            visNetwork(nodes_df, edges_df, width = "100%", height = "8000px") %>%
              visIgraphLayout() %>%
              visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
              visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
              
              visInteraction(dragNodes = TRUE) %>%
              visEvents(click = "function(nodes){
    Shiny.setInputValue('interactive_graph_nodes', nodes, {priority: 'event'});
  }")
            
          })
          shinyjs::hide("generate_button")
          
        }
        
      }
      else{
        showModal(
          modalDialog(
            title = "Error",
            "Please upload file.",
            easyClose = TRUE
          )
        ) 
      }
    }
    if ("diane" %in% input$tool_choice){
      if (fileUploaded()) {
        if (input$num_files == "one" ){
          expression_data <- file1
          
          
          regulators <- targets <- row_ids <- rownames(file1)
          print(abiotic_stresses$conditions)
          network_matrix <- network_inference(expression_data, 
                                              conds =abiotic_stresses$conditions,
                                              targets = targets, 
                                              regressors = regulators,
                                              nTrees = 10, 
                                              nCores = 4)
          
          # Step 4: Threshold the Network
          # Thresholding to keep significant edges
          network <- network_thresholding(network_matrix, n_edges = 100)  # Adjust n_edges as needed
          
          # Step 5: Prepare Data for Visualization
          network_data_for_vis <- network_data(network, regulators)
          
          # Step 6: Visualize the Network
          nodes <- network_data_for_vis$nodes
          edges <- network_data_for_vis$edges
          
          
          get_connected_nodes <- function(clicked_node, edges_df) {
            # Find rows in edges_df where the clicked node is either 'from' or 'to'
            connected <- subset(edges_df, from == clicked_node | to == clicked_node)
            # Return the unique node names connected to the clicked node
            unique(c(connected$from, connected$to))
          }
          observeEvent(input$interactive_graph_nodes2, {
            if(!identical(input$interactive_graph_nodes2$nodes, list())) {
              clicked_node <- input$interactive_graph_nodes2$nodes[1]
              print(clicked_node)
              connected_nodes_node <- get_connected_nodes(clicked_node, edges)
              
              showModal(modalDialog(
                title = "Node Information",
                paste("Clicked Node:", clicked_node, 
                      "\nConnected Nodes:", paste(connected_nodes_node, collapse = ", "))
              ))
            }
          }, ignoreNULL = TRUE)
          output$diane_network_plot <- renderVisNetwork({
            visNetwork(nodes, edges, width = "100%", height = "8000px") %>%
              visIgraphLayout() %>%
              visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
              visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
              
              visInteraction(dragNodes = TRUE) %>%
              visEvents(click = "function(nodes){
    Shiny.setInputValue('interactive_graph_nodes2', nodes, {priority: 'event'});
  }")
          })
          print("here")
        }
        else if (input$num_files == "two"){
          expression_data <- file1
          
          
          regulators <- targets <- row_ids <- rownames(file1)
          print(abiotic_stresses$conditions)
          network_matrix <- network_inference(expression_data, 
                                              conds =abiotic_stresses$conditions,
                                              targets = targets, 
                                              regressors = regulators,
                                              nTrees = 10, 
                                              nCores = 4)
          
          # Step 4: Threshold the Network
          # Thresholding to keep significant edges
          network <- network_thresholding(network_matrix, n_edges = 100)  # Adjust n_edges as needed
          
          # Step 5: Prepare Data for Visualization
          network_data_for_vis <- network_data(network, regulators)
          
          # Step 6: Visualize the Network
          nodes <- network_data_for_vis$nodes
          edges <- network_data_for_vis$edges
          get_connected_nodes <- function(clicked_node, edges_df) {
            # Find rows in edges_df where the clicked node is either 'from' or 'to'
            connected <- subset(edges_df, from == clicked_node | to == clicked_node)
            # Return the unique node names connected to the clicked node
            unique(c(connected$from, connected$to))
          }
          observeEvent(input$interactive_graph_nodes2, {
            if(!identical(input$interactive_graph_nodes2$nodes, list())) {
              clicked_node <- input$interactive_graph_nodes2$nodes[1]
              print(clicked_node)
              connected_nodes_node <- get_connected_nodes(clicked_node, edges)
              
              showModal(modalDialog(
                title = "Node Information",
                paste("Clicked Node:", clicked_node, 
                      "\nConnected Nodes:", paste(connected_nodes_node, collapse = ", "))
              ))
            }
          }, ignoreNULL = TRUE)
          output$diane_network_plot_1 <- renderVisNetwork({
            visNetwork(nodes, edges, width = "100%", height = "8000px") %>%
              visIgraphLayout() %>%
              visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
              visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
              
              visInteraction(dragNodes = TRUE) %>%
              visEvents(click = "function(nodes){
    Shiny.setInputValue('interactive_graph_nodes2', nodes, {priority: 'event'});
  }")
          })
          
          
          expression_data2 <- file2
          
          
          regulators2 <- targets2 <- row_ids <- rownames(file2)
          print(abiotic_stresses$conditions)
          network_matrix <- network_inference(expression_data2, 
                                              conds =abiotic_stresses$conditions,
                                              targets = targets2, 
                                              regressors = regulators2,
                                              nTrees = 10, 
                                              nCores = 4)
          
          # Step 4: Threshold the Network
          # Thresholding to keep significant edges
          network2 <- network_thresholding(network_matrix, n_edges = 100)  # Adjust n_edges as needed
          
          # Step 5: Prepare Data for Visualization
          network_data_for_vis2 <- network_data(network, regulators)
          
          # Step 6: Visualize the Network
          nodes2 <- network_data_for_vis$nodes
          edges2 <- network_data_for_vis$edges
          print(edges2)
          print(nodes2)
          
          get_connected_nodes <- function(clicked_node, edges_df) {
            # Find rows in edges_df where the clicked node is either 'from' or 'to'
            connected <- subset(edges_df, from == clicked_node | to == clicked_node)
            # Return the unique node names connected to the clicked node
            unique(c(connected$from, connected$to))
          }
          observeEvent(input$interactive_graph_nodes3, {
            if(!identical(input$interactive_graph_nodes3$nodes, list())) {
              clicked_node <- input$interactive_graph_nodes3$nodes[1]
              print(clicked_node)
              connected_nodes_node <- get_connected_nodes(clicked_node, edges2)
              
              showModal(modalDialog(
                title = "Node Information",
                paste("Clicked Node:", clicked_node, 
                      "\nConnected Nodes:", paste(connected_nodes_node, collapse = ", "))
              ))
            }
          }, ignoreNULL = TRUE)
          output$diane_network_plot_2 <- renderVisNetwork({
            visNetwork(nodes, edges, width = "100%", height = "8000px") %>%
              visIgraphLayout() %>%
              visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
              visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
              
              visInteraction(dragNodes = TRUE) %>%
              visEvents(click = "function(nodes){
    Shiny.setInputValue('interactive_graph_nodes3', nodes, {priority: 'event'});
  }")
          })
        }
      }
      else{
        showModal(
          modalDialog(
            title = "Error",
            "Please upload file.",
            easyClose = TRUE
          )
        )
      }
    }
  })
  
  
  
  
  
  
  #end of diane 
  
  
  
  observeEvent(input$change_color, {
    if ("seqnet" %in% input$tool_choice) {
      if (input$num_files == "one"){
        if (!is.null(nodes_df_reactive())) {
          print("hiding the button")
          shinyjs::hide("change_color")
          shinyjs::show("revert_color")
          updated_nodes_df <- nodes_df_reactive()
          updated_nodes_df$color <- "blue"
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
              visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
              
              visEvents(click = "function(nodes){
    Shiny.setInputValue('interactive_graph_nodes', nodes, {priority: 'event'});
  }")
            
          })
        }
      }
      else if (input$num_files == "two"){
        if (!is.null(nodes_df_reactive())) {
          print("hiding the button")
          shinyjs::hide("change_color")
          shinyjs::show("revert_color")
          updated_nodes_df <- nodes_df_reactive()
          updated_nodes_df$color <- "blue"
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
          
          output$seqnet_network_plot_2 <- renderVisNetwork({
            visNetwork(nodes_df_reactive(), edges_df_reactive(), width = "100%", height = "8000px") %>%
              visIgraphLayout() %>%
              visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
              visInteraction(dragNodes = TRUE) %>%
              visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
              
              visEvents(click = "function(nodes){
    Shiny.setInputValue('interactive_graph_nodes', nodes, {priority: 'event'});
  }")
            
          })
          output$seqnet_network_plot_1 <- renderVisNetwork({
            visNetwork(nodes_df_reactive(), edges_df_reactive(), width = "100%", height = "8000px") %>%
              visIgraphLayout() %>%
              visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
              visInteraction(dragNodes = TRUE) %>%
              visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
              
              visEvents(click = "function(nodes){
    Shiny.setInputValue('interactive_graph_nodes', nodes, {priority: 'event'});
  }")
            
          })
        }
      }
    }
  })
  
  observeEvent(input$revert_color, {
    if ("seqnet" %in% input$tool_choice) {
      if (input$num_files == "one"){
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
              visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
              
              visEvents(click = "function(nodes){
            Shiny.setInputValue('interactive_graph_nodes', nodes, {priority: 'event'});
          }")
          })
        }
      }
      else if (input$num_files == "two"){
        if (!is.null(nodes_df_reactive())) {
          updated_nodes_df <- nodes_df_reactive()
          updated_nodes_df$color <- "orange"  # Revert color to orange
          nodes_df_reactive(updated_nodes_df)
          shinyjs::show("change_color")
          shinyjs::hide("revert_color")
          
          # Update the network plot
          output$seqnet_network_plot_1 <- renderVisNetwork({
            visNetwork(nodes_df_reactive(), edges_df_reactive(), width = "100%", height = "8000px") %>%
              visIgraphLayout() %>%
              visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
              visInteraction(dragNodes = TRUE) %>%
              visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
              
              visEvents(click = "function(nodes){
            Shiny.setInputValue('interactive_graph_nodes', nodes, {priority: 'event'});
          }")
          })
          output$seqnet_network_plot_2 <- renderVisNetwork({
            visNetwork(nodes_df_reactive(), edges_df_reactive(), width = "100%", height = "8000px") %>%
              visIgraphLayout() %>%
              visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
              visInteraction(dragNodes = TRUE) %>%
              visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
              
              visEvents(click = "function(nodes){
            Shiny.setInputValue('interactive_graph_nodes', nodes, {priority: 'event'});
          }")
          })
        }
      }
    }
  })
  
  
  
  nodes_df <- reactive({
    data.frame(id = 1:10,
               label = paste("Label", 1:10),
               color = rep("orange", 10),
               shape = rep("circle", 10))
  })
  
  # Example reactive dataset for edges, replace with your actual code to generate edges dataframe
  edges_df <- reactive({
    data.frame(from = sample(1:10, 10, replace = TRUE),
               to = sample(1:10, 10, replace = TRUE))
  })
  
  # Example observer for changing node color intensity
  observeEvent(input$nodeColorIntensity, {
    updated_nodes <- nodes_df() # Get the current nodes data
    updated_nodes$color <- adjustcolor(updated_nodes$color, alpha.f = input$nodeColorIntensity)
    
    output$seqnet_network_plot <- renderVisNetwork({
      visNetwork(updated_nodes, edges_df()) %>%
        visIgraphLayout() %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
    })
  })
  
  network_density <- reactiveVal(0)
  
  # When the density button pressed The network density gives a quick overview of how connected the network is 
  observeEvent(input$showDensity, {
    
    # Calculate network density
    g <- graph_from_data_frame(edges_df(), directed = FALSE)
    density_val <- graph.density(g)
    network_density(density_val)
  })
  
  # Display network density when the button is clicked
  output$networkDensity <- renderText({
    input$showDensity
    paste("Network Density:", network_density())
  })
  
  average_degree <- reactiveVal(0)
  
  observeEvent(input$showAvgDegree, {
    
    # Calculate the average degree after generating the network The average degree is a measure of the average number of connections per node in your network
    g <- igraph::graph_from_data_frame(edges_df(), directed = FALSE)
    degree_vals <- igraph::degree(g)
    avg_degree <- mean(degree_vals)
    average_degree(avg_degree)  # Store the average degree
  })
  
  # Display the average degree when the 'Show Average Node Degree' button is clicked
  output$averageNodeDegree <- renderText({
    req(input$showAvgDegree)  # This ensures the calculation is done only after the button is clicked
    paste("Average Node Degree:", average_degree())
  })
  
  # ... (rest of your server logic)
  observeEvent(input$tabs, {
    if (input$tabs == "download") {
      print(input$tool_choice)
      print( class(input$tool_choice))
      
      
      
      # Code to execute when the Download tab is selected
      # Example: print a message to the R console
      shinyjs::hide("download_data")
      shinyjs::hide("download_data2")
      
      shinyjs::hide("download_data3")
      shinyjs::hide("download_data4")
      if (can_download2()==TRUE){
        shinyjs::show("download_data")
      }
      
      if (can_download3()==TRUE){
        shinyjs::show("download_data2")
      }
      
      if (can_download1()==TRUE){
        shinyjs::show("download_data3")
      }
      if (can_download4()==TRUE){
        shinyjs::show("download_data4")
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
