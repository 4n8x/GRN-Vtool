#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

server <- function(input, output, session) {
  
  # Initialize an empty igraph graph object
  g <- reactiveVal(NULL)
  
  
  # Generate tool selection input based on user choices
  output$tool_input <- renderUI({
    if (input$num_files == "one") {
      checkboxGroupInput("tool_choice", label = "Choose Network Generation Tool:",
                         choices = list("DIANE" = "diane", "SeqNet" = "seqnet"),
                         selected = c("diane", "seqnet"))
    } else {
      radioButtons("tool_choice", label = "Choose Network Generation Tool:",
                   choices = list("DIANE" = "diane", "SeqNet" = "seqnet"),
                   selected = "diane")
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
          plotOutput("seqnet_network_plot"),
          align = "center" # Center align the box
        )
      } else if (input$num_files == "two") {
        fluidRow(
          column(6, box(
            title = "Gene Network (SeqNet - File 1)",
            status = "primary",
            solidHeader = TRUE,
            width = 20, # Increase the width to 20 columns for bigger box
            plotOutput("seqnet_network_plot_1"),
            align = "center" # Center align the box
          )),
          column(6, box(
            title = "Gene Network (SeqNet - File 2)",
            status = "primary",
            solidHeader = TRUE,
            width = 20, # Increase the width to 20 columns for bigger box
            plotOutput("seqnet_network_plot_2"),
            align = "center" # Center align the box
          ))
        )
      }
    }
  })
  
  
  
  
  observeEvent(input$generate_button, {
    if ("seqnet" %in% input$tool_choice) {
      p <- 100  # Number of nodes (genes) in the network.
      n <- 100  # Number of samples to generate.
      network <- random_network(p)  # Create a random network of size p.
      network <- gen_partial_correlations(network)  # Generate weights.
      df_ref <- SeqNet::reference$rnaseq  # The reference dataset
      
      g(network)  # Assign the generated network to the 'g' object
      
      if (input$num_files == "one") {
        output$seqnet_network_plot <- renderPlot({
          if (!is.null(g())) {
            plot(g(), layout = layout_with_fr)
          }
        })
      } else if (input$num_files == "two") {
        output$seqnet_network_plot_1 <- renderPlot({
          if (!is.null(g())) {
            plot(g(), layout = layout_with_fr)
          }
        })
        output$seqnet_network_plot_2 <- renderPlot({
          if (!is.null(g())) {
            plot(g(), layout = layout_with_fr)
          }
        })
      }
    }
  })
  
  #DIANE code
  
  
  observeEvent(input$norm_button, {
    if ("diane" %in% input$tool_choice) {
      
      
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
      
      
      
      
    
    }
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
      
      
      
      
      output$network_plot<<-renderUI({DIANE::draw_network(data$nodes, data$edges)})
      
      
      
      
      
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
  
}
  
