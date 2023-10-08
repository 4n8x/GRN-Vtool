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




library(DIANE)




function(input, output, session) {
  

observeEvent(input$tool_choice ,{


  if(input$tool_choice == "seqNet"){
   
  

#---


  
  
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
    
    
    output$network_plot<-renderui({plot(g)})
    
    output$network_plot_2<-renderUI({plot(g)})
    
    
    
    
  })##end network buttonn
  
  
  
  
  


  }#end of seqnet

  if(input$tool_choice == "DIANE"){


  
  
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
    

    
    
    
  })##end norm button
  
  
  observeEvent(input$de_button, {
    
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
    
  })##end eda button
  
  observeEvent(input$ebc_button, {
    
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
    
    
    
  })##end clustering button
  
  observeEvent(input$network_button, {
    
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
    

    
  })##end network button
  
}#end of diane


}#end of tool choice


)}#end session function

