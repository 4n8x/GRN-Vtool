# Load required packages
if (!require("shiny"))
  install.packages("shiny")

if (!require("edgeR"))
  install.packages("edgeR")

library(shiny)
library(edgeR)

# Define UI for the Shiny app
ui <- fluidPage(
  titlePanel("Normalization App with edgeR"),
  sidebarLayout(
    sidebarPanel(
      fileInput("data_file", "Choose a CSV file containing gene expression data:"),
      actionButton("normalize_button", "Normalize Data"),
      downloadButton("download_normalized_data", "Download Normalized Data")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Before Normalization", plotOutput("before_norm_plot")),
        tabPanel("After Normalization", plotOutput("after_norm_plot")),
        tabPanel("Normalized Data", tableOutput("normalized_table"))
      )
    )
  )
)

# Define server for the Shiny app
server <- function(input, output) {
  data <- reactive({
    req(input$data_file)
    raw_data <- read.csv(input$data_file$datapath, row.names = 1, header = TRUE, sep = ",")
    return(raw_data)
  })
  
  normalized_data <- eventReactive(input$normalize_button, {
    # Load data
    counts <- data()
    
    # Perform normalization using edgeR's TMM (Trimmed Mean of M-values) method
    normalized_counts <- calcNormFactors(counts)
    
    return(normalized_counts)
  })
  
  output$before_norm_plot <- renderPlot({
    # Plot data before normalization (e.g., boxplot or heatmap)
    # Use 'data()' to access the uploaded data
  })
  
  output$after_norm_plot <- renderPlot({
    # Plot data after normalization (e.g., boxplot or heatmap)
    # Use 'normalized_data()' to access the normalized data
  })
  
  output$normalized_table <- renderTable({
    # Display the normalized data in a table format
    # Use 'normalized_data()' to access the normalized data
    normalized_data()
  })
  
  # Create a downloadable CSV file for the normalized data
  output$download_normalized_data <- downloadHandler(
    filename = function() {
      paste("normalized_data", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(normalized_data(), file)
    }
  )
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
