library(shiny)
library(Biostrings)
library(forstringr)

# Assume OVERHANG_FIDELITY and RAD_27 are loaded or defined within the app

ui <- fluidPage(
    titlePanel("TRexDAD Pipeline: Assembly Design for Tile Region
               Exchange Mutagenesis"),
    
    sidebarLayout(
        sidebarPanel(
            textInput("gene_sequence", "Enter your gene sequence", value = RAD_27),
            numericInput("max_iter", "Maximum Iterations", value = 10),
            numericInput("scan_rate", "Scan Rate", value = 7),
            actionButton("run_optimization", "Run Optimization")
        ),
        
        mainPanel(
            plotOutput("score_plot"),
            verbatimTextOutput("optimal_tiles_output"),
            verbatimTextOutput("codons_output"),
            verbatimTextOutput("oligo_cost_output")
        )
    )
)

server <- function(input, output) {
    observeEvent(input$run_optimization, {
        # Input validation can be added here
        
        
        # Convert the gene sequence into codons
        gene_codons <- split_into_codons(input$gene_sequence)
        num_codons <- length(gene_codons)
        
        # Calculate optimal tiles
        optimal_tiles_result <- calculate_optimal_tiles(num_codons)
        
        # Calculate oligo cost
        oligo_cost_result <- oligo_cost(optimal_tiles_result$optimal_num_tiles, 
                                        num_codons)
        
        # Execute the optimization process and plot results
        # This would be a reactive function if it affects outputs
        results <- execute_and_plot(input$gene_sequence, input$max_iter, 
                                    input$scan_rate)
        
        # Output the results to the UI
        output$codons_output <- renderPrint({
            print(paste("gene_codons=", list(gene_codons)))
        })
        
        output$oligo_cost_output <- renderPrint({
            print(paste("oligo_cost_result =", oligo_cost_result))
        })
        
        output$optimal_tiles_output <- renderPrint({
            print(paste("optimal_tiles_result =", optimal_tiles_result))
        })
        
        output$score_plot <- renderPlot({
            results$plot
        })
    })
}

shinyApp(ui = ui, server = server)