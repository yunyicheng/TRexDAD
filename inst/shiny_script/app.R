library(shiny)
library(shinythemes)
library(ggplot2)

# --- SECTION: Setup and Global Variable Initialization -----------------


# Read overhang fidelity chart
overhang_fidelity_path <- system.file("extdata", "overhang_fidelity.csv", 
                                      package = "TRexDAD")
OVERHANG_FIDELITY <- utils::read.csv(overhang_fidelity_path)

# Rad-27 gene for example usage
RAD_27 <- "AATATGGGTATTAAAGGTTTGAATGCAATTATATCGGAACATGTTCCCTCTGCTATCAGGAAAAGCGATATCAAGAGCTTTTTTGGCAGAAAGGTTGCCATCGATGCCTCTATGTCTCTATATCAGTTTTTAATTGCTGTAAGACAGCAAGACGGTGGGCAGTTGACCAATGAAGCCGGTGAAACAACGTCACACTTGATGGGTATGTTTTATAGGACACTGAGAATGATTGATAACGGTATCAAGCCTTGTTATGTCTTCGACGGCAAACCTCCAGATTTGAAATCTCATGAGTTGACAAAGCGGTCTTCAAGAAGGGTGGAAACAGAAAAAAAACTGGCAGAGGCAACAACAGAATTGGAAAAGATGAAGCAAGAAAGAAGATTGGTGAAGGTTTCAAAAGAGCATAATGAAGAAGCCCAAAAATTACTAGGACTAATGGGAATCCCATATATAATAGCGCCAACGGAAGCTGAGGCTCAATGTGCTGAGTTGGCAAAGAAGGGAAAGGTGTATGCCGCAGCAAGTGAAGATATGGACACACTCTGTTATAGAACACCCTTCTTGTTGAGACATTTGACTTTTTCAGAGGCCAAGAAGGAACCGATTCACGAAATAGATACTGAATTAGTTTTGAGAGGACTCGACTTGACAATAGAGCAGTTTGTTGATCTTTGCATAATGCTTGGTTGTGACTACTGTGAAAGCATCAGAGGTGTTGGTCCAGTGACAGCCTTAAAATTGATAAAAACGCATGGATCCATCGAAAAAATCGTGGAGTTTATTGAATCTGGGGAGTCAAACAACACTAAATGGAAAATCCCAGAAGACTGGCCTTACAAACAAGCAAGAATGCTGTTTCTTGACCCTGAAGTTATAGATGGTAACGAAATAAACTTGAAATGGTCGCCACCAAAGGAGAAGGAACTTATCGAGTATTTATGTGATGATAAGAAATTCAGTGAAGAAAGAGTTAAATCTGGTATATCAAGATTGAAAAAAGGCTTGAAATCTGGCATTCAGGGTAGGTTAGATGGGTTCTTCCAAGTGGTGCCTAAGACAAAGGAACAGCTGGCTGCTGCGGCGAAAAGAGCACAAGAAAATAAAAAATTGAACAAAAATAAGAATAAAGTCACAAAGGGAAGAAGATGAGGG"

# --- SECTION: Shiny App Implementation -----------------
# Define UI
ui <- fluidPage(
    theme = shinythemes::shinytheme("flatly"),
    
    tags$head(
        tags$style(HTML("
            body { 
                background-color: #f3f3f3;
                font-family: 'Arial', sans-serif;
            }
            .shiny-output-error { 
                color: #a94442; 
                font-weight: bold; 
            }
            .well {
                background-color: #fff; 
                border-radius: 15px; 
                padding: 20px; 
                box-shadow: 0 4px 8px 0 rgba(0,0,0,0.2);
            }
            .shiny-input-container {
                margin-bottom: 15px;
            }
            #run_optimization {
                width: 100%;
                font-size: 16px;
                font-weight: bold;
                border: none;
                color: white;
                background-color: #3498DB;
                padding: 14px 20px;
                margin: 8px 0;
                border-radius: 4px;
                cursor: pointer;
            }
            #run_optimization:hover {
                background-color: #286090;
            }
            .panel-title {
                font-size: 24px;
                font-weight: bold;
                color: #333;
                margin-bottom: 20px;
            }
            .panel-body {
                padding: 35px;
            }
        ")),
        tags$link(href = 
            "https://fonts.googleapis.com/css?family=Open+Sans:400,600,700", 
            rel = "stylesheet")
    ),
    
    titlePanel("TRexDAD Pipeline: 
               Assembly Design for Tile Region Exchange Mutagenesis", 
               windowTitle = "TRexDAD Assembly Design"),
    
    sidebarLayout(
        sidebarPanel(
            textInput("gene_sequence", "Enter your gene sequence", 
                      value = RAD_27),
            numericInput("max_iter", "Maximum Iterations", value = 10),
            numericInput("scan_rate", "Scan Rate", value = 7),
            actionButton("run_optimization", "Run Optimization", 
                         class = "btn-primary")
        ),
        
        mainPanel(
            tabsetPanel(
                tabPanel("Plot", plotOutput("score_plot")),
                tabPanel("Assembly Position", verbatimTextOutput("final_pos")),
                tabPanel("Optimal Tiles", 
                         verbatimTextOutput("optimal_tiles_output")),
                tabPanel("Gene Codons", verbatimTextOutput("codons_output")),
                tabPanel("Oligo Cost", verbatimTextOutput("oligo_cost_output")),
                tabPanel("Overhang Fidelity", 
                         verbatimTextOutput("overhang_fidelity_output"))
            )
        )
    )
)

# Define server logic
server <- function(input, output) {
    observeEvent(input$run_optimization, {
        # Input validation
        
        # Verify that the sequence only contains valid nucleotides
        if(!grepl("^[ACGT]+$", input$gene_sequence)) {
            stop("Invalid gene sequence: sequence should only contain A, C, G, and T.")
        }
        
        # Check if the length of the sequence is a multiple of 3
        if (nchar(input$gene_sequence) %% 3 != 0) {
            stop("Invalid gene sequence length: length should be a multiple of 3.")
        }
        
        # Ensure 'ATG' is the second codon
        if (substr(input$gene_sequence, 4, 6) != "ATG") {
            stop("Invalid sequence: 'ATG' is not the second codon.")
        }
        
        # Convert the gene sequence into codons
        gene_codons <- split_into_codons(input$gene_sequence)
        
        # Check if the number of codons is 60 or less
        # an ORF would typically need to be at least 60-90 nucleotides long to 
        # encode a functional protein of minimal size.
        if (length(gene_codons) < 60) {
            stop("Target gene produces 60 or fewer codons, which is insufficient 
             for the optimization process.")
        }
        
        # Calculate optimal tiles
        optimal_tiles_result <- calculate_optimal_tiles(length(gene_codons))
        
        # Calculate oligo cost
        oligo_cost_result <- oligo_cost(optimal_tiles_result$optimal_num_tiles, 
                                        length(gene_codons))
        
        
        # Execute the optimization process and plot results
        result <- execute_and_plot(input$gene_sequence, 
                                   input$max_iter, input$scan_rate)
        
        
        # Output the results to the UI
        output$overhang_fidelity_output <- renderPrint({
            OVERHANG_FIDELITY
        })
        output$codons_output <- renderPrint({
            gene_codons
        })
        
        output$oligo_cost_output <- renderPrint({
            oligo_cost_result
        })
        
        output$optimal_tiles_output <- renderPrint({
            optimal_tiles_result
        })
        
        output$score_plot <- renderPlot({
            result$plot
        })
        
        output$final_pos <- renderPrint({
            result$final_position
        })
    })
}

# Run the application
shinyApp(ui = ui, server = server)