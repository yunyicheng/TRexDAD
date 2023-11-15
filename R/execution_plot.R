source("R/calculate_scores.R")
source("R/optimize_position.R")
source("R/process_data.R")
library(readxl)
# utils::globalVariables(c(overhang_fidelity, target_gene, num_tiles_global, 
#                          length_tiles_global))
#' Execute the Optimization Process and Plot Results
#'
#' This function performs an optimization process to determine the optimal tile
#' positions in a gene sequence and plots the progression of the score over iterations.
#' It reads the overhang fidelity data, computes the optimal number of tiles, and iteratively
#' improves the positions, plotting the change in score with respect to number of iterations at the end.
#'
#' @param iteration_max The maximum number of iterations for the optimization process (default: 30).
#' @param scan_rate The range within which to scan for optimizing tile positions (default: 7).
#'
#' @return This function does not return a value; it generates a plot showing the score
#'         optimization over iterations.
#'
#' @import readxl
#' 
#' @export
execute_and_plot <- function(iteration_max=30, scan_rate=7) {
    
    # Read overhang fidelity chart
    overhang_fidelity <<- read_excel("inst/extdata/overhang_fidelity.xlsx",
                                    sheet = "S1 Table. BsaI-HFv2")
    
    target_gene <<- "AATATGGGTATTAAAGGTTTGAATGCAATTATATCGGAACATGTTCCCTCTGCTATCAGGAAAAGCGATATCAAGAGCTTTTTTGGCAGAAAGGTTGCCATCGATGCCTCTATGTCTCTATATCAGTTTTTAATTGCTGTAAGACAGCAAGACGGTGGGCAGTTGACCAATGAAGCCGGTGAAACAACGTCACACTTGATGGGTATGTTTTATAGGACACTGAGAATGATTGATAACGGTATCAAGCCTTGTTATGTCTTCGACGGCAAACCTCCAGATTTGAAATCTCATGAGTTGACAAAGCGGTCTTCAAGAAGGGTGGAAACAGAAAAAAAACTGGCAGAGGCAACAACAGAATTGGAAAAGATGAAGCAAGAAAGAAGATTGGTGAAGGTTTCAAAAGAGCATAATGAAGAAGCCCAAAAATTACTAGGACTAATGGGAATCCCATATATAATAGCGCCAACGGAAGCTGAGGCTCAATGTGCTGAGTTGGCAAAGAAGGGAAAGGTGTATGCCGCAGCAAGTGAAGATATGGACACACTCTGTTATAGAACACCCTTCTTGTTGAGACATTTGACTTTTTCAGAGGCCAAGAAGGAACCGATTCACGAAATAGATACTGAATTAGTTTTGAGAGGACTCGACTTGACAATAGAGCAGTTTGTTGATCTTTGCATAATGCTTGGTTGTGACTACTGTGAAAGCATCAGAGGTGTTGGTCCAGTGACAGCCTTAAAATTGATAAAAACGCATGGATCCATCGAAAAAATCGTGGAGTTTATTGAATCTGGGGAGTCAAACAACACTAAATGGAAAATCCCAGAAGACTGGCCTTACAAACAAGCAAGAATGCTGTTTCTTGACCCTGAAGTTATAGATGGTAACGAAATAAACTTGAAATGGTCGCCACCAAAGGAGAAGGAACTTATCGAGTATTTATGTGATGATAAGAAATTCAGTGAAGAAAGAGTTAAATCTGGTATATCAAGATTGAAAAAAGGCTTGAAATCTGGCATTCAGGGTAGGTTAGATGGGTTCTTCCAAGTGGTGCCTAAGACAAAGGAACAGCTGGCTGCTGCGGCGAAAAGAGCACAAGAAAATAAAAAATTGAACAAAAATAAGAATAAAGTCACAAAGGGAAGAAGATGAGGG"
    target_gene <<- split_into_codons(target_gene)
    # Determine optimal number of tiles and codons
    num_codons <<- length(target_gene) - 2
    num_tiles <<- 1:50
    
    # Define the range for the number of tiles
    lower_bound <- 1
    upper_bound <- 50
    
    # Perform optimization
    optimization_result <- optimize(f = function(x) oligo_cost(x, num_codons),
                                    interval = c(lower_bound, upper_bound))
    
    # Extract the results
    num_tiles_global <<- round(optimization_result$minimum)
    length_tiles_global <<- round(num_codons / num_tiles_global)
    
    # utils::globalVariables(c(overhang_fidelity, target_gene, num_tiles_global, 
    #                          length_tiles_global))
    # utils::suppressForeignCheck(c(overhang_fidelity, target_gene, num_tiles_global, 
    #                               length_tiles_global))
    # Print results
    cat("Optimal number of tiles =", num_tiles_global, "\n")
    cat("Optimal length of tiles =", length_tiles_global, "\n")
    
    # Initialize cassette positions
    pos <- seq(3, length(target_gene) - 3, by = length_tiles_global)
    pos <- c(pos, length(target_gene) - 3)
    print(pos)
    # Obtain initial scores
    curr_score <- calculate_scores(pos)
    
    print(paste("Initial positions =", toString(pos)))
    print(paste("Initial score =", curr_score))
    
    # Initialization for iterations and data for plotting
    num_iter <- 0
    x_data <- numeric()
    y_data <- numeric()
    # iteration_max <- 30
    # scan_rate <- 7
    
    while (num_iter < iteration_max) {
        # Pick a position to improve
        pick_result <- pick_position(pos)
        ind <- pick_result$index
        imp <- pick_result$position
        l_bound <- max(pos[ind - 1], imp - scan_rate)
        r_bound <- min(pos[ind + 1], imp + scan_rate)
        
        # Optimize picked position
        new_imp <- optimize_position(pos, pos[ind], l_bound, r_bound, 1)
        pos[ind] <- new_imp
        
        # Print optimization details
        print(paste("Optimized target =", new_imp))
        print(paste("Modified pos =", toString(pos)))
        
        # Calculate change in score and update parameters
        new_score <- calculate_scores(pos)
        curr_score <- new_score
        num_iter <- num_iter + 1
        
        # Printing iteration details
        print(paste("#iteration =", num_iter, ", current score =", curr_score))
        
        # Data for plotting
        x_data <- c(x_data, num_iter)
        y_data <- c(y_data, curr_score)
    }
    
    plot(x_data, y_data, type = "b", col = "blue", xlab = "Iteration", ylab = "Score", main = "Score Optimization Over Iterations")
    
}
#execute_and_plot()