
target_gene <- "AATATGGGTATTAAAGGTTTGAATGCAATTATATCGGAACATGTTCCCTCTGCTATCAGGAAAAGCGATATCAAGAGCTTTTTTGGCAGAAAGGTTGCCATCGATGCCTCTATGTCTCTATATCAGTTTTTAATTGCTGTAAGACAGCAAGACGGTGGGCAGTTGACCAATGAAGCCGGTGAAACAACGTCACACTTGATGGGTATGTTTTATAGGACACTGAGAATGATTGATAACGGTATCAAGCCTTGTTATGTCTTCGACGGCAAACCTCCAGATTTGAAATCTCATGAGTTGACAAAGCGGTCTTCAAGAAGGGTGGAAACAGAAAAAAAACTGGCAGAGGCAACAACAGAATTGGAAAAGATGAAGCAAGAAAGAAGATTGGTGAAGGTTTCAAAAGAGCATAATGAAGAAGCCCAAAAATTACTAGGACTAATGGGAATCCCATATATAATAGCGCCAACGGAAGCTGAGGCTCAATGTGCTGAGTTGGCAAAGAAGGGAAAGGTGTATGCCGCAGCAAGTGAAGATATGGACACACTCTGTTATAGAACACCCTTCTTGTTGAGACATTTGACTTTTTCAGAGGCCAAGAAGGAACCGATTCACGAAATAGATACTGAATTAGTTTTGAGAGGACTCGACTTGACAATAGAGCAGTTTGTTGATCTTTGCATAATGCTTGGTTGTGACTACTGTGAAAGCATCAGAGGTGTTGGTCCAGTGACAGCCTTAAAATTGATAAAAACGCATGGATCCATCGAAAAAATCGTGGAGTTTATTGAATCTGGGGAGTCAAACAACACTAAATGGAAAATCCCAGAAGACTGGCCTTACAAACAAGCAAGAATGCTGTTTCTTGACCCTGAAGTTATAGATGGTAACGAAATAAACTTGAAATGGTCGCCACCAAAGGAGAAGGAACTTATCGAGTATTTATGTGATGATAAGAAATTCAGTGAAGAAAGAGTTAAATCTGGTATATCAAGATTGAAAAAAGGCTTGAAATCTGGCATTCAGGGTAGGTTAGATGGGTTCTTCCAAGTGGTGCCTAAGACAAAGGAACAGCTGGCTGCTGCGGCGAAAAGAGCACAAGAAAATAAAAAATTGAACAAAAATAAGAATAAAGTCACAAAGGGAAGAAGATGAGGG"
target_gene <- split_into_codons(target_gene)
print(target_gene)

num_codons <- length(target_gene) - 2
num_tiles <- 1:50

# Initialize cassette positions
pos <- seq(2, length(target_gene) - 3, by = length_tiles_global)
pos <- c(pos, length(target_gene) - 3)

# Obtain initial scores
curr_score <- calculate_scores(pos)
print(paste("Initial positions =", toString(pos)))
print(paste("Initial score =", curr_score))

# Customization
iteration_max <- 75
scan_rate <- 7

# Initialization for iterations and data for plotting
num_iter <- 0
x_data <- numeric()
y_data <- numeric()

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
