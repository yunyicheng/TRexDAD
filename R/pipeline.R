# --- SECTION: Preprocess Data -----------------


#' Split Gene Sequence into Codons
#'
#' This function takes a gene sequence (as a string) and splits it into codons. 
#' A codon consists of three nucleotides, and this function organizes the 
#' sequence into these triplets. It verifies that the input is a valid DNA sequence,
#' that its length is a multiple of 3, contains 'ATG' as the second codon,
#' and has an extra codon both before the start and after the end codon of the ORF.
#'
#' @param gene_sequence A character string representing the gene sequence.
#' 
#' @return A character vector where each element is a codon (a sequence of 
#'         three nucleotides), starting from the first nucleotide, with 'ATG' 
#'         as the second codon.
#'
#' @examples
#' # Example gene sequence
#' gene_seq <- "GCTATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATC"
#' split_into_codons(gene_seq)
#'
#' @export
split_into_codons <- function(gene_sequence) {
    # Verify that the sequence only contains valid nucleotides
    if(!grepl("^[ACGT]+$", gene_sequence)) {
        stop("Invalid gene sequence: sequence should only contain A, C, G, and T.")
    }
    
    # Check if the length of the sequence is a multiple of 3
    if (nchar(gene_sequence) %% 3 != 0) {
        stop("Invalid gene sequence length: length should be a multiple of 3.")
    }
    
    # Ensure 'ATG' is the second codon
    if (substr(gene_sequence, 4, 6) != "ATG") {
        stop("Invalid sequence: 'ATG' is not the second codon.")
    }
    
    # Split the sequence into codons starting from the first nucleotide
    codons <- strsplit(gene_sequence, "(?<=.{3})", perl=TRUE)[[1]]
    return(codons)
}


#' Calculate Oligo Cost
#'
#' This function calculates the cost of oligonucleotides (oligos) based on 
#' the number of tiles and the number of codons. It considers the cost of 
#' ccdB (a cloning vector) and the cost of base replacement per tile.
#' The formula for calculation is obtained from ANNB lab in which I did my
#' research project.
#'
#' @param num_tiles The number of tiles in the oligo sequence.
#' @param num_codons The number of codons in the oligo sequence.
#'
#' @return The calculated price per codon, considering both ccdB cost and 
#'         base replacement cost.
#' 
#' @examples
#' # Calculate oligo cost for 10 tiles and 100 codons
#' oligo_cost(10, 100)
#' 
#' @export
oligo_cost <- function(num_tiles, num_codons) {
    ccdB_cost <- num_tiles * 16
    one_br_per_tile <- (num_codons * 3 / num_tiles + 30) * 0.01
    complete_replacement <- one_br_per_tile * num_codons
    total <- complete_replacement + ccdB_cost
    price_per_codon <- total / num_codons
    return(price_per_codon)
}


#' Calculate Optimal Number of Tiles and Cost
#'
#' This function determines the optimal number of tiles for a given number of
#' codons. It calculates the oligo cost for a range of tile numbers and finds 
#' the number of tiles that minimizes this cost. In addition, the function computes 
#' the global length of the tiles and the optimal cost associated with the optimal 
#' number of tiles. The function limits the number of tiles to 50 at maximum, based 
#' on the reference paper's indication that assembly fidelity decreases with an 
#' increasing number of overhang pairs.
#'
#' @param num_codons An integer representing the total number of codons.
#'
#' @return A list containing three elements: `optimal_tiles`, the number of tiles 
#'         that minimize the oligo cost, `length_tiles_global`, the average 
#'         length of these tiles, and `optimal_cost`, the cost associated with 
#'         the optimal number of tiles. Additionally, the function prints a 
#'         structured sentence summarizing these outputs.
#'
#' @references 
#' Pryor, J. M. et al. Enabling one-pot Golden Gate assemblies of 
#' unprecedented complexity using data-optimized assembly design. PLOS ONE 15, 
#' e0238592 (2020).
#'
#' @examples
#' # Calculate optimal tiles and cost for 100 codons
#' calculate_optimal_tiles(100)
#'
#' @export
calculate_optimal_tiles <- function(num_codons) {
    num_tiles_range <- 1:50
    costs <- sapply(num_tiles_range, function(x) oligo_cost(x, num_codons))
    optimal_tiles <- which.min(costs)
    optimal_cost <- min(costs)
    length_tiles <- round(num_codons / optimal_tiles)
    
    # Print structured sentence
    cat("For", num_codons, "codons, the optimal number of tiles is", optimal_tiles, 
        "with an average length of", length_tiles, "nucleotides per tile and an optimal cost of", 
        optimal_cost, "per codon.\n")
    
    return(list(optimal_tiles = optimal_tiles, length_tiles = length_tiles, optimal_cost = optimal_cost))
}


# --- SECTION: Calculate Scores -----------------


#' Get Overhang Sequences
#'
#' This function calculates the overhang sequences at the specified positions 
#' in a gene sequence. Depending on the flag, it returns these sequences as 
#' either `DNAString` objects or plain character strings.
#'
#' @param start An integer indicating the start position of the tile.
#' @param end An integer indicating the end position of the tile.
#' @param flag A binary flag (0 or 1) indicating the return type of the sequences. 
#'        If 1, sequences are returned as `DNAString` objects; if 0, as character strings.
#'
#' @return A list containing two elements: `head` and `tail`, representing the 
#'         overhang sequences at the specified positions.
#'
#' @import Biostrings
#' @import forstringr
#' @export
#' @examples
#' # Assuming target_gene is a predefined character vector of codons
#' target_gene <- c("GAG", "CTG", "TGT", "AGG", "TGC", "CGG", "CCA", 
#' "ATT", "TGA" "TAG" "GAA" "TAT" "AGC")
#' get_overhangs(2, 4, 1)
get_overhangs <- function(start, end, flag) {
    # Extract the overhangs
    head <- paste(forstringr::str_right(target_gene[start - 2], 1), 
                  forstringr::str_left(target_gene[start - 1], 3), sep = "")
    tail <- paste(forstringr::str_right(target_gene[end + 1], 3), 
                  forstringr::str_left(target_gene[end + 2], 1), sep = "")
    
    # If flag is true, return as DNAString objects; otherwise, as character strings
    if (flag == 1) {
        head <- Biostrings::DNAString(head)
        tail <- Biostrings::DNAString(tail)
    } else {
        head <- as.character(head)
        tail <- as.character(tail)
    }
    
    return(list(head = head, tail = tail))
}

#' Get All Overhang Sequences for a List of Positions
#'
#' This function computes overhang sequences for a given list of positions 
#' in a gene sequence. It iteratively calls the `get_overhangs` function to 
#' calculate the head and tail overhangs for each position and accumulates 
#' them in a list.
#'
#' @param pos_lst A numeric vector representing positions in the gene sequence.
#'
#' @return A list of overhang sequences corresponding to each position in 
#'         `pos_lst`. Each position contributes two elements to the list: the 
#'         head and tail overhangs.
#'
#' @importFrom Biostrings DNAString
#' @export
get_all_overhangs <- function(pos_lst) {
    all_overhangs <- list()
    
    for (i in 1:(length(pos_lst) - 1)) {
        overhangs <- get_overhangs(pos_lst[i], pos_lst[i + 1] - 1, 1)
        all_overhangs <- c(all_overhangs, overhangs$head, overhangs$tail)
    }
    
    return(all_overhangs)
}

#' Obtain Score for a Given Tile
#'
#' This function calculates a score for a specific tile within a gene sequence.
#' The score is based on palindromicity, length variation from a global standard,
#' and on-target reactivity, using a pre-defined overhang fidelity dataframe.
#'
#' @param start An integer representing the start position of the tile in the gene sequence.
#' @param end An integer representing the end position of the tile in the gene sequence.
#'
#' @return The calculated score for the specified tile, which takes into account
#'         the palindromicity, length variation, and on-target reactivity.
#'
#' @importFrom Biostrings reverseComplement DNAString
#' @examples
#' # Assuming predefined variables and setup:
#' # start = 1, end = 3, overhang_fidelity (dataframe), length_tiles_global (value)
#' obtain_score(1, 3)
#'
#' @export
obtain_score <- function(start, end) {
    overhang_fidelity_df <- as.data.frame(overhang_fidelity)
    # Set the first column as row names
    rownames(overhang_fidelity_df) <- overhang_fidelity_df[[1]]
    overhang_fidelity_df <- overhang_fidelity_df[-1]
    
    overhangs <- get_overhangs(start, end, 1)
    head <- overhangs$head
    tail <- overhangs$tail
    
    palindromicity <- 0
    if (Biostrings::reverseComplement(head) == head) {
        palindromicity <- palindromicity + 1
    } else {}
    if (Biostrings::reverseComplement(tail) == tail) {
        palindromicity <- palindromicity + 1
    } else {}
    delta_len <- end - start + 1 - length_tiles_global
    
    head_comp <- Biostrings::reverseComplement(head)
    tail_comp <- Biostrings::reverseComplement(tail)
    
    head <- as.character(head)
    head_comp <- as.character(head_comp)
    tail <- as.character(tail)
    tail_comp <- as.character(tail_comp)
    head_reactivity <- overhang_fidelity_df[head, head_comp]
    tail_reactivity<- overhang_fidelity_df[tail, tail_comp]
    
    on_target <- head_reactivity + tail_reactivity
    local_score <- on_target - palindromicity * 1000 - 5 * delta_len^2
    return(local_score)
}

#' Calculate Global Score for a List of Positions
#'
#' This function calculates a global score for a list of positions in a gene sequence. 
#' The score takes into account off-target reactions, repetitions, and the sum of 
#' local scores. It utilizes a pre-defined overhang fidelity dataframe for calculations.
#'
#' @param pos_lst A numeric vector representing positions in the gene sequence.
#'
#' @return The global score, which is a composite measure considering the off-target 
#'         reactions, repetitions in overhangs, and the sum of local scores for each position.
#'
#' @importFrom Biostrings reverseComplement DNAString
#' @export
calculate_scores <- function(pos_lst) {
    overhang_fidelity_df <- as.data.frame(overhang_fidelity)
    # Set the first column as row names
    rownames(overhang_fidelity_df) <- overhang_fidelity_df[[1]]
    overhang_fidelity_df <- overhang_fidelity_df[-1]
    
    overhangs <- get_all_overhangs(pos_lst)
    overhangs_chr <- lapply(overhangs, as.character)
    complements <- lapply(overhangs, Biostrings::reverseComplement)
    complements_chr <- lapply(complements, as.character)
    
    off_reaction <- 0
    overhang_comb <- combn(overhangs_chr, 2, simplify = FALSE)
    complement_comb <- combn(complements_chr, 2, simplify = FALSE)
    
    for (a in (1:length(overhang_comb))) {
        oc <- overhang_comb[[a]]
        # off_reaction <- off_reaction + (overhang_fidelity %>% 
        #     filter(Overhang == oc[[1]]) %>% pull(oc[[2]]))
        off_reaction <- overhang_fidelity_df[oc[[1]], oc[[2]]]
    }
    for (b in (1:length(complement_comb))) {
        cc <- complement_comb[[b]]
        # off_reaction <- off_reaction + (overhang_fidelity %>% 
        #     filter(Overhang == cc[[1]] %>% pull(cc[[2]])))
        off_reaction <- overhang_fidelity_df[cc[[1]], cc[[2]]]
    }
    
    for (i in seq_along(overhangs_chr)) {
        for (j in seq_along(complements_chr)) {
            if (i != j) {
                off_reaction <- off_reaction + 
                    (overhang_fidelity_df[overhangs_chr[[i]], overhangs_chr[[j]]])
            } 
        }
    }
    
    repetition <- length(overhangs) - length(unique(overhangs))
    local_score <- sum(sapply(seq_along(pos_lst)[-length(pos_lst)], function(i) obtain_score(pos_lst[i], pos_lst[i + 1] - 1)))
    global_score <- 50 * local_score - 100 * off_reaction - 1000 * repetition
    
    return(global_score)
}


# --- SECTION: Optimize Positions -----------------


#' Pick a Position for Optimization
#'
#' This function selects a position from a list of positions for optimization 
#' based on their scores. Positions with lower scores are more likely to be picked.
#' It uses a weighted sampling approach where the weights are inversely proportional 
#' to the scores of the positions.
#'
#' @param pos_lst A numeric vector representing positions in a gene sequence.
#'
#' @return A list with two elements: `index`, the index in `pos_lst` of the picked 
#'         position (shifted by 1), and `position`, the actual value of the picked position.
#'
#' @examples
#' # Example usage assuming pos_lst and obtain_score function are defined
#' pos_lst <- c(1, 5, 10, 15)
#' picked <- pick_position(pos_lst)
#' print(picked)
#'
#' @export
pick_position <- function(pos_lst) {
    score_weights <- numeric(length(pos_lst) - 2)
    indexes <- numeric(length(pos_lst) - 2)
    
    for (i in 1:(length(pos_lst) - 2)) {
        new_score <- obtain_score(pos_lst[i], pos_lst[i + 1] - 1)
        score_weights[i] <- -1 * new_score
        indexes[i] <- i
    }
    
    offset <- min(score_weights)
    score_weights <- score_weights - offset + 1
    pick_index <- sample(indexes, size = 1, prob = score_weights)
    position <- pos_lst[pick_index + 1]
    
    return(list(index = pick_index + 1, position = position))
}

#' Optimize a Single Position in a Position List
#'
#' This function optimizes a single position within a list of positions. It adjusts
#' the specified position to maximize the overall score (obtained via `calculate_scores`).
#' The function supports different optimization modes, including greedy and Markov Chain
#' Monte Carlo (MCMC) approaches.
#'
#' @param pos_lst A numeric vector representing the current list of positions.
#' @param curr_pos The current position value that needs optimization.
#' @param left The left boundary of the scanning range for position optimization.
#' @param right The right boundary of the scanning range for position optimization.
#' @param mode The optimization mode: 1 for greedy, 2 for MCMC.
#'
#' @return The optimized position value within the specified range. If the function
#'         does not find a better position, it returns the original position value.
#'
#' @examples
#' # Example usage assuming pos_lst and calculate_scores function are defined
#' pos_lst <- c(1, 5, 10, 15)
#' optimized_pos <- optimize_position(pos_lst, pos_lst[2], 3, 7, 1)
#' print(optimized_pos)
#'
#' @export
optimize_position <- function(pos_lst, curr_pos, left, right, mode) {
    score_weights <- numeric()
    indexes <- numeric()
    copy_pos_lst <- pos_lst
    original_score <- calculate_scores(copy_pos_lst)
    
    for (n in left:right) {
        curr_pos_ind <- match(curr_pos, copy_pos_lst)
        copy_pos_lst[curr_pos_ind] <- n
        s <- calculate_scores(copy_pos_lst)
        
        if (s >= original_score) {
            indexes <- c(indexes, n)
            score_weights <- c(score_weights, s)
        }
        curr_pos <- n
    }
    
    offset <- min(score_weights)
    score_weights <- score_weights - offset + 1
    
    # GREEDY
    if (mode == 1) {
        result_index <- indexes[which.max(score_weights)]
        # MCMC
    } else if (mode == 2) {
        result_index <- sample(indexes, size = 1, prob = score_weights)
    }
    
    return(result_index)
}


# --- SECTION: Execution and Plot -----------------


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
#' @references 
#' Vladimir Potapov, Jennifer L. Ong, Rebecca B. Kucera, Bradley W. Langhorst, Katharina Bilotti, John M. Pryor, Eric J. Cantor, Barry Canton, Thomas F. Knight, Thomas C. Jr. Evans, and Gregory J. S. Lohman. 
#'      Comprehensive Profiling of Four Base Overhang Ligation Fidelity by T4 DNA Ligase and Application to DNA Assembly. 
#'      ACS Synthetic Biology, 7(11):2665–2674, 2018.
#' Pryor JM, Potapov V, Kucera RB, Bilotti K, Cantor EJ, et al. 
#'      Enabling one-pot Golden Gate assemblies of unprecedented complexity using data-optimized assembly design. 
#'      PLOS ONE 15(9): e0238592, 2020.
#' @import readxl
#' 
#' @export
execute_and_plot <- function(iteration_max=30, scan_rate=7) {
    
    # Read overhang fidelity chart
    overhang_fidelity <<- read.csv("inst/extdata/overhang_fidelity.csv")
    
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
execute_and_plot()