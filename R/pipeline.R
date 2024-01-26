# Purpose: Pipeline Implementation
# Author: Yunyi Cheng
# Date: 01.07.2024
# Version: 1.0.0
# Bugs and Issues: N/A


# --- SECTION: Setup and Global Variable Initialization -----------------


# Read overhang fidelity chart
OVERHANG_FIDELITY <- utils::read.csv("inst/extdata/overhang_fidelity.csv")

# Rad-27 gene for example usage
RAD_27 <- "AATATGGGTATTAAAGGTTTGAATGCAATTATATCGGAACATGTTCCCTCTGCTATCAGGAAAAGCGATATCAAGAGCTTTTTTGGCAGAAAGGTTGCCATCGATGCCTCTATGTCTCTATATCAGTTTTTAATTGCTGTAAGACAGCAAGACGGTGGGCAGTTGACCAATGAAGCCGGTGAAACAACGTCACACTTGATGGGTATGTTTTATAGGACACTGAGAATGATTGATAACGGTATCAAGCCTTGTTATGTCTTCGACGGCAAACCTCCAGATTTGAAATCTCATGAGTTGACAAAGCGGTCTTCAAGAAGGGTGGAAACAGAAAAAAAACTGGCAGAGGCAACAACAGAATTGGAAAAGATGAAGCAAGAAAGAAGATTGGTGAAGGTTTCAAAAGAGCATAATGAAGAAGCCCAAAAATTACTAGGACTAATGGGAATCCCATATATAATAGCGCCAACGGAAGCTGAGGCTCAATGTGCTGAGTTGGCAAAGAAGGGAAAGGTGTATGCCGCAGCAAGTGAAGATATGGACACACTCTGTTATAGAACACCCTTCTTGTTGAGACATTTGACTTTTTCAGAGGCCAAGAAGGAACCGATTCACGAAATAGATACTGAATTAGTTTTGAGAGGACTCGACTTGACAATAGAGCAGTTTGTTGATCTTTGCATAATGCTTGGTTGTGACTACTGTGAAAGCATCAGAGGTGTTGGTCCAGTGACAGCCTTAAAATTGATAAAAACGCATGGATCCATCGAAAAAATCGTGGAGTTTATTGAATCTGGGGAGTCAAACAACACTAAATGGAAAATCCCAGAAGACTGGCCTTACAAACAAGCAAGAATGCTGTTTCTTGACCCTGAAGTTATAGATGGTAACGAAATAAACTTGAAATGGTCGCCACCAAAGGAGAAGGAACTTATCGAGTATTTATGTGATGATAAGAAATTCAGTGAAGAAAGAGTTAAATCTGGTATATCAAGATTGAAAAAAGGCTTGAAATCTGGCATTCAGGGTAGGTTAGATGGGTTCTTCCAAGTGGTGCCTAAGACAAAGGAACAGCTGGCTGCTGCGGCGAAAAGAGCACAAGAAAATAAAAAATTGAACAAAAATAAGAATAAAGTCACAAAGGGAAGAAGATGAGGG"


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
    # Check for negative input values
    if (num_tiles <= 0) {
        stop("The number of tiles should be a positive integer.")
    } else {}
    if (num_codons <= 0) {
        stop("The number of codons should be a positive integer.")
    } else {}
    
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
#' the optimal length of the tiles and the optimal cost associated with the optimal 
#' number of tiles. The function limits the number of tiles to 50 at maximum, based 
#' on the reference paper's indication that assembly fidelity decreases with an 
#' increasing number of overhang pairs.
#'
#' @param num_codons An integer representing the total number of codons.
#'
#' @return A list containing three elements: `optimal_tiles`, the number of tiles 
#'         that minimize the oligo cost, `optimal_tile_length`, the average 
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
    # Validate num_codons input
    if (num_codons <= 0) {
        stop("The number of codons should be a positive integer.")
    } else {}
    
    num_tiles_range <- 1:50
    costs <- sapply(num_tiles_range, function(x) oligo_cost(x, num_codons))
    optimal_num_tiles <- which.min(costs)
    optimal_cost <- min(costs)
    optimal_tile_length <- round(num_codons / optimal_num_tiles)
    
    # Print structured sentence
    cat("For", num_codons, "codons, the optimal number of tiles is", 
        optimal_num_tiles, "with an average length of", optimal_tile_length, 
        "nucleotides per tile and an optimal cost of", 
        optimal_cost, "per codon.\n")
    
    return(list(optimal_num_tiles = optimal_num_tiles, 
                optimal_tile_length = optimal_tile_length, 
                optimal_cost = optimal_cost))
}


# --- SECTION: Calculate Scores -----------------


#' Get Overhang Sequences for Single Tile
#'
#' This function calculates the overhang sequences at the specified positions 
#' in a gene sequence. Depending on the flag, it returns these sequences as 
#' either `DNAString` objects or plain character strings.
#'
#' @param gene_codons A list of splitted codons of target genes
#' @param start An integer indicating the start position of the tile.
#' @param end An integer indicating the end position of the tile.
#' @param as_dna_strings A boolean indicating the return type of the sequences. 
#'        If TRUE, sequences are returned as `DNAString` objects; 
#'        if FALSE, as character strings.
#'
#' @return A list containing two elements: `head` and `tail`, representing the 
#'         overhang sequences at the specified positions.
#'
#' @import forstringr
#' @importFrom Biostrings DNAString
#' @export
#' 
#' @references 
#' Pryor, J. M. et al. Enabling one-pot Golden Gate assemblies of 
#' unprecedented complexity using data-optimized assembly design. PLOS ONE 15, 
#' e0238592 (2020).
#' 
#' @examples
#' # Overhangs for sample gene_codons at position(start, end)
#' get_overhangs(c("GAG", "CTG", "TGT", "AGG", "TGC", "CGG", "CCA", 
#' "ATT", "TGA", "TAG", "GAA", "TAT", "AGC"), 3, 4, TRUE)
get_overhangs <- function(gene_codons, start, end, as_dna_strings) {
    # Check if start and end positions are valid
    if ((start <= 2) || (end >= length(gene_codons) - 1) || (end < start)) {
        stop("Invalid start or end positions.")
    } else {}
    
    # Check for invalid input codons
    for(codon in gene_codons) {
        if (!grepl("^[ACGT]{3}$", codon)) {
            stop("Invalid codon detected.")
        }
    }
    
    # Extract the overhangs
    head <- paste(forstringr::str_right(gene_codons[start - 2], 1), 
                  forstringr::str_left(gene_codons[start - 1], 3), sep = "")
    tail <- paste(forstringr::str_right(gene_codons[end + 1], 3), 
                  forstringr::str_left(gene_codons[end + 2], 1), sep = "")
    
    # If flag is true, return as DNAString objects; otherwise, as character strings
    if (as_dna_strings) {
        head <- Biostrings::DNAString(head)
        tail <- Biostrings::DNAString(tail)
    } else {
        head <- as.character(head)
        tail <- as.character(tail)
    }
    # print(list(head = head, tail = tail))
    return(list(head = head, tail = tail))
}


#' Get All Overhang Sequences for a List of Positions
#'
#' This function computes overhang sequences for a given list of positions 
#' in a gene sequence. It iteratively calls the `get_overhangs` function to 
#' calculate the head and tail overhangs for each position and accumulates 
#' them in a list.
#'
#' @param gene_codons A list of splitted codons of target genes
#' @param pos_lst A numeric vector representing positions in the gene sequence. Each element in pos_lst is the start of a gene tile.
#'
#' @return A list of overhang sequences corresponding to each position in 
#'         `pos_lst`. Each position contributes two elements to the list: the 
#'         head and tail overhangs.
#'
#' @export
#' 
#' @references 
#' Pryor, J. M. et al. Enabling one-pot Golden Gate assemblies of 
#' unprecedented complexity using data-optimized assembly design. PLOS ONE 15, 
#' e0238592 (2020).
#' 
#' @examples
#' gene_codons <- c("GAG", "CTG", "TGT", "AGG", "TGC", "CGG", "CCA", 
#'                     "ATT", "TGA", "TAG", "GAA", "TGA", "AGC")
#' pos_lst <- c(3, 5, length(gene_codons) - 2)
#' get_all_overhangs(gene_codons, pos_lst)
get_all_overhangs <- function(gene_codons, pos_lst) {
    if (any(pos_lst <= 0)) {
        stop("pos_lst must contain only positive integers.")
    } else {}
    
    if (!(length(pos_lst >= 3))) {
        stop("pos_lst must contain at least three positions.")
    } else {}
    
    # Check if the last position in pos_lst is equal to the length of gene_codons - 2
    if (tail(pos_lst, 1) != length(gene_codons) - 2) {
        stop("The last position in pos_lst must equal length(gene_codons) - 2.")
    } else {}
    
    # Check if pos_lst is in increasing order
    if (!(all(diff(pos_lst) > 0))){
        stop("pos_lst must be in increasing order.")
    } else {}
    
    all_overhangs <- list()
    for (i in 1:(length(pos_lst) - 2)) {
        overhangs <- get_overhangs(gene_codons, pos_lst[i], pos_lst[i + 1] - 1, TRUE)
        all_overhangs <- c(all_overhangs, overhangs$head, overhangs$tail)
    }
    
    last_start <- pos_lst[length(pos_lst) - 1]
    last_end <- pos_lst[length(pos_lst)]
    last_pair <- get_overhangs(gene_codons, pos_lst[length(pos_lst) - 1], last_end, 1)
    all_overhangs <- c(all_overhangs, last_pair$head, last_pair$tail)
    return(all_overhangs)
}

gene_codons <- c("GAG", "CTG", "TGT", "AGG", "TGC", "CGG", "CCA", 
                 "ATT", "TGA", "TAG", "GAA", "TGA", "AGC")
pos_lst <- c(3, 5, length(gene_codons) - 2)
length(pos_lst)
get_all_overhangs(gene_codons, pos_lst)

#' Obtain Score for a Given Tile
#'
#' This function calculates a score for a specific tile within a gene sequence.
#' The score is based on palindromicity, length variation from a global standard,
#' and on-target reactivity, using a pre-defined overhang fidelity dataframe.
#'
#' @param gene_codons A vector of DNA codons.
#' @param tile_length An integer representing the expected length of the tile.
#' @param start An integer representing the start position of the tile in the gene sequence.
#' @param end An integer representing the end position of the tile in the gene sequence.
#'
#' @return An integer representing the calculated score for the specified tile.
#'         The score takes into account the palindromicity, length variation,
#'         and on-target reactivity of the tile's overhang sequences.
#'
#' @importFrom Biostrings reverseComplement
#' @export
#' 
#' @references
#' Pryor, J. M. et al. Enabling one-pot Golden Gate assemblies of 
#' unprecedented complexity using data-optimized assembly design. PLOS ONE 15, 
#' e0238592 (2020).
#' 
#' @examples
#' # Example usage with predefined variables and setup:
#' gene_codons <- c("ATG", "CAG", "TAC", "GGA", "TGA", "TCA")
#' tile_length <- 4
#' start <- 3
#' end <- 4
#' calculate_local_score(gene_codons, tile_length, start, end)
calculate_local_score <- function(gene_codons, tile_length, start, end) {
    # Validate inputs
    if (tile_length <= 0) {
        stop("Tile length must be a positive integer.")
    }
    
    # Load from global variable
    overhang_fidelity_df <- as.data.frame(OVERHANG_FIDELITY)
    # Set the first column as row names
    rownames(overhang_fidelity_df) <- overhang_fidelity_df[[1]]
    overhang_fidelity_df <- overhang_fidelity_df[-1]
    
    overhangs <- get_overhangs(gene_codons,start, end, TRUE)
    head <- overhangs$head
    tail <- overhangs$tail
    
    palindromicity <- 0
    if (Biostrings::reverseComplement(head) == head) {
        palindromicity <- palindromicity + 1
    } else {}
    if (Biostrings::reverseComplement(tail) == tail) {
        palindromicity <- palindromicity + 1
    } else {}
    delta_len <- end - start + 1 - tile_length
    
    head_comp <- Biostrings::reverseComplement(head)
    tail_comp <- Biostrings::reverseComplement(tail)
    
    head <- as.character(head)
    head_comp <- as.character(head_comp)
    tail <- as.character(tail)
    tail_comp <- as.character(tail_comp)
    head_reactivity <- overhang_fidelity_df[head, head_comp]
    tail_reactivity<- overhang_fidelity_df[tail, tail_comp]
    
    on_target <- head_reactivity + tail_reactivity
    local_score <- (on_target - (palindromicity * 1000) - (5 * delta_len^2))
    return(local_score)
}

#' Calculate Global Score for a List of Positions
#'
#' This function calculates a global score for a list of positions in a gene sequence. 
#' The score takes into account off-target reactions, repetitions, and the sum of 
#' local scores. It utilizes a pre-defined overhang fidelity dataframe for calculations.
#'
#' @param gene_codons A vector of DNA codons.
#' @param tile_length An integer representing the expected length of the tile.
#' @param pos_lst A numeric vector representing positions in the gene sequence.
#'
#' @return The global score, which is a composite measure considering the off-target 
#'         reactions, repetitions in overhangs, and the sum of local scores for each position.
#'
#' @importFrom Biostrings reverseComplement
#' @import utils
#' @export
#' 
#' @references
#' Pryor, J. M. et al. Enabling one-pot Golden Gate assemblies of 
#' unprecedented complexity using data-optimized assembly design. PLOS ONE 15, 
#' e0238592 (2020).
#' 
#' @examples
#' gene_codons <- c("GAG", "ATG", "TGT", "AGG", "TGC", "CGG", "CCA", 
#'                     "ATT", "TGA", "TAG", "GAA", "TGA", "AGC")
#' tile_length <- 5
#' pos_lst <- c(3, 6, 8, 11)
#' calculate_global_score(gene_codons, tile_length, pos_lst)
calculate_global_score <- function(gene_codons, tile_length, pos_lst) {
    overhang_fidelity_df <- as.data.frame(OVERHANG_FIDELITY)
    # Set the first column as row names
    rownames(overhang_fidelity_df) <- overhang_fidelity_df[[1]]
    overhang_fidelity_df <- overhang_fidelity_df[-1]
    
    overhangs <- get_all_overhangs(gene_codons, pos_lst)
    overhangs_chr <- lapply(overhangs, as.character)
    # print(overhangs_chr)
    complements <- lapply(overhangs, Biostrings::reverseComplement)
    complements_chr <- lapply(complements, as.character)
    
    off_reaction <- 0
    overhang_comb <- utils::combn(overhangs_chr, 2, simplify = FALSE)
    complement_comb <- utils::combn(complements_chr, 2, simplify = FALSE)
    
    for (a in seq_along(overhang_comb)) {
        oc <- overhang_comb[[a]]
        # off_reaction <- off_reaction + (overhang_fidelity %>% 
        #     filter(Overhang == oc[[1]]) %>% pull(oc[[2]]))
        off_reaction <- overhang_fidelity_df[oc[[1]], oc[[2]]]
    }
    for (b in seq_along(complement_comb)) {
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
    local_score <- sum(sapply(seq_along(pos_lst)[-length(pos_lst)], 
                              function(i) calculate_local_score(
                                  gene_codons, tile_length, 
                                  pos_lst[i], pos_lst[i + 1] - 1)))
    global_score <- ((50 * local_score) 
                     - (100 * off_reaction) - (1000 * repetition))
    
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
#' @param gene_codons A vector of DNA codons.
#' @param tile_length An integer representing the expected length of the tile.
#' @param pos_lst A numeric vector representing positions in a gene sequence.
#'
#' @return A list with two elements: `index`, the index in `pos_lst` of the picked 
#'         position (shifted by 1), and `position`, the actual value of the picked position.
#'
#' @export
#' 
#' @examples
#' gene_codons <- c("GAG", "ATG", "TGT", "AGG", "TGC", "CGG", "CCA", 
#'                   "ATT", "TGA", "TAG", "GAA", "TGA", "AGC")
#' tile_length <- 5
#' pos_lst <- c(3, 6, 8, 11)
#' pick_position(gene_codons, tile_length, pos_lst)
pick_position <- function(gene_codons, tile_length, pos_lst) {
    score_weights <- numeric(length(pos_lst) - 2)
    indexes <- numeric(length(pos_lst) - 2)
    
    for (i in 1:(length(pos_lst) - 2)) {
        new_score <- calculate_local_score(gene_codons, tile_length, 
                                           pos_lst[i], pos_lst[i + 1] - 1)
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
#' This function optimizes a single position within a list of positions. 
#' It adjusts the specified position to maximize the overall score 
#' (obtained via `calculate_global_score`).
#' The function supports different optimization modes, 
#' including greedy and Markov Chain Monte Carlo (MCMC) approaches.
#'
#' @param gene_codons A vector of DNA codons.
#' @param tile_length An integer representing the expected length of the tile.
#' @param pos_lst A numeric vector representing the current list of positions.
#' @param curr_pos The current position value that needs optimization.
#' @param left The left boundary of the scanning range for position optimization.
#' @param right The right boundary of the scanning range for position optimization.
#' @param use_greedy The optimization mode: TRUE for greedy, FALSE for MCMC.
#'
#' @return The optimized position value within the specified range. 
#'         If the function does not find a better position, 
#'         it returns the original position value.
#' @export 
#'      
#' @references Hastings, W. K. Monte Carlo sampling methods using Markov chains
#'             and their applications (1970). Biometrika 57, 97–109.
#'             
#' @examples
#' gene_codons <- c("GAG", "ATG", "TGT", "AGG", "TGC", "CGG", "CCA", 
#' "ATT", "TGA", "TAG", "GAA", "TGA", "AGC")
#' tile_length <- 5
#' pos_lst <- c(3, 6, 11)
#' curr_pos <- pos_lst[2]
#' left <- 7
#' right <- 10
#' optimize_position(gene_codons, tile_length, pos_lst, curr_pos, 
#'                     left, right, TRUE)
optimize_position <- function(gene_codons, tile_length, pos_lst, 
                              curr_pos, left, right, use_greedy) {
    # Validate inputs
    if (!is.vector(gene_codons) || any(!grepl("^[ACGT]+$", gene_codons))) {
        stop("gene_codons must be a vector of valid DNA codons.")
    }
    if (!is.numeric(tile_length) || tile_length <= 0) {
        stop("tile_length must be a positive integer.")
    }
    if (!is.numeric(pos_lst) || any(pos_lst <= 0)) {
        stop("pos_lst must contain only positive integers.")
    }
    if (!curr_pos %in% pos_lst) {
        stop("curr_pos must be an element in pos_lst.")
    }
    if (!is.numeric(left) || !is.numeric(right) || left > right || left <= 0 || right > length(gene_codons)) {
        stop("left and right must be valid bounds within the gene_codons range.")
    }
    
    # Early return if left and right are equal
    if (left == right) {
        return(curr_pos)
    } else {}
    
    score_weights <- numeric()
    indexes <- numeric()
    original_score <- calculate_global_score(gene_codons, tile_length, pos_lst)
    
    for (n in left:right) {
        curr_pos_ind <- match(curr_pos, pos_lst)
        pos_lst[curr_pos_ind] <- n
        s <- calculate_global_score(gene_codons, tile_length, pos_lst)
        
        if (s >= original_score) {
            indexes <- c(indexes, n)
            score_weights <- c(score_weights, s)
        }
        curr_pos <- n
    }
    
    offset <- min(score_weights)
    score_weights <- score_weights - offset + 1
    
    # GREEDY
    if (use_greedy) {
        result_index <- indexes[which.max(score_weights)]
    # MCMC
    } else {
        result_index <- sample(indexes, size = 1, prob = score_weights)
    }
    
    return(result_index)
}


# --- SECTION: Execute and Plot -----------------


#' Execute the Optimization Process and Plot Results
#'
#' This function performs an optimization process to determine the optimal tile
#' positions in a gene sequence and plots the progression of the score 
#' over iterations. It reads the overhang fidelity data, computes the optimal 
#' number of tiles, and iteratively improves the positions, plotting the change 
#' in score with respect to number of iterations at the end.
#' 
#' @param target_gene The target gene for assembly (default: Rad27), 
#'                    must contain one codon before and after ORF. 
#' @param max_iter The maximum number of iterations for the optimization process
#'                 (default: 30).
#' @param scan_rate The range within which to scan for optimizing tile positions 
#'                 (default: 7).
#'
#' @return This function does not return a value; it generates a plot showing 
#' the score optimization over iterations.
#' 
#' @import ggplot2 scales
#' @export
#' 
#' @references 
#' Vladimir Potapov, Jennifer L. Ong, Rebecca B. Kucera, Bradley W. Langhorst, 
#'      Katharina Bilotti, John M. Pryor, Eric J. Cantor, Barry Canton, 
#'      Thomas F. Knight, Thomas C. Jr. Evans, and Gregory J. S. Lohman. 
#'      Comprehensive Profiling of Four Base Overhang Ligation Fidelity by 
#'      T4 DNA Ligase and Application to DNA Assembly. 
#'      ACS Synthetic Biology, 7(11):2665–2674, 2018.
#' Pryor JM, Potapov V, Kucera RB, Bilotti K, Cantor EJ, et al. 
#'      Enabling one-pot Golden Gate assemblies of unprecedented 
#'      complexity using data-optimized assembly design. 
#'      PLOS ONE 15(9): e0238592, 2020.
#' Alberts B., Johnson A., Lewis J., et al. 
#'      "Molecular Biology of the Cell." 6th edition. 
#'      Garland Science, 2014.
#' 
#' @examples
#' execute_and_plot(max_iter = 5, scan_rate = 2)
execute_and_plot <- function(target_gene = RAD_27, 
                             max_iter = 30, scan_rate = 7) {
    
    # Split target gene into codons
    gene_codons <- split_into_codons(target_gene)
    
    # Check if the number of codons is 60 or less (from citation 3 of function)
    # an ORF would typically need to be at least 60-90 nucleotides long to 
    # encode a functional protein of minimal size.
    if (length(gene_codons) < 60) {
        stop("Target gene produces 60 or fewer codons, which is insufficient 
             for the optimization process.")
    }
    
    # Calculate optimal number of tiles
    num_codons <- length(gene_codons)
    tile_result <- (calculate_optimal_tiles(num_codons))
    num_tile <- tile_result$optimal_num_tiles
    tile_length <- tile_result$optimal_tile_length
    
    # Initialize tile positions
    pos <- seq(3, length(gene_codons) - 2, by = tile_length)
    pos <- c(pos, length(gene_codons) - 2)
    print(paste("Initial positions =", toString(pos)))
    # Obtain initial scores
    curr_score <- calculate_global_score(gene_codons, tile_length, pos)
    print(paste("Initial score =", curr_score))
    
    # Initialization for iterations and data for plotting
    num_iter <- 0
    x_data <- numeric(max_iter)
    y_data <- numeric(max_iter)
    
    while (num_iter < max_iter) {
        # Pick a position to improve
        pick_result <- pick_position(gene_codons, tile_length, pos)
        ind <- pick_result$index
        imp <- pick_result$position
        l_bound <- max(pos[ind - 1], imp - scan_rate)
        r_bound <- min(pos[ind + 1], imp + scan_rate)
        
        # Optimize picked position
        new_imp <- optimize_position(
            gene_codons, tile_length, pos, pos[ind], l_bound, r_bound, TRUE)
        pos[ind] <- new_imp
        
        # Print optimization details
        print(paste("Optimized target =", new_imp))
        print(paste("Modified pos =", toString(pos)))
        
        # Calculate change in score and update parameters
        new_score <- calculate_global_score(gene_codons, tile_length, pos)
        curr_score <- new_score
        num_iter <- num_iter + 1
        
        # Printing iteration details
        print(paste("#iteration =", num_iter, ", current score =", curr_score))
        
        # Data for plotting
        x_data <- c(x_data, num_iter)
        y_data <- c(y_data, curr_score)
    }
    
    # Create a ggplot  with the collected data
    df <- data.frame(Iteration = x_data, Score = y_data)
    
    plot <- ggplot2::ggplot(df, ggplot2::aes(x = Iteration, y = Score)) +
        ggplot2::geom_line(color = "blue") +
        ggplot2::geom_point(color = "blue") +
        ggplot2::labs(
            x = "Iteration", 
            y = "Score",
            title = "Score Optimization Over Iterations"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10)
        ) +
        scale_y_continuous(labels = scales::label_number(big.mark = ","))
    
    # Return the ggplot object
    return(list(plot = plot,
                final_position = pos))
}