#' Split Gene Sequence into Codons
#'
#' This function takes a gene sequence (as a string) and splits it into codons. 
#' A codon consists of three nucleotides, and this function organizes the 
#' sequence into these triplets.
#'
#' @param gene_sequence A character string representing the gene sequence.
#' 
#' @return A character vector where each element is a codon (a sequence of 
#'         three nucleotides).
#'
#' @examples
#' # Example gene sequence
#' gene_seq <- "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
#' split_into_codons(gene_seq)
#'
#' @export
split_into_codons <- function(gene_sequence) {
    codons <- strsplit(gene_sequence, "")[[1]]
    codons <- matrix(codons, nrow = 3)
    codons <- apply(codons, 2, paste, collapse = "")
    return(codons)
}


#' Calculate Oligo Cost
#'
#' This function calculates the cost of oligonucleotides (oligos) based on 
#' the number of tiles and the number of codons. It considers the cost of 
#' ccdB (a cloning vector) and the cost of base replacement per tile.
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

#' Calculate Optimal Number of Tiles
#'
#' This function determines the optimal number of tiles for a given number of codons. 
#' It calculates the oligo cost for a range of tile numbers and finds the number of 
#' tiles that minimizes this cost. It also computes the global length of the tiles.
#'
#' @param num_codons An integer representing the total number of codons.
#'
#' @return A list containing two elements: `optimal_tiles`, the number of tiles 
#'         that minimize the oligo cost, and `length_tiles_global`, the average 
#'         length of these tiles.
#'
#' @examples
#' # Calculate optimal tiles for 100 codons
#' calculate_optimal_tiles(100)
#'
#' @export
calculate_optimal_tiles <- function(num_codons) {
    num_tiles_range <- 1:50
    costs <- sapply(num_tiles_range, function(x) oligo_cost(x, num_codons))
    optimal_tiles <- which.min(costs)
    length_tiles_global <- round(num_codons / optimal_tiles)
    return(list(optimal_tiles = optimal_tiles, length_tiles = length_tiles_global))
}
