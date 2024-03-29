library(testthat)

# --- SECTION: Test Preprocess Data -----------------

## Test split_into_codons function ----

test_that("split_into_codons returns correct codons for valid input", {
    gene_seq <- "GCTATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATC"
    expected_codons <- c("GCT", "ATG", "GCC", "ATT", "GTA", "ATG", "GGC", 
                         "CGC", "TGA", "AAG", "GGT", "GCC", "CGA", "TAG", "ATC")
    expect_equal(split_into_codons(gene_seq), expected_codons)
})

test_that("split_into_codons throws an error for invalid nucleotides", {
    invalid_seq <- "GCTXATG"
    expect_error(split_into_codons(invalid_seq), 
        "Invalid gene sequence: sequence should only contain A, C, G, and T.")
})

test_that("split_into_codons throws an error for sequences not multiple of 3", {
    invalid_length_seq <- "CGGATGA"
    expect_error(split_into_codons(invalid_length_seq), 
            "Invalid gene sequence length: length should be a multiple of 3.")
})

test_that("split_into_codons throws an error if 'ATG' is not the second codon", {
    invalid_atg_seq <- "GCAGCCATG"
    expect_error(split_into_codons(invalid_atg_seq), 
                 "Invalid sequence: 'ATG' is not the second codon.")
})


## Test oligo_cost function ----

test_that("oligo_cost calculates correctly for normal input", {
    expect_equal(oligo_cost(10, 100), 2.2)
})

test_that("oligo_cost handles edge cases", {
    # Example edge case with 1 tile and 1 codon
    expect_equal(oligo_cost(1, 1), 16.33)
})

test_that("oligo_cost throws an error for zero or negative tiles", {
    expect_error(oligo_cost(0, 100), 
                 "The number of tiles should be a positive integer.")
    expect_error(oligo_cost(-1, 100), 
                 "The number of tiles should be a positive integer.")
})

test_that("oligo_cost throws an error for zero or negative codons", {
    expect_error(oligo_cost(10, 0), 
                 "The number of codons should be a positive integer.")
    expect_error(oligo_cost(10, -5), 
                 "The number of codons should be a positive integer.")
})


## Test calculate_optimal_tiles function ----

test_that("calculate_optimal_tiles returns correct structure", {
    result <- calculate_optimal_tiles(100)
    
    # Check if it's a list and has the correct elements
    expect_type(result, "list")
    expect_true(all(
        c("optimal_num_tiles", "optimal_tile_length", "optimal_cost")
        %in% names(result)))
})

test_that("calculate_optimal_tiles returns valid values", {
    result <- calculate_optimal_tiles(100)
    
    # Validate the types of the list elements
    expect_type(result$optimal_tile_length, "double")
    expect_type(result$optimal_cost, "double")
    
    # Check if the values make sense
    expect_true(result$optimal_num_tiles > 0)
    expect_true(result$optimal_tile_length > 0)
    expect_true(result$optimal_cost > 0)
})

test_that("calculate_optimal_tiles handles invalid input", {
    # Example with zero codons
    expect_error(calculate_optimal_tiles(0), 
                 "The number of codons should be a positive integer.")
    
    # Example with negative codons
    expect_error(calculate_optimal_tiles(-100),
                 "The number of codons should be a positive integer.")
})


# --- SECTION: Test Calculate Scores -----------------

## Test get_overhangs function ----

test_that("get_overhangs calculates correct sequences", {
    gene_codons <- c("GAG", "CTG", "TGT", "AGG", "TGC", "CGG", "CCA", 
                     "ATT", "TGA", "TAG", "GAA", "TAT", "AGC")
    result <- get_overhangs(gene_codons, 3, 4, FALSE)
    
    # Check if head and tail sequences are correct
    # head: last 1 of GAG (start - 2) + first 3 of CTG (start - 1) = "GCTG"
    # tail: last 3 of TGC (end + 1) + first 1 of CGG (end + 2) = "TGCC"
    expect_equal(result$head, "GCTG")
    expect_equal(result$tail, "TGCC")
})

test_that("get_overhangs returns DNAString objects when as_dna_strings is TRUE", {
    gene_codons <- c("GAG", "CTG", "TGT", "AGG", "TGC", "CGG", "CCA", 
                     "ATT", "TGA", "TAG", "GAA", "TAT", "AGC")
    result <- get_overhangs(gene_codons, 3, 4, TRUE)
    
    # Check if head and tail are DNAString objects
    expect_s4_class(result$head, "DNAString")
    expect_s4_class(result$tail, "DNAString")
})

test_that("get_overhangs returns character strings when as_dna_strings is FALSE", {
    gene_codons <- c("GAG", "CTG", "TGT", "AGG", "TGC", "CGG", "CCA", 
                     "ATT", "TGA", "TAG", "GAA", "TAT", "AGC")
    result <- get_overhangs(gene_codons, 3, 4, FALSE)
    
    # Check if head and tail are character strings
    expect_type(result$head, "character")
    expect_type(result$tail, "character")
})

test_that("get_overhangs handles invalid input", {
    # 'X' is an invalid nucleotide
    gene_codons <- c("GAG", "CTG", "XTG", "AGT", "GAA", "TGG")
    expect_error(get_overhangs(gene_codons, 3, 4, FALSE), 
                 "Invalid codon detected.")
    
    # Test for out of bounds start and end positions
    gene_codons <- c("GAG", "CTG", "TGT")
    # start <= 2
    expect_error(get_overhangs(gene_codons, 0, 4, FALSE), 
                 "Invalid start or end positions.")
    # end + 2 > length(gene_codons)
    expect_error(get_overhangs(gene_codons, 2, 4, FALSE), 
                 "Invalid start or end positions.")  
})


## Test get_all_overhangs function ----

test_that("get_all_overhangs returns correct overhang sequences", {
    gene_codons <- c("GAG", "CTG", "TGT", "AGG", "TGC", "CGG", "CCA", 
                     "ATT", "TGA", "TAG", "GAA", "TGA", "AGC")
    pos_lst <- c(3, 5, length(gene_codons) - 2)
    
    result <- get_all_overhangs(gene_codons, pos_lst)
    
    # Check the size of the result list
    expect_length(result, (length(pos_lst) - 1) * 2)
    
    # Check if the overhang sequences are DNAString objects
    expect_true(all(sapply(result, class) == "DNAString"))
    
    # Checks for the overhang sequences based on gene_codons and pos_lst
    expect_equal(result[[1]], DNAString("GCTG"))
    expect_equal(result[[2]], DNAString("TGCC"))
    expect_equal(result[[3]], DNAString("TAGG"))
    expect_equal(result[[4]], DNAString("TGAA"))
})

test_that("get_all_overhangs handles invalid input", {
    gene_codons <- c("GAG", "CTG", "TGT", "AGG", "TGC", "CGG", "CCA", 
                     "ATT", "TGA", "TAG", "GAA", "TGA", "AGC")
    
    # Invalid because it contains non-positive integer
    expect_error(get_all_overhangs(gene_codons, c(-3, 8, 11)), 
                 "pos_lst must contain only positive integers.")
    
    # Invalid because the last position is not the length of gene_codons - 2
    expect_error(get_all_overhangs(gene_codons, c(3, 8, 12), 
                                   "The last position in pos_lst must 
                                   equal length(gene_codons) - 2."))
    
    # Invalid because pos_lst is not increasing
    expect_error(get_all_overhangs(gene_codons, c(11, 11, 11),
                                   "pos_lst must be in increasing order."))
    
    # Invalid because pos_lst has only two positions
    expect_error(get_all_overhangs(gene_codons, c(8, 11),
                            "pos_lst must contain at least three positions."))
})

## Test calculate_local_score function ----

test_that("calculate_local_score calculates correct score", {
    gene_codons <- c("ATG", "CAG", "TAC", "GGA", "TGA", "TCA")
    tile_length <- 4
    start <- 3
    end <- 4
    result <- calculate_local_score(gene_codons, tile_length, start, end)
    expected_score <- 988
    
    expect_equal(result, expected_score)
})

test_that("calculate_local_score handles invalid codons", {
    gene_codons <- c("ATG", "CAG", "XXX", "GGA", "XXY", "TCA")
    tile_length <- 4
    start <- 3
    end <- 4
    
    expect_error(calculate_local_score(gene_codons, tile_length, start, end))
})

test_that("calculate_local_score handles invalid parameters", {
    gene_codons <- c("ATG", "CAG", "TAC", "GGA", "TGA", "TCA")
    tile_length <- 4
    start <- 3
    end <- 4
    
    # Invalid tile_length
    expect_error(calculate_local_score(gene_codons, -1, start, end))
    # Invalid start
    expect_error(calculate_local_score(gene_codons, tile_length, 0, end))
    # Invalid end
    expect_error(calculate_local_score(gene_codons, tile_length, start, 6))
})


## Test calculate_global_score function

# library(testthat)
# library(Biostrings)
# library(utils)

## Test calculate_global_score function----

test_that("calculate_global_score calculates correct global score", {
    gene_codons <- c("GAG", "ATG", "TGT", "AGG", "TGC", "CGG", "CCA", 
                     "ATT", "TGA", "TAG", "GAA", "TGA", "AGC")
    tile_length <- 5
    pos_lst <- c(3, 6, 8, 11)
    
    expected_global_score <- 122700
    result <- calculate_global_score(gene_codons, tile_length, pos_lst)
    
    expect_equal(result, expected_global_score)
})

test_that("calculate_local_score handles invalid codons", {
    gene_codons <- c("ATG", "CAG", "XXX", "GGA", "XXY", "TCA")
    tile_length <- 5
    expect_error(calculate_global_score(gene_codons, tile_length, c(3, 8, 11)))
})

test_that("calculate_global_score handles invalid pos_lst and tile length", {
    gene_codons <- c("ATG", "CAG", "TAC", "GGA", "TCA")
    tile_length <- 5
    
    # Invalid start
    expect_error(calculate_global_score(gene_codons, tile_length, c(0, 2, 6)))
    # Invalid end
    expect_error(calculate_global_score(gene_codons, tile_length, c(3, 8, 12)))
    # Not increasing
    expect_error(calculate_global_score(gene_codons, tile_length, c(11, 8, 3)))
    
    # Invalid tile_length
    expect_error(calculate_global_score(gene_codons, -1, c(1, 3, 5)))
})


# --- SECTION: Test Optimize Positions -----------------

## Test pick_position function ----

test_that("pick_position selects a valid position", {
    gene_codons <- c("GAG", "ATG", "TGT", "AGG", "TGC", "CGG", "CCA", 
                     "ATT", "TGA", "TAG", "GAA", "TGA", "AGC")
    tile_length <- 5
    pos_lst <- c(3, 6, 8, 11)
    
    # Run the function
    picked <- pick_position(gene_codons, tile_length, pos_lst)
    
    # Check if the picked position is within the range of pos_lst
    expect_true(picked$index %in% seq_along(pos_lst))
    expect_true(picked$position %in% pos_lst)
})

test_that("pick_position handles edge cases and invalid input", {
    gene_codons <- c("GAG", "ATG", "TGT", "AGG", "TGC", "CGG", "CCA", 
                     "ATT", "TGA", "TAG", "GAA", "TGA", "AGC")
    tile_length <- 5
    
    # Test with minimal length pos_lst
    pos_lst_minimal <- c(3, 6, 11)
    expect_silent(pick_position(gene_codons, tile_length, pos_lst_minimal))
    
    # Test with invalid pos_lst (too short)
    pos_lst_short <- c(1)  # Too short for the function to work
    expect_error(pick_position(gene_codons, tile_length, pos_lst_short))
    
})

# Test optimize_position function ----

test_that("optimize_position returns a value within specified range", {
    gene_codons <- c("GAG", "ATG", "TGT", "AGG", "TGC", "CGG", "CCA", 
                     "ATT", "TGA", "TAG", "GAA", "TGA", "AGC")
    tile_length <- 5
    pos_lst <- c(3, 6, 11)
    curr_pos <- pos_lst[2]
    left <- 7
    right <- 10
    
    optimized_pos_greedy <- optimize_position(gene_codons, tile_length,
                                              pos_lst, curr_pos, 
                                              left, right, TRUE)
    optimized_pos_mcmc <- optimize_position(gene_codons, tile_length, 
                                            pos_lst, curr_pos, 
                                            left, right, FALSE)
    
    # Check if the optimized positions are within the range
    expect_true(optimized_pos_greedy %in% left:right)
    expect_true(optimized_pos_mcmc %in% left:right)
})

test_that("optimize_position handles edge cases", {
    gene_codons <- c("GAG", "ATG", "TGT", "AGG", "TGC", "CGG", "CCA", 
                     "ATT", "TGA", "TAG", "GAA", "TGA", "AGC")
    tile_length <- 5
    pos_lst <- c(3, 6, 11)
    curr_pos <- pos_lst[2]
    left <- 6  # Left boundary same as current position
    right <- 6  # Right boundary same as current position
    
    # Test should pass without throwing an error and return the same position
    expect_equal(optimize_position(gene_codons, tile_length, pos_lst, 
                                   curr_pos, left, right, TRUE), curr_pos)
    expect_equal(optimize_position(gene_codons, tile_length, pos_lst, 
                                   curr_pos, left, right, FALSE), curr_pos)
})


# Test execute_and_plot function ----

# Helper function to randomly generate testing examples
generate_random_dna_sequence <- function(num_codons) {
    if (num_codons <= 1) {
        stop("Number of codons must be greater than 1 to include
             'ATG' as the second codon.")
    }
    
    # Generate the first codon randomly
    first_codon <- paste(sample(c("A", "T", "G", "C"), 
                                size = 3, replace = TRUE),
                         collapse = "")
    
    # Generate the rest of the sequence randomly, excluding the first and 
    # second codon
    remaining_sequence <- paste(sample(c("A", "T", "G", "C"), 
                                       size = (num_codons - 2) * 3, 
                                       replace = TRUE), 
                                collapse = "")
    
    # Concatenate the first codon, 'ATG', and the remaining sequence
    full_sequence <- paste(first_codon, "ATG", remaining_sequence, sep = "")
    
    return(full_sequence)
}

test_that("execute_and_plot completes without error", {
    # Generate a random gene sequence with a length of at least 60 codons
    test_gene_1 <- generate_random_dna_sequence(num_codons = 60)
    test_gene_2 <- generate_random_dna_sequence(num_codons = 60)
    test_gene_3 <- generate_random_dna_sequence(num_codons = 60)
    max_iter <- 5  # Reduced number for test efficiency
    scan_rate <- 2
    
    expect_no_error(execute_and_plot(test_gene_1, max_iter, scan_rate))
    expect_no_error(execute_and_plot(test_gene_2, max_iter, scan_rate))
    expect_no_error(execute_and_plot(test_gene_3, max_iter, scan_rate))
})

test_that("execute_and_plot handles invalid input", {
    # Generate an invalid gene sequence
    invalid_gene <- paste(generate_random_dna_sequence(num_codons = 60), 
                          "XXX", sep = "")
    max_iter <- 5
    scan_rate <- 2
    
    expect_error(execute_and_plot(invalid_gene, max_iter, scan_rate))
})