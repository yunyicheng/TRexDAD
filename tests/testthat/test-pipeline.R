library(testthat)

# Test split_into_codons function

test_that("split_into_codons returns correct codons for valid input", {
    gene_seq <- "GCTATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATC"
    expected_codons <- c("GCT", "ATG", "GCC", "ATT", "GTA", "ATG", "GGC", "CGC", "TGA", "AAG", "GGT", "GCC", "CGA", "TAG", "ATC")
    expect_equal(split_into_codons(gene_seq), expected_codons)
})

test_that("split_into_codons throws an error for invalid nucleotides", {
    invalid_seq <- "GCTXATG"
    expect_error(split_into_codons(invalid_seq), "Invalid gene sequence: sequence should only contain A, C, G, and T.")
})

test_that("split_into_codons throws an error for sequences not multiple of 3", {
    invalid_length_seq <- "CGGATGA"
    expect_error(split_into_codons(invalid_length_seq), "Invalid gene sequence length: length should be a multiple of 3.")
})

test_that("split_into_codons throws an error if 'ATG' is not the second codon", {
    invalid_atg_seq <- "GCAGCCATG"
    expect_error(split_into_codons(invalid_atg_seq), "Invalid sequence: 'ATG' is not the second codon.")
})


# Test oligo_cost function

test_that("oligo_cost calculates correctly for normal input", {
    expect_equal(oligo_cost(10, 100), 2.2)
})

test_that("oligo_cost handles edge cases", {
    # Example edge case with 1 tile and 1 codon
    expect_equal(oligo_cost(1, 1), 16.33)
})

test_that("oligo_cost throws an error for zero or negative tiles", {
    expect_error(oligo_cost(0, 100), "The number of tiles should be a positive integer.")
    expect_error(oligo_cost(-1, 100), "The number of tiles should be a positive integer.")
})

test_that("oligo_cost throws an error for zero or negative codons", {
    expect_error(oligo_cost(10, 0), "The number of codons should be a positive integer.")
    expect_error(oligo_cost(10, -5), "The number of codons should be a positive integer.")
})