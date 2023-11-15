test_that("split_into_codons splits sequence correctly", {
    gene_seq <- "ATGGCC"
    result <- split_into_codons(gene_seq)
    expect_equal(result, c("ATG", "GCC"))
})

test_that("oligo_cost returns correct values", {
    expect_equal(oligo_cost(10, 100), (100 * 3 / 10 + 30) * 0.01 * 100 + 10 * 16 / 100)
    expect_equal(oligo_cost(20, 200), (200 * 3 / 20 + 30) * 0.01 * 200 + 20 * 16 / 200)
    expect_equal(oligo_cost(5, 50), (50 * 3 / 5 + 30) * 0.01 * 50 + 5 * 16 / 50)
})

test_that("calculate_optimal_tiles returns correct values", {
    test_codons <- 100
    test_result <- calculate_optimal_tiles(test_codons)
    
    # Manually calculate expected values
    num_tiles_range <- 1:50
    costs <- sapply(num_tiles_range, function(x) oligo_cost(x, test_codons))
    expected_optimal_tiles <- which.min(costs)
    expected_length_tiles_global <- round(test_codons / expected_optimal_tiles)
    
    expect_equal(test_result$optimal_tiles, expected_optimal_tiles)
    expect_equal(test_result$length_tiles, expected_length_tiles_global)

})

