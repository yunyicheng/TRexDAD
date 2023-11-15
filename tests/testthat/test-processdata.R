test_that("split_into_codons splits sequence correctly", {
    gene_seq <- "ATGGCC"
    result <- split_into_codons(gene_seq)
    expect_equal(result, c("ATG", "GCC"))
})