
#' @export
split_into_codons <- function(gene_sequence) {
    codons <- strsplit(gene_sequence, "")[[1]]
    codons <- matrix(codons, nrow = 3)
    codons <- apply(codons, 2, paste, collapse = "")
    return(codons)
}

#' @export
oligo_cost <- function(num_tiles, num_codons) {
    ccdB_cost <- num_tiles * 16
    one_br_per_tile <- (num_codons * 3 / num_tiles + 30) * 0.01
    complete_replacement <- one_br_per_tile * num_codons
    total <- complete_replacement + ccdB_cost
    price_per_codon <- total / num_codons
    return(price_per_codon)
}

calculate_optimal_tiles <- function(num_codons) {
    num_tiles_range <- 1:50
    costs <- sapply(num_tiles_range, function(x) oligo_cost(x, num_codons))
    optimal_tiles <- which.min(costs)
    length_tiles_global <- round(num_codons / optimal_tiles)
    return(list(optimal_tiles = optimal_tiles, length_tiles = length_tiles_global))
}

# Assuming the file is already downloaded or accessible locally
read_overhang_fidelity_chart <- function(file_path) {
    if (!file.exists(file_path)) {
        stop("File does not exist")
    }
    df <- readxl::read_excel(file_path, sheet = "S1 Table. BsaI-HFv2")
    rownames(df) <- df$Overhang
    df$Overhang <- NULL
    return(df)
}

