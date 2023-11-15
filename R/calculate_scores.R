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
#' target_gene <- c("ATG", "GCT", "TGA", "AGT", "CTA")
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
#' @examples
#' # Assuming pos_lst and target_gene are predefined
#' pos_lst <- c(1, 5, 10)
#' target_gene <- c("ATG", "GCT", "TGA", "AGT", "CTA")
#' get_all_overhangs(pos_lst)
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

