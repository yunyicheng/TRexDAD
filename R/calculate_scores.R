#' @export
get_overhangs <- function(start, end, flag, target_gene) {
    head <- paste0(substr(target_gene[start - 2], nchar(target_gene[start - 2]), nchar(target_gene[start - 2])), target_gene[start - 1])
    tail <- paste0(target_gene[end + 1], substr(target_gene[end + 2], 1, 1))
    if (flag == 1) {
        head <- Biostrings::DNAString(head)
        tail <- Biostrings::DNAString(tail)
    }
    return(list(head = head, tail = tail))
}

#' @export
obtain_score <- function(start, end, length_tiles_global, target_gene, df) {
    overhangs <- get_overhangs(start, end, 1, target_gene)
    head <- overhangs$head
    tail <- overhangs$tail
    
    palindromicity <- 0
    if (Biostrings::reverseComplement(head) == head) palindromicity <- palindromicity + 1
    if (Biostrings::reverseComplement(tail) == tail) palindromicity <- palindromicity + 1
    
    delta_len <- end - start + 1 - length_tiles_global
    
    head_comp <- Biostrings::reverseComplement(head)
    tail_comp <- Biostrings::reverseComplement(tail)
    
    head_reactivity <- df[as.character(head), as.character(head_comp)]
    tail_reactivity <- df[as.character(tail), as.character(tail_comp)]
    
    on_target <- head_reactivity + tail_reactivity
    local_score <- on_target - palindromicity * 1000 - 5 * delta_len^2
    
    return(local_score)
}

calculate_scores <- function(pos_lst, target_gene, df, length_tiles_global) {
    overhangs <- get_all_overhangs(pos_lst, target_gene)
    complements <- lapply(overhangs, Biostrings::reverseComplement)
    
    off_reaction <- 0
    overhang_comb <- combn(overhangs, 2, simplify = FALSE)
    complement_comb <- combn(complements, 2, simplify = FALSE)
    
    for (oc in overhang_comb) {
        off_reaction <- off_reaction + df[as.character(oc[1]), as.character(oc[2])]
    }
    for (cc in complement_comb) {
        off_reaction <- off_reaction + df[as.character(cc[1]), as.character(cc[2])]
    }
    
    for (i in seq_along(overhangs)) {
        for (j in seq_along(complements)) {
            if (i != j) {
                off_reaction <- off_reaction + df[as.character(overhangs[[i]]), as.character(complements[[j]])]
            }
        }
    }
    
    repetition <- length(overhangs) - length(unique(overhangs))
    local_score <- sum(sapply(seq_along(pos_lst)[-length(pos_lst)], function(i) obtain_score(pos_lst[i], pos_lst[i + 1] - 1, length_tiles_global, target_gene, df)))
    
    global_score <- 50 * local_score - 100 * off_reaction - 1000 * repetition
    
    return(global_score)
}


