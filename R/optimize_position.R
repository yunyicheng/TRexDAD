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
