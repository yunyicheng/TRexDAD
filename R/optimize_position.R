#' Pick a Position for Optimization
#'
#' This function selects a position from a list of positions for optimization 
#' based on their scores. Positions with lower scores are more likely to be picked.
#' It uses a weighted sampling approach where the weights are inversely proportional 
#' to the scores of the positions.
#'
#' @param pos_lst A numeric vector representing positions in a gene sequence.
#'
#' @return A list with two elements: `index`, the index in `pos_lst` of the picked 
#'         position (shifted by 1), and `position`, the actual value of the picked position.
#'
#' @examples
#' # Example usage assuming pos_lst and obtain_score function are defined
#' pos_lst <- c(1, 5, 10, 15)
#' picked <- pick_position(pos_lst)
#' print(picked)
#'
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

#' Optimize a Single Position in a Position List
#'
#' This function optimizes a single position within a list of positions. It adjusts
#' the specified position to maximize the overall score (obtained via `calculate_scores`).
#' The function supports different optimization modes, including greedy and Markov Chain
#' Monte Carlo (MCMC) approaches.
#'
#' @param pos_lst A numeric vector representing the current list of positions.
#' @param curr_pos The current position value that needs optimization.
#' @param left The left boundary of the scanning range for position optimization.
#' @param right The right boundary of the scanning range for position optimization.
#' @param mode The optimization mode: 1 for greedy, 2 for MCMC.
#'
#' @return The optimized position value within the specified range. If the function
#'         does not find a better position, it returns the original position value.
#'
#' @examples
#' # Example usage assuming pos_lst and calculate_scores function are defined
#' pos_lst <- c(1, 5, 10, 15)
#' optimized_pos <- optimize_position(pos_lst, pos_lst[2], 3, 7, 1)
#' print(optimized_pos)
#'
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

