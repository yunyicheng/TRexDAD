% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimize_position.R
\name{optimize_position}
\alias{optimize_position}
\title{Optimize a Single Position in a Position List}
\usage{
optimize_position(pos_lst, curr_pos, left, right, mode)
}
\arguments{
\item{pos_lst}{A numeric vector representing the current list of positions.}

\item{curr_pos}{The current position value that needs optimization.}

\item{left}{The left boundary of the scanning range for position optimization.}

\item{right}{The right boundary of the scanning range for position optimization.}

\item{mode}{The optimization mode: 1 for greedy, 2 for MCMC.}
}
\value{
The optimized position value within the specified range. If the function
does not find a better position, it returns the original position value.
}
\description{
This function optimizes a single position within a list of positions. It adjusts
the specified position to maximize the overall score (obtained via \code{calculate_scores}).
The function supports different optimization modes, including greedy and Markov Chain
Monte Carlo (MCMC) approaches.
}
\examples{
# Example usage assuming pos_lst and calculate_scores function are defined
pos_lst <- c(1, 5, 10, 15)
optimized_pos <- optimize_position(pos_lst, pos_lst[2], 3, 7, 1)
print(optimized_pos)

}