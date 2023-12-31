% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline.R
\name{pick_position}
\alias{pick_position}
\title{Pick a Position for Optimization}
\usage{
pick_position(gene_codons, tile_length, pos_lst)
}
\arguments{
\item{pos_lst}{A numeric vector representing positions in a gene sequence.}
}
\value{
A list with two elements: \code{index}, the index in \code{pos_lst} of the picked
position (shifted by 1), and \code{position}, the actual value of the picked position.
}
\description{
This function selects a position from a list of positions for optimization
based on their scores. Positions with lower scores are more likely to be picked.
It uses a weighted sampling approach where the weights are inversely proportional
to the scores of the positions.
}
\examples{
# Example usage assuming pos_lst and obtain_score function are defined
pos_lst <- c(1, 5, 10, 15)
picked <- pick_position(pos_lst)
print(picked)

}
