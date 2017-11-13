#' @docType data
#' @keywords datasets
#' @title Example from Henderson (1988)
#' @description Pedigree file including diploids of uncertain parentage
#' @usage Henderson.1988
#' @format Data frame with 4 rows and 10 columns
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Henderson CR (1988) Use of an average numerator relationship matrix for multiple-sire joining. Journal of Animal Science 66, 1614-1621.
#' @references Hamilton MG, Kerr RJ. Computation of the inverse additive relationship matrix for autopolyploid and multiple-ploidy populations

#4-column example
Henderson.1988 <- data.frame(
  INDIV.ID    = c(1,   2,   3,   4,   5,   6,   7,   7,   8,   8,   9,   9,   9,   10),
  SIRE.ID     = c(0,   0,   1,   1,   3,   3,   3,   5,   1,   5,   1,   4,   5,   1),
  DAM.ID      = c(0,   0,   2,   2,   4,   0,   6,   6,   4,   4,   6,   6,   6,   4),
  PROBABILITY = c(1,   1,   1,   1,   1,   1,   0.6, 0.4, 0.3, 0.7, 0.3, 0.6, 0.1, 1)
)
devtools::use_data(Henderson.1988, overwrite = TRUE)
