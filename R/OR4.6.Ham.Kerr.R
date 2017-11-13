#' @docType data
#' @keywords datasets
#' @title Example from Hamilton and Kerr
#' @description Pedigree file including diploid and tetraploid pedigree with normal meiosis, first division restitution and somatic chromosome doubling.
#' @usage OR4.6.Ham.Kerr
#' @format Data frame with 8 rows and 10 columns
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton MG, Kerr RJ. Computation of the inverse additive relationship matrix for autopolyploid and multiple-ploidy populations
#' @references Hamilton MG, Kerr RJ. Computation of the inverse additive relationship matrix for autopolyploid and multiple-ploidy populations

#7-column example with diploids, tetraploids, first division restitution and somatic chromosome doubling
OR4.6.Ham.Kerr <- data.frame(
  INDIV.ID           = c(1,	    2,	    3,	    4,	    5,	    6,	    7,	    8),
  SIRE.ID            = c(0,	    0,	    0,	    1,	    1,	    0,	    0,	    5),
  DAM.ID             = c(0,	    0,	    2,	    0,	    3,	    3,	    4,	    0),
  SIRE.GAMETE.PLOIDY = c(1,	    2,	    0,	    2,	    1,	    0,	    0,	    4),
  SIRE.LAMBDA        = c(0,	    0.167,  0,      0.041,  0,      0,      0,      0.333),
  DAM.GAMETE.PLOIDY  = c(1,	    2,	    2,	    0,	    1,	    4,	    4,	    0),
  DAM.LAMBDA         = c(0,	    0.167,  0.167,  0,      0,      0.333,  0.333,  0),
  SIRE.SEGREGATION   = c("Normal", "Normal", "NA", "First division restitution", "Normal", "NA", "NA", "Somatic chromosome doubling"),
  DAM.SEGREGATION    = c("Normal", "Normal", "Normal", "NA", "Normal", "Somatic chromosome doubling", "Somatic chromosome doubling", "NA"),
  INDIV.PLOIDY       = c("Diploid", "Tetraploid", "Diploid", "Diploid", "Diploid", "Tetraploid", "Tetraploid", "Tetraploid")
)
devtools::use_data(OR4.6.Ham.Kerr, overwrite = TRUE)
