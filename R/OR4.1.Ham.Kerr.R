#' @docType data
#' @keywords datasets
#' @title Example from Hamilton and Kerr
#' @description Pedigree file including diploid and triploid pedigree with normal meiosis and second division restitution.
#' @usage OR4.1.Ham.Kerr
#' @format Data frame with 8 rows and 10 columns
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton MG, Kerr RJ. Computation of the inverse additive relationship matrix for autopolyploid and multiple-ploidy populations
#' @references Hamilton MG, Kerr RJ. Computation of the inverse additive relationship matrix for autopolyploid and multiple-ploidy populations


#7-column example with diploids, triploids and second division restitution
OR4.1.Ham.Kerr <- data.frame(
  INDIV.ID           = c(1,	    2,	    3,	    4,	    5,	    6,	    7,	    8),
  SIRE.ID            = c(0,	    0,	    1,	    1,	    1,	    1,	    1,	    1),
  DAM.ID             = c(0,	    0,	    2,	    2,	    3,	    2,	    2,	    3),
  SIRE.GAMETE.PLOIDY = c(1,	    1,	    1,	    1,	    1,	    1,	    1,	    1),
  SIRE.LAMBDA        = c(0,	    0,	    0,	    0,	    0,	    0,	    0,	    0),
  DAM.GAMETE.PLOIDY  = c(1,	    1,	    1,	    1,	    1,	    2,	    2,	    2),
  DAM.LAMBDA         = c(0,	    0,	    0,	    0,	    0,	    0.918,  0.918,  0.918),
  SIRE.SEGREGATION   = c("Normal", "Normal", "Normal", "Normal", "Normal", "Normal", "Normal", "Normal"),
  DAM.SEGREGATION    = c("Normal", "Normal", "Normal", "Normal", "Normal", "Second division restitution", "Second division restitution", "Second division restitution"),
  INDIV.PLOIDY       = c("Diploid", "Diploid", "Diploid", "Diploid", "Diploid", "Triploid", "Triploid", "Triploid")
)

devtools::use_data(OR4.1.Ham.Kerr, overwrite = TRUE)

