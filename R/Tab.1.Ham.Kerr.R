#' @docType data
#' @keywords datasets
#' @title Example from Hamilton and Kerr 2018
#' @description Pedigree file including diploid and tetraploid pedigree with normal meiosis, first division restitution, second division restitution and parthenogenesis.
#' @usage Tab.1.Ham.Kerr
#' @format Data frame with 8 rows and 10 columns
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton MG, Kerr RJ (2018) Computation of the inverse additive relationship matrix for autopolyploid and multiple-ploidy populations. 131:851-860. doi: 10.1007/s00122-017-3041-y
#' @references Hamilton MG, Kerr RJ (2018) Computation of the inverse additive relationship matrix for autopolyploid and multiple-ploidy populations. 131:851-860. doi: 10.1007/s00122-017-3041-y

#7-column example with diploids, tetraploids, first division restitution, second division restitution and parthenogenesis
Tab.1.Ham.Kerr <- data.frame(
  INDIV.ID           = c(1,	    2,	    3,	    4,	    5,	    6,	    7,	    8),
  SIRE.ID            = c(0,	    0,	    0,	    1,	    1,	    1,	    6,	    6),
  DAM.ID             = c(0,	    0,	    2,	    0,	    3,	    3,	    2,	    2),
  SIRE.GAMETE.PLOIDY = c(1,	    2,	    0,	    2,	    1,	    2,	    2,	    2),
  SIRE.LAMBDA        = c(0,	    0.167,  0,      0.041,  0,      0.918,  0.167,  0.167),
  DAM.GAMETE.PLOIDY  = c(1,	    2,	    2,	    0,	    1,	    2,	    2,	    2),
  DAM.LAMBDA         = c(0,	    0.167,  0.167,  0,      0,      0.041,  0.167,  0.167),
  SIRE.SEGREGATION   = c("Normal", "Normal", "NA", "First division restitution", "Normal", "Second division restitution", "Normal", "Normal"),
  DAM.SEGREGATION    = c("Normal", "Normal", "Normal", "NA", "Normal", "First division restitution", "Normal", "Normal"),
  INDIV.PLOIDY       = c("Diploid", "Tetraploid", "Diploid", "Diploid", "Diploid", "Tetraploid", "Tetraploid", "Tetraploid")
)

usethis::use_data(Tab.1.Ham.Kerr, overwrite = TRUE)
