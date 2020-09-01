#' @docType data
#' @keywords datasets
#' @title Example from Kerr et al 2012
#' @description Pedigree file for tetraploid with double reduction in Table 4
#' @usage Tab.4.Kerr.et.al.2012
#' @format Data frame with 8 rows and 10 columns
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Kerr RJ, Li L, Tier B, Dutkowski GW, McRae TA (2012) Use of the numerator relationship matrix in genetic analysis of autopolyploid species. Theor Appl Genet 124:1271-1282. doi: DOI 10.1007/s00122-012-1785-y
#' @references Hamilton MG, Kerr RJ (2018) Computation of the inverse additive relationship matrix for autopolyploid and multiple-ploidy populations. 131:851-860. doi: 10.1007/s00122-017-3041-y


#7-column example with diploids, triploids and second division restitution

Tab.4.Kerr.et.al.2012 <- data.frame(
  INDIV.ID           = c(1,	    2,	    3,	    4,	    5,	    6,	    7,	  8,    9,    10),
  SIRE.ID            = c(0,	    0,	    0,	    1,	    1,	    3,	    5,	  6,    6,    6),
  DAM.ID             = c(0,	    0,	    0,	    2,	    4,	    5,	    6,	  7,    8,    8),
  SIRE.GAMETE.PLOIDY = c(2,	    2,	    2,	    2,	    2,	    2,	    2,	  2,    2,    2),
  SIRE.LAMBDA        = c(0.1,	  0.1,	  0.1,	  0.1,	  0.1,	  0.1,	  0.1,	0.1,  0.1,  0.1),
  DAM.GAMETE.PLOIDY  = c(2,	    2,	    2,	    2,	    2,	    2,	    2,	  2,    2,    2),
  DAM.LAMBDA         = c(0.1,	  0.1,	  0.1,	  0.1,	  0.1,	  0.1,	  0.1,	0.1,  0.1,  0.1),
  SIRE.SEGREGATION   = c("Normal", "Normal", "Normal", "Normal", "Normal", "Normal", "Normal", "Normal", "Normal", "Normal"),
  DAM.SEGREGATION    = c("Normal", "Normal", "Normal", "Normal", "Normal", "Normal", "Normal", "Normal", "Normal", "Normal"),
  INDIV.PLOIDY       = c("Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid")
)

usethis::use_data(Tab.4.Kerr.et.al.2012, overwrite = TRUE)

