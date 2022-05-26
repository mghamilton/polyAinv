# polyAinv #

polyAinv is a R package (R Core Team 2016) comprised of a single function also named polyAinv.

### What is polyAinv ###

* The polyAinv function implements the methods of Hamilton and Kerr (2018) and Henderson (1988).  Specifically, it computes the additive relationship matrix (A), A inverse, coancestry matrix (K), K inverse and inbreeding coefficients (F) for autopolyploid and multiple-ploidy populations (Hamilton and Kerr 2018), and populations in which individuals have uncertain pedigree (Henderson 1988).
* For further details refer to the help page once you have installed polyAinv.

### To install polyAinv in R ###

* install.packages("devtools")
* library(devtools)
* install_github("mghamilton/polyAinv")
* library(polyAinv)
* help(polyAinv)

### Contact details ###

* <M.Hamilton@cgiar.org>

### References ###

* Hamilton MG, Kerr RJ (2018) Computation of the inverse additive relationship matrix for autopolyploid and multiple-ploidy populations. Theoretical and Applied Genetics 131:851–860. doi: 10.1007/s00122-017-3041-y
* Henderson CR (1988) Use of an average numerator relationship matrix for multiple-sire joining. Journal of Animal Science 66:1614-1621. doi: 10.2527/jas1988.6671614x
* R Core Team (2016) R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria

### Note ###
* Version 0.1.2 contained a bug causing identifiers in the 'INDIV.ID' column of the 'F' data frame to be erroneous in some circumstances. This has been rectified in version 0.1.3.  Thanks to Tim Millar for identifying this issue. 