#' polyAinv
#'
#' Computes the additive relationship matrix (A), A inverse,
#' coancestry matrix (K), K inverse and inbreeding coefficients (F) for autopolyploid and multiple-ploidy populations
#' (Hamilton and Kerr 2018), and populations in which individuals have uncertain pedigree (Henderson 1988).
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @param ped A 3-, 4-, 7- or 8-column dataframe containing data of class numeric and/or integer.
#' @details Columns in the 3-, 4-, 7- or 8-column 'ped' dataframe must be in the following order:
#' @details - Individual identifier (required)
#' @details - Sire identifier (required)
#' @details - Dam identifier (required)
#' @details - Sire gamete ploidy level (lowercase tau in Hamilton and Kerr 2018; required in 7- and 8-column dataframes only)
#' @details - Sire lambda - the probability that two genes, drawn at random from a random locus in the sire gamete, are identical
#' by decent as a result of two copies of a gene being inherited from a single parental chromosome.  This may
#' be greater than zero in non-monoploid gametes through double reduction, irregular meiosis (e.g. the formation
#' of unreduced gametes through first division restitution or second division restitution), double fertilisation,
#' or irregular somatic chromosome doubling. (Hamilton and Kerr 2018; required in 7- and 8-column dataframes only)
#' @details - Dam gamete ploidy level (lowercase tau in Hamilton and Kerr 2018; required in 7- and 8-column dataframes only)
#' @details - Dam lambda - the probability that two genes, drawn at random from a random locus in the dam gamete, are identical
#' by decent as a result of two copies of a gene being inherited from a single parental chromosome.  This may
#' be greater than zero in non-monoploid gametes for the reasons outlined for 'Sire lambda' above (Hamilton and Kerr 2018;
#' required in 7- and 8-column dataframes only)
#' @details - Probability that the individual is the progeny of the sire and dam (Henderson 1988; required in 4- and 8-column dataframes only)
#' @details Note:
#' @details - column headings in 'ped' are ignored
#' @details - if 'ped' contains 3-columns it is assumed to be a standard diploid pedigree file
#' @details - if 'ped' contains 4-columns it is assumed to be a diploid pedigree file in which individuals have uncertain pedigree (Henderson 1988)
#' @details - if 'ped' contains 7-columns it is assumed to be a polyploid pedigree (Hamilton and Kerr 2018)
#' @details - if 'ped' contains 8-columns it is assumed to be a polyploid pedigree in which individuals have uncertain pedigree (Hamilton and Kerr 2018; Henderson 1988)
#' @param ASReml.giv.only Logical.  If TRUE only ASReml.giv is returned.  This can reduce memory requirements.
#' @return A list with the following objects (Hamilton and Kerr 2018):
#' @return - K.parents:  Coancestry (K) matrix including parents only
#' @return - K.inv:      Inverse coancestry (K) triangular matrix as a 3-column (triplet) dataframe with zero elements removed
#' @return - A.parents:  Numerator relationship (A) matrix including parents only
#' @return - A.inv:      Inverse numerator relationship (A) triangular matrix as a 3-column (triplet) dataframe with zero elements removed
#' @return - ASReml.giv: A.inv as an ASReml-R 'ginverse' object (see 'asreml.Ainverse' in Butler et al. 2009 and '.giv' file in Gilmour et al. 2014)
#' @return - F:          Dataframe with 5 columns: INDIV.ID (Individual identifier), PLOIDY (Individual ploidy level, e.g. 2 for diploids), A.DIAG (diagonal elements of the A matrix), K.DIAG (diagonal elements of the K matrix), F (Coefficient of inbreeding)
#' @references Amadeu RR, Cellon C, Olmstead JW, Garcia AAF, Resende MFR, Munoz PR (2016) AGHmatrix: R package to construct relationship matrices for autotetraploid and diploid species. A blueberry example. The Plant Genome 9
#' @references Butler DG, Cullis BR, Gilmour AR, Gogel BJ (2009) ASReml-R reference manual. Queensland Department of Primary Industries and Fisheries, Brisbane
#' @references Gilmour AR, Gogel BJ, Cullis BR, Welham SJ, Thompson R (2014) ASReml user guide release 4.1 functional specification. VSN International Ltd, Hemel Hempstead, UK
#' @references Hamilton MG, Kerr RJ (2018) Computation of the inverse additive relationship matrix for autopolyploid and multiple-ploidy populations. 131:851-860. doi: 10.1007/s00122-017-3041-y
#' @references Henderson CR (1988) Use of an average numerator relationship matrix for multiple-sire joining. Journal of Animal Science 66, 1614-1621.
#' @examples
#' #7-column example with diploids, tetraploids, first division restitution,
#' #second division restitution and parthenogenesis
#' data(Tab.1.Ham.Kerr)
#' polyAinv(ped = Tab.1.Ham.Kerr[,1:7])
#'
#' #7-column example with diploids, triploids and second division restitution
#' data(OR4.1.Ham.Kerr)
#' polyAinv(ped = OR4.1.Ham.Kerr[,1:7])
#'
#' #7-column example with diploids, tetraploids, first division restitution and
#' #somatic chromosome doubling.  Returns an ASReml ginverse object only.
#' data(OR4.6.Ham.Kerr)
#' polyAinv(ped = OR4.6.Ham.Kerr[,1:7], ASReml.giv.only = TRUE)
#'
#' #4-column example from Henderson (1988)
#' data(Henderson.1988)
#' polyAinv(ped = Henderson.1988[,1:4])
#'
#' #8-column example (OR4.1.Ham.Kerr example pedigree modified)
#' data(OR4.1.Ham.Kerr.uncert)
#' polyAinv(ped = OR4.1.Ham.Kerr.uncert[,c(1:7,11)])
#'
#' #All tetraploid example with double reduction in founders (Tab.4.Kerr.et.al.2012).  
#' #Note that results differ from Table 4 of Kerr et al. (2012) and those using the 
#' #approach of Amadeu (2016), as inbreeding in founders is accounted for.
#' data(Tab.4.Kerr.et.al.2012)
#' polyAinv(ped = Tab.4.Kerr.et.al.2012[,c(1:7)])
#' @export

polyAinv <- function(ped, ASReml.giv.only = FALSE) {

  print(Sys.time())

  #########################################################################################
  # Define sub functions
  #########################################################################################

  # Equation 12 of Hamilton and Kerr 2018
  get.founder.F <- function(indiv.ploidy, sire.lambda, dam.lambda) {
    F <- (indiv.ploidy/2 - 1) * (sire.lambda + dam.lambda) /
      (2 * indiv.ploidy/2 + indiv.ploidy/2 *
         (sire.lambda + dam.lambda) - (sire.lambda + dam.lambda))
    return(F)
  }

  get.p.and.q <- function(ped, indiv.parents, K.parents, indiv.count) {

    p.q.dims <- min(nrow(K.parents), indiv.count - 1)

    # Compute p.q, F.parent and parent.pr.ibd ###################################

    tmp.p.q <- aggregate(indiv.parents[, "PROBABILITY"], by = list(indiv.parents[, "PARENT.ID"]), na.rm = TRUE, FUN = "sum")
    colnames(tmp.p.q) <- c("PARENT.ID", "PROBABILITY.PARENT")
    indiv.parents <- merge(indiv.parents, tmp.p.q, by = "PARENT.ID", all.x = TRUE)

    tmp.p.q <- unique(indiv.parents[, c("PARENT.ID", "PARENT.GAMETE.PLOIDY", "OTHER.PARENT.GAMETE.PLOIDY", "PROBABILITY.PARENT")])
    p.q.multiplier <- tmp.p.q[1, "PARENT.GAMETE.PLOIDY"] /
      (tmp.p.q[1, "PARENT.GAMETE.PLOIDY"] + tmp.p.q[1, "OTHER.PARENT.GAMETE.PLOIDY"])
    tmp.p.q <- tmp.p.q[, "PROBABILITY.PARENT"]
    p.q <- data.frame(PARENT.ID = unique(indiv.parents[, "PARENT.ID"]), p.q = tmp.p.q)
    rm(tmp.p.q)

    if (indiv.parents[1, "PARENT.ID"] == 0) {
      # if parent is a founder
      F.parent <- get.founder.F(indiv.ploidy  = indiv.parents[1, "PARENT.GAMETE.PLOIDY"] * 2,
                                sire.lambda = founder.lambda[founder.lambda[, "INDIV.PLOIDY"] ==
                                                               indiv.parents[1, "PARENT.GAMETE.PLOIDY"] * 2, "SIRE.LAMBDA"],
                                dam.lambda    = founder.lambda[founder.lambda[, "INDIV.PLOIDY"] ==
                                                                 indiv.parents[1, "PARENT.GAMETE.PLOIDY"] * 2, "DAM.LAMBDA"])
    } else {
      # F of parent 8v.p p.i' K.(i-1) p.i^-1 / 2v.p-1 from rhs of Equation 15-16 of Hamilton and Kerr 2018
      v.p.q <- ped[ped[, "INDIV.ID"] %in%
                     indiv.parents[1, "PARENT.ID"], "INDIV.PLOIDY"][1]/2 #assumes all possible parents of same ploidy

      F.parent <- (2 * v.p.q * as.numeric(t(as.matrix(p.q[, "p.q"])) %*% K.parents[p.q[, "PARENT.ID"], p.q[, "PARENT.ID"]] %*%
                                            as.matrix(p.q[, "p.q"])) - 1)/(2 * v.p.q - 1)
    }  #END if(indiv.parents[1,'PARENT.ID'] == 0) {

    if (sum(indiv.parents[, "PARENT.GAMETE.PLOIDY"]) == 0) {
      F.parent <- 0 # if there is no gametic contribution then let F.parent = 0
    }

    # Equation 15-16 of Hamilton and Kerr 2018
    parent.pr.ibd <- indiv.parents[1, "PARENT.LAMBDA"] + (1 - indiv.parents[1, "PARENT.LAMBDA"]) * F.parent

    return(list(parent.pr.ibd=parent.pr.ibd, p.q = p.q, p.q.multiplier = p.q.multiplier))
  }

  #########################################################################################
  #Check pedigree file
  #########################################################################################

  print("Checking pedigree file")

  orig.col.names <- colnames(ped)

  colnames(ped) <- c("INDIV.ID", "SIRE.ID", "DAM.ID")

  indiv.ids <- unique(ped[, "INDIV.ID"])

  parent.ids <- unique(c(ped[, "SIRE.ID"], ped[, "DAM.ID"]))
  parent.ids <- parent.ids[order(parent.ids, decreasing = FALSE)]  #sort
  parent.ids <- parent.ids[parent.ids != 0]  #remove 0

  #Check that all parents are present in the indiv.ids
  if (sum(!(parent.ids %in% indiv.ids)) > 0) {
    stop("NOT ALL PARENTS ARE PRESENT AS INDIVIDUALS IN THE PEDIGREE FILE")
  }
  rm(indiv.ids)

  #Check that all 'founders' are parents in pedigree
  if(sum((!ped[,"INDIV.ID"] %in% parent.ids) & ped[,"SIRE.ID"] == 0 & ped[,"DAM.ID"] == 0)  > 0) {
    stop("NOT ALL FOUNDERS ARE LISTED AS PARENTS IN THE PEDIGREE FILE.  REMOVE FOUNDERS THAT ARE NOT PARENTS")
  }

  # Check that there are 3, 4, 7 or 8 columns in ped
  if (!(ncol(ped) == 3 | ncol(ped) == 4 | ncol(ped) == 7 | ncol(ped) == 8)) {
    stop("THE ped FILE DOES NOT CONTAIN 3,4, 7 or 8 COLUMNS")
  }

  #Generate 8-column pedigree file from 3-, 4- or 7- column files
  if(ncol(ped) == 3) { #assume standard diploid pedigree file
    ped[,"SIRE.GAMETE.PLOIDY"] <- 1
    ped[,"SIRE.LAMBDA"]        <- 0
    ped[,"DAM.GAMETE.PLOIDY"]  <- 1
    ped[,"DAM.LAMBDA"]         <- 0
    ped[,"PROBABILITY"]        <- 1
  }

  if(ncol(ped) == 4) { #assume diploid pedigree file with aggregates (Henderson 1988)
    ped[,"SIRE.GAMETE.PLOIDY"] <- 1
    ped[,"SIRE.LAMBDA"]        <- 0
    ped[,"DAM.GAMETE.PLOIDY"]  <- 1
    ped[,"DAM.LAMBDA"]         <- 0
    ped <- ped[,c(1, 2, 3, 5, 6, 7, 8, 4)] #reorder to ensure PROBABILITY is last column
  }

  if(ncol(ped) == 7) { #assume polyploid pedigree without aggregates
    ped[,"PROBABILITY"]        <- 1
  }

  colnames(ped) <- c("INDIV.ID", "SIRE.ID", "DAM.ID", "SIRE.GAMETE.PLOIDY", "SIRE.LAMBDA", "DAM.GAMETE.PLOIDY",
                     "DAM.LAMBDA", "PROBABILITY")
  # Add INDIV.PLOIDY column
  ped$INDIV.PLOIDY <- ped$SIRE.GAMETE.PLOIDY + ped$DAM.GAMETE.PLOIDY
  # reorder
  ped <- ped[, c("INDIV.ID", "SIRE.ID", "DAM.ID", "SIRE.GAMETE.PLOIDY", "SIRE.LAMBDA",
                 "DAM.GAMETE.PLOIDY", "DAM.LAMBDA", "PROBABILITY", "INDIV.PLOIDY")]

  #convert NA to 0 in sire and dam columns of pedigree file
  ped[is.na(ped[,"SIRE.ID"]), "SIRE.ID"] <- 0
  ped[is.na(ped[,"DAM.ID"]), "DAM.ID"] <- 0

  #convert ped fields to numeric/integer - intergers not used for identifiers because they are sometimes very large numbers
  ped[,"INDIV.ID"] <- as.numeric(ped[,"INDIV.ID"])
  ped[,"SIRE.ID"] <- as.numeric(ped[,"SIRE.ID"])
  ped[,"DAM.ID"] <- as.numeric(ped[,"DAM.ID"])
  ped[,"SIRE.GAMETE.PLOIDY"] <- as.integer(ped[,"SIRE.GAMETE.PLOIDY"])
  ped[,"SIRE.LAMBDA"] <- as.numeric(ped[,"SIRE.LAMBDA"])
  ped[,"DAM.GAMETE.PLOIDY"] <- as.integer(ped[,"DAM.GAMETE.PLOIDY"])
  ped[,"DAM.LAMBDA"] <- as.numeric(ped[,"DAM.LAMBDA"])
  ped[,"PROBABILITY"] <- as.numeric(ped[,"PROBABILITY"])
  ped[,"INDIV.PLOIDY"] <- as.integer(ped[,"INDIV.PLOIDY"])

  # Check that individuals are in order
  if(sum(!(ped[, "INDIV.ID"] == ped[order(ped[, "INDIV.ID"], decreasing = FALSE), "INDIV.ID"])) > 0) {
    stop("PEDIGREE FILE MUST BE SORTED BY THE INDIVIDUAL IDENTIFIER")
  }

  #Check that SIRE.LAMBDA and DAM.LAMBDA are not less than 0 or greater than 1
  if (max(c(ped[,"SIRE.LAMBDA"], ped[,"DAM.LAMBDA"])) > 1) {
    stop("SIRE AND DAM LAMBDA VALUES IN ped INPUT MUST NOT BE GREATER THAN ONE")
  }

  if (min(c(ped[,"SIRE.LAMBDA"], ped[,"DAM.LAMBDA"])) < 0) {
    stop("SIRE AND DAM LAMBDA VALUES IN ped INPUT MUST NOT BE LESS THAN ZERO")
  }

  # sort by 'INDIV.ID'
  ped <- ped[order(ped[, "INDIV.ID"], decreasing = FALSE), ]

  # CHECK PROBABILITY SUMS TO ONE FOR EACH INDIVIDUAL
  sum.prop <- aggregate(ped[, "PROBABILITY"], by = list(ped[, "INDIV.ID"]), na.rm = T, FUN = "sum")
  colnames(sum.prop) <- c("INDIV.ID", "PROBABILITY.SUM")
  sum.prop[,"PROBABILITY.SUM"] <- round(sum.prop[,"PROBABILITY.SUM"],4)
  if (sum(sum.prop[, "PROBABILITY.SUM"] != 1)) {
    print(sum.prop[sum.prop[, "PROBABILITY.SUM"] != 1, ])
    stop("PROBABILITY DOES NOT SUM TO ONE FOR ALL INDIVIDUALS")
  }
  rm(sum.prop)

  # ENSURE THAT INDIVIDUAL SOMATIC PLOIDY LEVEL IS 2 OR MORE
  if (sum(!rowSums(ped[,c("SIRE.GAMETE.PLOIDY", "DAM.GAMETE.PLOIDY")]) > 1) > 0) {
    stop("THE SOMATIC PLOIDY LEVEL OF INDIVIDUALS MUST BE 2 OR MORE.  polyAinv CANNOT CURRENTLY DEAL WITH MONOPLOID INDIVIDUALS.")
  }

  # ENSURE THAT HALF-SIB INDIVIDUALS ARE NOT PART OF AN AGGREGATE
  half.sibs <- ped[(ped[, "SIRE.ID"] == 0 & ped[, "DAM.ID"] != 0 &
                      ped[, "SIRE.GAMETE.PLOIDY"] != 0) |
                     (ped[, "SIRE.ID"] != 0 & ped[, "DAM.ID"] == 0 &
                        ped[, "DAM.GAMETE.PLOIDY"] != 0), ]
  half.sibs.prop <- ped[(ped[, "SIRE.ID"] == 0 & ped[, "DAM.ID"] != 0 &
                         ped[, "SIRE.GAMETE.PLOIDY"] != 0)  |
                          (ped[, "SIRE.ID"] != 0 & ped[, "DAM.ID"] == 0 &
                             ped[, "DAM.GAMETE.PLOIDY"]  != 0), "PROBABILITY"]
  if (sum(half.sibs.prop != 1) != 0) {
    print(half.sibs)
    stop("HALF-SIB INDIVIDUALS MUST HAVE A PROBABILITY OF 1 (i.e. MUST NOT BE PART OF AN AGGREGATE) IN THE PEDIGREE FILE")
  }
  rm(half.sibs, half.sibs.prop)

  # test if sire gametic ploidy levels equal within individual
  sire.gamete.ploidy.sd <- aggregate(ped[, "SIRE.GAMETE.PLOIDY"], by = list(ped[, "INDIV.ID"]), na.rm=T, FUN = "sd")
  colnames(sire.gamete.ploidy.sd) <- c("INDIV.ID", "SIRE.GAMETE.PLOIDY.SD")
  if(sum(sire.gamete.ploidy.sd[,2] > 0, na.rm = TRUE) != 0) {
    print(sire.gamete.ploidy.sd)
    stop("SIRE GAMETE PLOIDY LEVEL NOT THE SAME FOR ALL POSSIBLE SIRES")
  }
  rm(sire.gamete.ploidy.sd)

  # test if dam gametic ploidy levels equal within individual
  dam.gamete.ploidy.sd <- aggregate(ped[, "DAM.GAMETE.PLOIDY"], by = list(ped[, "INDIV.ID"]), na.rm=T, FUN = "sd")
  colnames(dam.gamete.ploidy.sd) <- c("INDIV.ID", "DAM.GAMETE.PLOIDY.SD")
  if(sum(dam.gamete.ploidy.sd[,2] > 0, na.rm = TRUE) != 0) {
    print(dam.gamete.ploidy.sd)
    stop("DAM GAMETE PLOIDY LEVEL NOT THE SAME FOR ALL POSSIBLE DAMS")
  }
  rm(dam.gamete.ploidy.sd)

  # test if sire lambda levels equal within individual
  sire.lambda.sd <- aggregate(ped[, "SIRE.LAMBDA"], by = list(ped[, "INDIV.ID"]), na.rm=T, FUN = "sd")
  colnames(sire.lambda.sd) <- c("INDIV.ID", "SIRE.LAMBDA.SD")
  if(sum(sire.lambda.sd[,2] > 0, na.rm = TRUE) != 0) {
    print(sire.lambda.sd)
    stop("SIRE LAMBDA NOT THE SAME FOR ALL POSSIBLE SIRES")
  }
  rm(sire.lambda.sd)

  # test if dam lambda levels equal within individual
  dam.lambda.sd <- aggregate(ped[, "DAM.LAMBDA"], by = list(ped[, "INDIV.ID"]), na.rm=T, FUN = "sd")
  colnames(dam.lambda.sd) <- c("INDIV.ID", "DAM.LAMBDA.SD")
  if(sum(dam.lambda.sd[,2] > 0, na.rm = TRUE) != 0) {
    print(dam.lambda.sd)
    stop("DAM LAMBDA NOT THE SAME FOR ALL POSSIBLE DAMS")
  }
  rm(dam.lambda.sd)

  #########################################################################################
  # Generate founder.lambda
  #########################################################################################

  founders.dam   <- unique(ped[ped[,"DAM.ID"] == 0, c("DAM.GAMETE.PLOIDY", "DAM.LAMBDA")])
  founders.dam$INDIV.PLOIDY   <- founders.dam$DAM.GAMETE.PLOIDY * 2
  founders.dam  <- founders.dam[,c("INDIV.PLOIDY", "DAM.LAMBDA")]

  founders.sire   <- unique(ped[ped[,"SIRE.ID"] == 0, c("SIRE.GAMETE.PLOIDY", "SIRE.LAMBDA")])
  founders.sire$INDIV.PLOIDY   <- founders.sire$SIRE.GAMETE.PLOIDY * 2
  founders.sire  <- founders.sire[,c("INDIV.PLOIDY", "SIRE.LAMBDA")]

  founder.lambda <- merge(founders.sire, founders.dam, by = "INDIV.PLOIDY", all = TRUE)
  rm(founders.dam, founders.sire)

  #convert founder.lambda fields to numeric or integer
  colnames(founder.lambda) <- c("INDIV.PLOIDY", "SIRE.LAMBDA", "DAM.LAMBDA")

  founder.lambda[,"INDIV.PLOIDY"] <- as.integer(founder.lambda[,"INDIV.PLOIDY"])
  founder.lambda[,"SIRE.LAMBDA"]  <- as.numeric(founder.lambda[,"SIRE.LAMBDA"])
  founder.lambda[,"DAM.LAMBDA"]   <- as.numeric(founder.lambda[,"DAM.LAMBDA"])

  founder.lambda <- founder.lambda[founder.lambda[,"INDIV.PLOIDY"] != 0,]
  founder.lambda[is.na(founder.lambda[,"SIRE.LAMBDA"]),"SIRE.LAMBDA"] <- 0
  founder.lambda[is.na(founder.lambda[,"DAM.LAMBDA"]),"DAM.LAMBDA"]   <- 0

  #Check that SIRE.LAMBDA and DAM.LAMBDA are not less than 0 or greater than 1
  if (max(c(founder.lambda[,"SIRE.LAMBDA"], founder.lambda[,"DAM.LAMBDA"])) > 1) {
    stop("FOUNDER SIRE AND DAM LAMBDA VALUES MUST NOT BE GREATER THAN ONE")
  }

  if (min(c(founder.lambda[,"SIRE.LAMBDA"], founder.lambda[,"DAM.LAMBDA"])) < 0) {
    stop("FOUNDER SIRE AND DAM LAMBDA VALUES MUST NOT BE LESS THAN ZERO")
  }

  count.by.ploidy <- aggregate(!is.na(founder.lambda[, "INDIV.PLOIDY"]),
                               by = list(founder.lambda[, "INDIV.PLOIDY"]),
                               na.rm=T, FUN = "sum")

  if(sum(count.by.ploidy[,2] != 1) > 0) {
    stop("LAMBDA VALUES ARE NOT THE SAME WITHIN SEX AND PLOIDY LEVEL FOR FOUNDERS")
  }
  rm(count.by.ploidy)
  if (sum(founder.lambda[,"INDIV.PLOIDY"] == 2) > 0) {
    if (founder.lambda[founder.lambda[,"INDIV.PLOIDY"] == 2 ,"SIRE.LAMBDA"] != 0 |
        founder.lambda[founder.lambda[,"INDIV.PLOIDY"] == 2 ,"DAM.LAMBDA"]  != 0) {
      stop("DIPLOID FOUNDERS MUST HAVE A SIRE.LAMBDA AND DAM.LAMBDA OF ZERO")
    }
  }
  #########################################################################################
  # Reorder pedigree so that parents are at the top
  #########################################################################################

  print("Rearranging pedigree")

  ped.parents     <- ped[ped[,"INDIV.ID"] %in% parent.ids,]
  ped.parents$PARENT <- TRUE
  ped.non.parents <- ped[!ped[,"INDIV.ID"] %in% parent.ids,]
  ped.non.parents$PARENT <- FALSE
  ped <- rbind(ped.parents, ped.non.parents)
  rm(ped.parents, ped.non.parents, parent.ids)

  #Recode individuals in pedigree file to be sequential numbers from one
  col.order <- colnames(ped)
  col.order <- c(col.order,"INDIV.ID.ORIG", "SIRE.ID.ORIG", "DAM.ID.ORIG")
  colnames(ped)[1:3] <- c("INDIV.ID.ORIG", "SIRE.ID.ORIG", "DAM.ID.ORIG")

  tmp <- data.frame(DAM.ID.ORIG = unique(ped[,"INDIV.ID.ORIG"]),
                    DAM.ID      = 1:length(unique(ped[,"INDIV.ID.ORIG"])))
  tmp <- rbind(tmp,c(0,0))
  ped <- merge(ped, tmp, by = "DAM.ID.ORIG", all.x = TRUE)
  colnames(tmp) <- c("SIRE.ID.ORIG", "SIRE.ID")
  ped <- merge(ped, tmp, by = "SIRE.ID.ORIG", all.x = TRUE)
  colnames(tmp) <- c("INDIV.ID.ORIG", "INDIV.ID")
  ped <- merge(ped, tmp, by = "INDIV.ID.ORIG", all.x = TRUE)
  ped <- ped[,col.order]
  rm(tmp, col.order)

  #convert ped to matrix
  ped <- as.matrix(ped)
  ped <- ped[order(ped[, "INDIV.ID"], decreasing = FALSE),]

  #Non-parents ###################################

  fam.ped <- ped[ped[,"PARENT"] == 0,!colnames(ped) %in% c("INDIV.ID.ORIG", "SIRE.ID.ORIG", "DAM.ID.ORIG")]

  #Identify families
  if(is.matrix(fam.ped)) {
  fam.ped2 <- data.frame(INDIV.ID = fam.ped[,1],
                         UNIQUE.STRING = paste(fam.ped[,2], fam.ped[,3], fam.ped[,4], fam.ped[,5], fam.ped[,6],
                                               fam.ped[,7], fam.ped[,8], fam.ped[,9], fam.ped[,10], sep = "."))
  } else {
    stop("polyAinv REQUIRES AT LEAST TWO NON-PARENTS TO BE PRESENT IN THE PEDIGREE FILE.  ADD DUMMY INDIVIDUALS IF NECESSARY")
  }
  fam.ped2$UNIQUE.STRING <- as.character(fam.ped2$UNIQUE.STRING)
  ped.unique.rows <- unique(fam.ped2[,"UNIQUE.STRING"])
  ped.unique.rows <- data.frame(UNIQUE.STRING = ped.unique.rows,
                                UNIQUE.ROW.ID = 1:length(ped.unique.rows)) #Unique rows
  ped.unique.rows$UNIQUE.STRING <- as.character(ped.unique.rows$UNIQUE.STRING)

  fam.ped2 <- merge(fam.ped2, ped.unique.rows, by = "UNIQUE.STRING", all.x=TRUE) #Merge unique rows

  fam.ped2 <- fam.ped2[order(fam.ped2[,"UNIQUE.ROW.ID"], decreasing = FALSE),]
  fam.ped2 <- fam.ped2[order(fam.ped2[,"INDIV.ID"], decreasing = FALSE),]

  fam.ped2 <- aggregate(UNIQUE.ROW.ID ~ INDIV.ID, data = fam.ped2, c)
  fam.ped2$UNIQUE.ROW.ID <- paste(fam.ped2$UNIQUE.ROW.ID)

  ped.unique.rows <- unique(fam.ped2$UNIQUE.ROW.ID) #Idenify unique rows in pivot table (i.e. families)
  ped.unique.rows <- data.frame(FAM.ID = fam.ped2[1,"INDIV.ID"]:(fam.ped2[1,"INDIV.ID"] + length(ped.unique.rows) - 1),
                                UNIQUE.ROW.ID = ped.unique.rows)

  indiv.fam.lookup <- merge(fam.ped2, ped.unique.rows, by = "UNIQUE.ROW.ID", all.x=TRUE, sort = TRUE) #Merge so that individual and family identifiers in one matrix
  indiv.fam.lookup <- indiv.fam.lookup[,c("INDIV.ID", "FAM.ID")] #retain only individual identifiers and family identifiers
  indiv.fam.lookup <- indiv.fam.lookup[order(indiv.fam.lookup[,"INDIV.ID"], decreasing = FALSE),]

  rm(fam.ped2, ped.unique.rows)

  #Get count of individuals by family
  n.by.family <- aggregate(!is.na(indiv.fam.lookup$FAM.ID), by = list(indiv.fam.lookup$FAM.ID), na.rm=T, FUN = "sum")
  colnames(n.by.family) <- c("FAM.ID", "COUNT")

  #remove duplicated families from fam.ped
  unique.fams <- indiv.fam.lookup[!duplicated(indiv.fam.lookup[,"FAM.ID"]),]
  fam.ped <- merge(fam.ped, unique.fams, by = "INDIV.ID", all.y = TRUE)
  fam.ped$INDIV.ID <- fam.ped$FAM.ID
  fam.ped <- fam.ped[,1:10]
  rm(unique.fams)

  #Add parents
  fam.ped <- rbind(ped[ped[,"PARENT"] == 1,!colnames(ped) %in% c("INDIV.ID.ORIG", "SIRE.ID.ORIG", "DAM.ID.ORIG")], fam.ped)
  fam.ped <- merge(fam.ped, n.by.family, by = 1 , all.x = TRUE)
  fam.ped[is.na(fam.ped[,"COUNT"]),"COUNT"] <- 1

  fam.ped <- fam.ped[order(fam.ped[,"INDIV.ID"], decreasing = FALSE),
                     c("INDIV.ID", "SIRE.ID", "DAM.ID", "SIRE.GAMETE.PLOIDY", "SIRE.LAMBDA",
                       "DAM.GAMETE.PLOIDY", "DAM.LAMBDA",  "PROBABILITY", "INDIV.PLOIDY",  "PARENT", "COUNT")]

  #########################################################################################
  #Generate K and K.inverse
  #########################################################################################

  print("Generating K and K inverse matrices")
  print(paste("0% complete.",Sys.time()))

  parent.ids <- unique(fam.ped[fam.ped[,"PARENT"] == 1,"INDIV.ID"])

  indiv.ids =  unique(fam.ped[,"INDIV.ID"])

  # Define matrices
  K.parents <- matrix(0, ncol = length(parent.ids), nrow = length(parent.ids))

  K.inv.mat <- data.frame(FAM.ID = NA,
                          INDIV.1 = NA,
                          INDIV.2 = NA,
                          K.INV = NA)[-1,]

  K.inv.vec <- data.frame(INDIV.1 = NA,
                          INDIV.2 = NA,
                          K.INV = NA)[-1,]

  K.diag <- matrix(0,ncol=1, max(indiv.ids))

  # Loop by each INDIV.ID in fam.ped and generate K and K.inv ###################################
  for (indiv.count in indiv.ids) {

    indiv.parents           <- fam.ped[fam.ped[, "INDIV.ID"] == indiv.count, ]
    if(is.null(nrow(indiv.parents))) {
      indiv.parents           <- matrix(indiv.parents, ncol = 11) #one row
    }
    colnames(indiv.parents) <- colnames(fam.ped)

    indiv.sire <- indiv.parents
    colnames(indiv.sire)[colnames(indiv.sire) == "SIRE.ID"] <- "PARENT.ID"
    colnames(indiv.sire)[colnames(indiv.sire) == "DAM.ID"] <- "OTHER.PARENT.ID"
    colnames(indiv.sire)[colnames(indiv.sire) == "SIRE.GAMETE.PLOIDY"] <- "PARENT.GAMETE.PLOIDY"
    colnames(indiv.sire)[colnames(indiv.sire) == "DAM.GAMETE.PLOIDY"] <- "OTHER.PARENT.GAMETE.PLOIDY"
    colnames(indiv.sire)[colnames(indiv.sire) == "SIRE.LAMBDA"] <- "PARENT.LAMBDA"
    tmp <- get.p.and.q(ped = fam.ped, indiv.parents = indiv.sire, K.parents = K.parents, indiv.count = indiv.count)
    rm(indiv.sire)
    sire.pr.ibd <- tmp$parent.pr.ibd
    p <- tmp$p.q
    colnames(p) <-  c("INDIV.ID", "p")
    p.multiplier <- tmp$p.q.multiplier
    rm(tmp)

    indiv.dam <- indiv.parents
    colnames(indiv.dam)[colnames(indiv.dam) == "DAM.ID"] <- "PARENT.ID"
    colnames(indiv.dam)[colnames(indiv.dam) == "SIRE.ID"] <- "OTHER.PARENT.ID"
    colnames(indiv.dam)[colnames(indiv.dam) == "DAM.GAMETE.PLOIDY"] <- "PARENT.GAMETE.PLOIDY"
    colnames(indiv.dam)[colnames(indiv.dam) == "SIRE.GAMETE.PLOIDY"] <- "OTHER.PARENT.GAMETE.PLOIDY"
    colnames(indiv.dam)[colnames(indiv.dam) == "DAM.LAMBDA"] <- "PARENT.LAMBDA"
    tmp <- get.p.and.q(ped = fam.ped, indiv.parents = indiv.dam, K.parents = K.parents, indiv.count = indiv.count)
    rm(indiv.dam)
    dam.pr.ibd <- tmp$parent.pr.ibd
    q <- tmp$p.q
    colnames(q) <-  c("INDIV.ID", "q")
    q.multiplier <- tmp$p.q.multiplier
    rm(tmp)

    # Compute s and sire.dam.pr.ibd ###################################
    # Equation 17 of Hamilton and Kerr 2018
    if(p[1, "INDIV.ID"] == 0 | q[1, "INDIV.ID"] == 0) { #can't have half-sibs in aggreagates
      sire.dam.pr.ibd <- 0
    } else {
      sire.dam.pr.ibd <- t(as.matrix(p[, "p"])) %*% K.parents[p[, "INDIV.ID"], q[, "INDIV.ID"]] %*% as.matrix(q[,"q"])
    }

    s <- merge(p, q, by = "INDIV.ID", all = TRUE)
    s[is.na(s)] <- 0
    s[, "s"] <- s[, "p"] * p.multiplier + s[, "q"] * q.multiplier
    s <- s[, c("INDIV.ID", "s")]

    s <- s[s[,"INDIV.ID"] != 0,] #remove founders

    # Founders in K ###################################

    if ((indiv.parents[1, "SIRE.ID"] == 0 & indiv.parents[1, "DAM.ID"] == 0)) { # If founder in K

      # Equation 12 of Hamilton and Kerr 2018
      Fi <- get.founder.F(indiv.ploidy = indiv.parents[1, "INDIV.PLOIDY"],
                          sire.lambda = indiv.parents[1, "SIRE.LAMBDA"],
                          dam.lambda = indiv.parents[1, "DAM.LAMBDA"])

    } else {

      # Equation 11 of Hamilton and Kerr 2018
      Fi <- (indiv.parents[1, "SIRE.GAMETE.PLOIDY"] *
               (indiv.parents[1, "SIRE.GAMETE.PLOIDY"] - 1) * sire.pr.ibd +
               indiv.parents[1, "DAM.GAMETE.PLOIDY"] *
               (indiv.parents[1, "DAM.GAMETE.PLOIDY"] - 1) * dam.pr.ibd +
               2 * indiv.parents[1, "SIRE.GAMETE.PLOIDY"] *
               indiv.parents[1, "DAM.GAMETE.PLOIDY"] * sire.dam.pr.ibd)/
        (indiv.parents[1, "SIRE.GAMETE.PLOIDY"] * (indiv.parents[1, "SIRE.GAMETE.PLOIDY"] - 1) +
           indiv.parents[1, "DAM.GAMETE.PLOIDY"] *
           (indiv.parents[1, "DAM.GAMETE.PLOIDY"] - 1) + 2 *
           indiv.parents[1, "SIRE.GAMETE.PLOIDY"] *
           indiv.parents[1,"DAM.GAMETE.PLOIDY"])

    }  #END else if((indiv.parents[1,'SIRE.ID'] == 0 & indiv.parents[1,'DAM.ID'] == 0)) {   #If founder in K
    rm(sire.pr.ibd, dam.pr.ibd, sire.dam.pr.ibd)

    # Equation 5 of Hamilton and Kerr 2018
    kii <- (1 + (indiv.parents[1, "INDIV.PLOIDY"] - 1) * Fi)/indiv.parents[1, "INDIV.PLOIDY"]
    K.diag[indiv.count-min(indiv.ids)+1] <- kii

    # Generate F vector
    if (indiv.count == 1) {
      F.vector <- data.frame(INDIV.ID = indiv.parents[1, "INDIV.ID"], F = as.matrix(Fi))
    } else {
      F.vector <- rbind(F.vector, data.frame(INDIV.ID = indiv.parents[1, "INDIV.ID"], F = as.matrix(Fi)))
    }

    if(indiv.count <= length(parent.ids)) {

      # Add diag to K.parents
      K.parents[indiv.parents[1, "INDIV.ID"], indiv.parents[1, "INDIV.ID"]] <- kii

      # Off-diagonal elements of non-founders in K ###################################

      # Coancestry between INDIV.IDs and parents
      if (indiv.count > 1) {
        tmp.s <- merge(c(1:(indiv.count-1)), s, by = 1, all.x = TRUE)
        tmp.s[is.na(tmp.s)] <- 0
        s.by.K.parents <- K.parents[1:(indiv.count-1), 1:(indiv.count-1)] %*% tmp.s[,"s"]
        rm(tmp.s)
        K.parents[1:(indiv.count-1), indiv.count] <- s.by.K.parents
        K.parents[indiv.count, 1:(indiv.count-1)] <- t(s.by.K.parents)
        rm(s.by.K.parents)
      }
    }

    # K.INV ###################################

    # Equation 14 of Hamilton and Kerr 2018
    if(nrow(s) == 0) {
      rhs.left <- 1/as.numeric(kii)
    } else {
      rhs.left <- 1/as.numeric(kii - t(s[, "s"]) %*% K.parents[s[, "INDIV.ID"], s[, "INDIV.ID"]] %*% s[, "s"])
    }

    rhs.right.mat <- matrix(0, nrow = nrow(s), ncol = nrow(s))
    rhs.right.vec <- matrix(0, nrow = 1, ncol = nrow(s)+1)

    if(indiv.count > 1) {
      rhs.right.mat <- s[, "s"] %*% t(s[, "s"])
      rhs.right.mat <- rhs.right.mat[lower.tri(rhs.right.mat, diag = TRUE)]
      rhs.mat <- rhs.left * rhs.right.mat

      indiv.1 <- matrix(rep(s[,"INDIV.ID"], nrow(s)), ncol = nrow(s))
      indiv.1 <- indiv.1[lower.tri(indiv.1, diag = TRUE)]
      indiv.2 <- matrix(rep(s[,"INDIV.ID"],each = nrow(s)), ncol = nrow(s))
      indiv.2 <- indiv.2[lower.tri(indiv.2, diag = TRUE)]
      rhs.mat <- data.frame(FAM.ID = rep(indiv.count,length(indiv.1)),
                            INDIV.1 = indiv.1,
                            INDIV.2 = indiv.2,
                            K.INV = rhs.mat)
      rhs.mat <- rhs.mat[rhs.mat[,"K.INV"] != 0,] #remove zero elements
      K.inv.mat <- rbind(K.inv.mat, rhs.mat)
    }

    rhs.right.vec <- c(-s[, "s"],1)
    rhs.vec <- rhs.left * rhs.right.vec
    rhs.vec <- data.frame(INDIV.1 = rep(indiv.count,length(s[,"INDIV.ID"])+1),
                          INDIV.2 = c(s[,"INDIV.ID"], indiv.count),
                          K.INV = as.vector(rhs.vec))
    rhs.vec <- rhs.vec[rhs.vec[,"K.INV"] != 0,] #remove zero elements
    K.inv.vec <- rbind(K.inv.vec, rhs.vec)

    #Print progress update - 10% increments
    if(indiv.count %in% seq(0,round(max(indiv.ids),-1)-round(max(indiv.ids),-1)/10,round(max(indiv.ids),-1)/10)) {
      print(paste(indiv.count / round(max(indiv.ids),-1) * 100 , "% complete. ", Sys.time(), sep = ""))
    }

  }  #END Loop by each INDIV.ID in fam.ped and generate K and K.inv

  #remove unnecessary values
  rm(fam.ped, indiv.1, indiv.2, indiv.count, indiv.ids, p.multiplier, parent.ids, q.multiplier,
     rhs.left, rhs.right.vec, rhs.right.mat)

  #remove unnecessary data
  rm(Fi, indiv.parents, kii, p, q, rhs.mat, rhs.vec, s)

  if(ASReml.giv.only) {
    rm(F.vector, K.diag, K.parents)
  }

  print(paste("100% complete.", Sys.time()))

  #########################################################################################
  #Rearrange matrices
  #########################################################################################

  print("Rearranging matrices")

  #########################################################################################
  #Duplicate according to the number of individuals represented by a family
  #########################################################################################
  K.inv.mat <- merge(K.inv.mat, n.by.family, by = "FAM.ID", all.x = TRUE)
  K.inv.mat[is.na(K.inv.mat[,"COUNT"]),"COUNT"] <- 1
  K.inv.mat[,"K.INV"] <- K.inv.mat[,"K.INV"] * K.inv.mat[,"COUNT"]
  K.inv.mat <- K.inv.mat[!is.na(K.inv.mat[,"K.INV"]),c("INDIV.1", "INDIV.2", "K.INV")]

  #########################################################################################
  #Expand K.inv.vec to individuals
  #########################################################################################

  #Add parents to indiv.fam.lookup
  indiv.fam.lookup <- rbind(data.frame(INDIV.ID = 1:(min(indiv.fam.lookup$INDIV.ID)-1),
                                       FAM.ID= 1:(min(indiv.fam.lookup$FAM.ID)-1)),
                            indiv.fam.lookup)

  K.inv.vec <- merge(indiv.fam.lookup, K.inv.vec, by.x = "FAM.ID", by.y = "INDIV.2", all.x = TRUE)
  colnames(K.inv.vec)[1] <- "FAM.2"
  K.inv.vec <- merge(indiv.fam.lookup, K.inv.vec, by.x = "FAM.ID", by.y = "INDIV.1", all.x = TRUE)
  colnames(K.inv.vec)[1] <- "FAM.1"
  colnames(K.inv.vec) <- c("FAM.1", "INDIV.1", "FAM.2", "INDIV.2", "K.INV")
  #remove rows where FAM.1 == FAM.2 but INDIV.1 != INDIV.2
  K.inv.vec <- K.inv.vec[!(K.inv.vec[,"FAM.1"] == K.inv.vec[,"FAM.2"] &
                             K.inv.vec[,"INDIV.1"] != K.inv.vec[,"INDIV.2"]),  ]
  K.inv.vec <- K.inv.vec[,c("INDIV.1", "INDIV.2", "K.INV")]

  #########################################################################################
  # Merge, aggregate and sort K.inv.mat, K.inv.vec
  #########################################################################################

  K.inv.matrix <- rbind(K.inv.mat, K.inv.vec)
  rm(K.inv.mat, K.inv.vec)
  K.inv.matrix <- aggregate(K.inv.matrix$K.INV, by=list(K.inv.matrix$INDIV.1,K.inv.matrix$INDIV.2), FUN = sum)
  colnames(K.inv.matrix) <- c("INDIV.1", "INDIV.2", "K.INV")
  K.inv.matrix <- K.inv.matrix[order(K.inv.matrix[,"INDIV.2"]),]
  K.inv.matrix <- K.inv.matrix[order(K.inv.matrix[,"INDIV.1"]),]

  #########################################################################################
  # Generate A.inv matrix
  #########################################################################################

  print("Generating A inverse matrix")

  # Equation 4 of Hamilton and Kerr 2018
  ploidy <- unique(ped[, c("INDIV.ID", "INDIV.PLOIDY")])
  A.inv.matrix <- K.inv.matrix
  if(ASReml.giv.only) {rm(K.inv.matrix)}
  colnames(A.inv.matrix) <- c("INDIV.1", "INDIV.2", "A.INV")
  A.inv.matrix <- merge(A.inv.matrix, ploidy, by.x = "INDIV.1", by.y = "INDIV.ID", all.x = TRUE)
  A.inv.matrix <- merge(A.inv.matrix, ploidy, by.x = "INDIV.2", by.y = "INDIV.ID", all.x = TRUE)

  A.inv.matrix[, "A.INV"] <- A.inv.matrix[, "A.INV"]/(2 * sqrt(A.inv.matrix[, "INDIV.PLOIDY.x"]/2 *
                                                                 A.inv.matrix[, "INDIV.PLOIDY.y"]/2))
  A.inv.matrix            <- A.inv.matrix[, c("INDIV.1", "INDIV.2", "A.INV")]

  ###################################################################################
  #ASReml .giv
  ###################################################################################

  ASReml.giv <- A.inv.matrix
  if(ASReml.giv.only) {rm(A.inv.matrix)}

  indiv.index <- unique(ped[, "INDIV.ID.ORIG"])
  indiv.index <- indiv.index[order(indiv.index, decreasing = FALSE)]
  indiv.index <- data.frame(INDIV.ID.ORIG = indiv.index,
                            INDIV.ID.INDEX = 1:length(indiv.index))
  ped         <- merge(ped, indiv.index, by = "INDIV.ID.ORIG", all.x = TRUE)

  ASReml.giv <- merge(ASReml.giv,
                        unique(ped[,c("INDIV.ID", "INDIV.ID.INDEX")]),
                        by.x = "INDIV.2",
                        by.y = "INDIV.ID",
                        all.x = TRUE)

  ASReml.giv <- merge(ASReml.giv,
                        unique(ped[,c("INDIV.ID", "INDIV.ID.INDEX")]),
                        by.x = "INDIV.1",
                        by.y = "INDIV.ID",
                        all.x = TRUE)
  colnames(ASReml.giv) <- c("INDIV.1.REMOVE", "INDIV.2.REMOVE", "Ainverse", "Column", "Row")
  ASReml.giv <- ASReml.giv[,c("Row", "Column", "Ainverse")]

  ASReml.giv <- ASReml.giv[order(ASReml.giv[, "Column"], decreasing = FALSE),]
  ASReml.giv <- ASReml.giv[order(ASReml.giv[, "Row"], decreasing = FALSE),]

  #ASReml-R ginv object attributes (see 'asreml.Ainverse' 'in Butler et al. 2009)

  rownames(ASReml.giv) <- 1:nrow(ASReml.giv)
  attr(ASReml.giv, "rowNames") <- as.character(indiv.index[,"INDIV.ID.ORIG"])
  attr(ASReml.giv, "geneticGroups") <- c(0,0)

  #########################################################################################
  #Generate A matrix
  #########################################################################################

  if(!ASReml.giv.only) {

    print("Generating A matrix")

    # A.matrix
    parent.ids <- unique(c(ped[, "SIRE.ID"], ped[, "DAM.ID"]))
    parent.ids <- parent.ids[order(parent.ids, decreasing = FALSE)]  #sort
    parent.ids <- parent.ids[parent.ids != 0]  #remove 0

    lhs <- 2 * sqrt(ploidy[parent.ids, "INDIV.PLOIDY"]/2) %*% t(sqrt(ploidy[parent.ids, "INDIV.PLOIDY"]/2))

    A.parents <- lhs * K.parents

    #########################################################################################
    #Generate diagonals
    #########################################################################################

    diagonals <- cbind(F.vector,K.diag)
    rm(F.vector, K.diag)

    #Expand diagonals to individuals
    diagonals <- merge(indiv.fam.lookup, diagonals, by.x = "FAM.ID", by.y = "INDIV.ID", all.x = TRUE)
    diagonals <- diagonals[order(diagonals[,"INDIV.ID"]),]

    #A.DIAG and PLOIDY
    diagonals$A.DIAG <- ploidy[, "INDIV.PLOIDY"] * diagonals$K.diag
    diagonals$PLOIDY <- ploidy[, "INDIV.PLOIDY"]
    colnames(diagonals)[colnames(diagonals) == "K.diag"] <- "K.DIAG"
    diagonals <- diagonals[,c("INDIV.ID", "PLOIDY", "F", "K.DIAG", "A.DIAG" )]

    #########################################################################################
    #Revert to original ids
    #########################################################################################

    print("Reverting to original identifiers")

    parent.ids.orig <-  (unique(ped[ped[,"PARENT"] == 1,"INDIV.ID.ORIG"]))
    colnames(K.parents) <- parent.ids.orig
    rownames(K.parents) <- parent.ids.orig

    colnames(A.parents) <- parent.ids.orig
    rownames(A.parents) <- parent.ids.orig

    K.inv.matrix <- merge(K.inv.matrix,
                          unique(ped[,c("INDIV.ID", "INDIV.ID.ORIG")]),
                          by.x = "INDIV.2",
                          by.y = "INDIV.ID",
                          all.x = TRUE)

    K.inv.matrix <- merge(K.inv.matrix,
                          unique(ped[,c("INDIV.ID", "INDIV.ID.ORIG")]),
                          by.x = "INDIV.1",
                          by.y = "INDIV.ID",
                          all.x = TRUE)
    colnames(K.inv.matrix) <- c("INDIV.1.REMOVE", "INDIV.2.REMOVE", "K.INV", "INDIV.2", "INDIV.1")
    K.inv.matrix <- K.inv.matrix[,c("INDIV.1", "INDIV.2", "K.INV")]

    A.inv.matrix <- merge(A.inv.matrix,
                          unique(ped[,c("INDIV.ID", "INDIV.ID.ORIG")]),
                          by.x = "INDIV.2",
                          by.y = "INDIV.ID",
                          all.x = TRUE)

    A.inv.matrix <- merge(A.inv.matrix,
                          unique(ped[,c("INDIV.ID", "INDIV.ID.ORIG")]),
                          by.x = "INDIV.1",
                          by.y = "INDIV.ID",
                          all.x = TRUE)
    colnames(A.inv.matrix) <- c("INDIV.1.REMOVE", "INDIV.2.REMOVE", "A.INV", "INDIV.2", "INDIV.1")
    A.inv.matrix <- A.inv.matrix[,c("INDIV.1", "INDIV.2", "A.INV")]

    diagonals[order(diagonals[,"INDIV.ID"]),"INDIV.ID"] <- unique(ped[order(ped[,"INDIV.ID"]),"INDIV.ID.ORIG"])

    #reorder
    K.parents <- K.parents[order(as.numeric(rownames(K.parents)), decreasing = FALSE),
                           order(as.numeric(colnames(K.parents)), decreasing = FALSE)]

    A.parents <- A.parents[order(as.numeric(rownames(A.parents)), decreasing = FALSE),
                           order(as.numeric(colnames(A.parents)), decreasing = FALSE)]

    K.inv.matrix <- K.inv.matrix[order(K.inv.matrix[, "INDIV.2"], decreasing = FALSE),]
    K.inv.matrix <- K.inv.matrix[order(K.inv.matrix[, "INDIV.1"], decreasing = FALSE),]

    A.inv.matrix <- A.inv.matrix[order(A.inv.matrix[, "INDIV.2"], decreasing = FALSE),]
    A.inv.matrix <- A.inv.matrix[order(A.inv.matrix[, "INDIV.1"], decreasing = FALSE),]

    diagonals <- diagonals[order(diagonals[, "INDIV.ID"], decreasing = FALSE),]

    #revert to orginal column names
    colnames(K.inv.matrix) <- c(paste(orig.col.names[1], "1", sep = ""),
                                paste(orig.col.names[1], "2", sep = ""),
                                "K.INV")

    colnames(A.inv.matrix) <- c(paste(orig.col.names[1], "1", sep = ""),
                                paste(orig.col.names[1], "2", sep = ""),
                                "A.INV")

    return(list(K.parents  = K.parents,
                K.inv  = K.inv.matrix,
                A.parents  = A.parents,
                A.inv  = A.inv.matrix,
                ASReml.giv = ASReml.giv,
                F = diagonals,
                founder.lambda = founder.lambda))

  } else {
    return(ASReml.giv)
  } #END   if(!ASReml.giv.only) {

  print(Sys.time())


}  #END polyAinv
