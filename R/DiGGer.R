# ------------------------------------------------------------------------------------------
# Description: These function was lifted from the DiGGer package version 1.0.5 which 
#              will be used to generate randomization for different designs.
# Source: DiGGer 1.0.5
# ------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------
#' @name prDiGGer 
#' @aliases prDiGGer 
#' @title prDiGGer 
#' 
#' @description Generate a DiGGer search for designs with some unreplicated treatments
#' @param numberOfTreatments Number of distinct treatments in the design.
#' @param rowsInDesign Number of rows in the design.
#' @param columnsInDesign Number of columns in the design.
#' @param blockSequence List of dimension pairs of blocks to be optimised in sequence.
#' @param betweenRowCorrModel The correlation pattern between rows may be "AR" AutoRegressive, "MA" Moving Average, or "ID" no correlation.
#' @param betweenRowCorr Valid parameter values p, for "AR" -1<p<1, for "MA" -0.5<p<0.5 and for "ID" p=0.
#' @param betweenColumnCorrModel As for betweenRowCorrModel.
#' @param betweenColumnCorr As for betweenRowCorr.
#' @param treatName Vector of texts giving the names of the treatment.
#' @param treatNumber Vector of numbers associated with the treatments.
#' @param treatRepPerRep Vector of replication levels for each treatment over the whole design.
#' @param treatGroup Vector of group codes (up to 200 distinct values) associated with treatments. Group codes may be used to modify the A-efficiency measure in the optimisation.
#' @param maxInterchanges Maximum number of interchanges used in the search to spatially optimise the block design created by the call to prDiGGer.
#' @param targetGroup Test treatment group used for optimisation in the spatial search phase.
#' @param runSearch Logical value, whether to run the final spatial search immediately.
#' @param splitOpt Dimension pair used for sub-optimal spatial distribution of replicated entries.
#' @param rngSeeds Seeds c(s1,s2) to control the DiGGer search. s1 must be in the range [0,31328], s2 must be in the range [0,30081].
#' @source DiGGer 1.0.5
# ------------------------------------------------------------------------------------------

prDiGGer <- function (numberOfTreatments, rowsInDesign, columnsInDesign, 
                      blockSequence = NULL, betweenRowCorrModel = "AR", betweenRowCorr = 0.5, 
                      betweenColumnCorrModel = "AR", betweenColumnCorr = 0.5, 
                      treatName = NULL, treatNumber = NULL, treatRepPerRep = NULL, 
                      treatGroup = NULL, maxInterchanges = 20000, targetGroup = 1, 
                      runSearch = FALSE, splitOpt = NULL, rngSeeds = NULL) 
{
  if (is.na(rowsInDesign)) 
    stop("rowsInDesign must be given")
  if (is.na(columnsInDesign)) 
    stop("columnsInDesign must be given")
  if (is.na(numberOfTreatments)) 
    stop("numberOfTreatments must be given")
  if (is.null(treatRepPerRep)) 
    stop("treatRepPerRep must be given")
  if (sum(treatRepPerRep) != rowsInDesign * columnsInDesign) 
    stop("Treatment repeats do not equal number of plots in design")
  if (is.null(blockSequence)) {
    stop("Blocking sequence must be specified for prDiGGer")
  }
  else {
    mbs <- matrix(unlist(blockSequence), ncol = 2, byrow = TRUE)
    rowsInBlockSequence = mbs[, 1]
    columnsInBlockSequence = mbs[, 2]
  }
  if (is.null(rngSeeds)) {
    rngSeeds <- c(sample(seq(31329), 1) - 1, sample(seq(30082), 
                                                    1) - 1)
  }
  else {
    if (length(rngSeeds) == 2) {
      if (rngSeeds[1] == -1) {
        rngSeeds[1] <- sample(seq(31329), 1) - 1
      }
      else {
        if (rngSeeds[1] < 0 | rngSeeds[1] > 31328) 
          stop("rngSeeds[1] must be 0 <= seed 1 <= 31328")
      }
      if (rngSeeds[2] == -1) {
        rngSeeds[2] <- sample(seq(30082), 1) - 1
      }
      else {
        if (rngSeeds[2] < 0 | rngSeeds[2] > 30081) 
          stop("rngSeeds[2] must be 0 <= seed 2 <= 30081")
      }
    }
    else stop("error in specifying rngSeeds")
  }
  rng <- rmarin(rngSeeds[1], rngSeeds[2])
  nrp <- sum(treatRepPerRep[treatRepPerRep > 1])
  nrtb <- rev(rowsInBlockSequence)[1]
  nctb <- rev(columnsInBlockSequence)[1]
  ntrb <- ceiling(rowsInDesign/nrtb)
  ntcb <- ceiling(columnsInDesign/nctb)
  ntb <- ntrb * ntcb
  reptb <- ceiling(nrp/ntb)
  nrrd <- reptb * ntrb
  ncrd <- ntcb
  ndummy <- nrrd * ncrd - nrp
  nttr <- treatRepPerRep[treatRepPerRep > 1 & treatGroup == 
                           targetGroup]
  ntt <- length(nttr)
  nttp <- sum(nttr)
  nctr <- treatRepPerRep[treatRepPerRep > 1 & treatGroup != 
                           targetGroup]
  nct <- length(nctr)
  nctp <- sum(nctr)
  ntrep0 <- c(nttp, nctp, ndummy)
  rtrep <- ntrep0[ntrep0 > 0]
  ntype <- length(rtrep)
  bsr <- lapply(blockSequence, function(x, n1, r1, n2, c1) {
    c(ceiling(n1/ceiling(r1/x[1])), ceiling(n2/ceiling(c1/x[2])))
  }, n1 = nrrd, r1 = rowsInDesign, n2 = ncrd, c1 = columnsInDesign)
  ribsr <- ceiling(nrrd/ceiling(rowsInDesign/rowsInBlockSequence))
  cibsr <- ceiling(ncrd/ceiling(columnsInDesign/columnsInBlockSequence))
  mrt <- matrix(1, nrrd, ncrd)
  nblkf <- length(ribsr)
  if (ntype > 1) {
    dtt <- corDiGGer(ntype, nrrd, ncrd, treatRep = rtrep, 
                     blockSequence = bsr, spatial = FALSE, rowColumn = FALSE, 
                     maxInterchanges = maxInterchanges, rngSeeds = rngSeeds, 
                     rngState = rng, runSearch = FALSE, aType = "A++")
    dtt <- run(dtt)
    mrt <- getDesign(dtt)
    rng <- dtt$.rng
  }
  else {
  }
  srt <- mrt
  if (ndummy > 0) 
    mrt[srt == ntype] <- 0
  if (ntt > 0) {
    ri <- ranmar(length(rep(seq(ntt), nttr)), rng)
    mrt[srt == 1] <- rep(seq(ntt), nttr)[order(ri$rvec)]
    rng <- ri$rng
    if (nct > 0) {
      ri <- ranmar(length(rep(seq(nct), nctr)), rng)
      mrt[srt == 2] <- rep(seq(nct), nctr)[order(ri$rvec)] + 
        ntt
      rng <- ri$rng
    }
  }
  else {
    ri <- ranmar(length(rep(seq(nct), nctr)), rng)
    mrt[srt == 1] <- rep(seq(nct), nctr)[order(ri$rvec)]
    rng <- ri$rng
  }
  trepr <- c(table(mrt[mrt != 0]))
  rtrt <- ntt + nct
  if (rtrt == 1) {
    mrto <- mrt
  }
  else {
    drt <- corDiGGer(rtrt, nrrd, ncrd, treatRep = trepr, 
                     initialDesign = mrt, initialSwap = srt, blockSequence = bsr, 
                     spatial = FALSE, rowColumn = FALSE, rngState = rng)
    mrto <- getDesign(drt)
    rngS <- drt$.rng
  }
  mrto[mrto == 0] <- rtrt + 1
  trtno <- c(seq(treatRepPerRep)[treatRepPerRep > 1 & treatGroup == 
                                   1], seq(treatRepPerRep)[treatRepPerRep > 1 & treatGroup != 
                                                             1], 0)
  mrtf <- matrix(trtno[mrto], nrrd)
  mfull <- matrix(0, rowsInDesign, columnsInDesign)
  fblockid <- findblks(rowsInDesign, columnsInDesign, nrtb, 
                       nctb)
  tblockid <- findblks(nrrd, ncrd, reptb, 1)
  for (i in seq(ntrb * ntcb)) {
    tmp <- mfull[fblockid == i]
    tmp[seq(reptb)] <- mrtf[tblockid == i]
    mfull[fblockid == i] <- tmp
  }
  utrt <- seq(treatRepPerRep)[treatRepPerRep == 1]
  ri <- ranmar(length(utrt), rng)
  mfull[mfull == 0] <- utrt[order(ri$rvec)]
  rng <- ri$rng
  for (i in seq(ntrb * ntcb)) {
    ri <- ranmar(length(mfull[fblockid == i]), rng)
    mfull[fblockid == i] <- mfull[fblockid == i][order(ri$rvec)]
    rng <- ri$rng
  }
  dfull <- DiGGer(numberOfTreatments, rowsInDesign, columnsInDesign, 
                  treatRep = treatRepPerRep, treatGroup = treatGroup, 
                  treatName = treatName, treatNumber = treatNumber, initialDesign = mfull, 
                  initialSwap = c(nrtb, nctb), rngState = rng)
  c1 <- Correlation(betweenRowCorrModel, betweenRowCorr, betweenColumnCorrModel, 
                    betweenColumnCorr, 1)
  addPhase(dfull, newphase = Phase(rowsInSwapBlock = nrtb, 
                                   columnsInSwapBlock = nctb, rowsInCorrelationBlock = rowsInDesign, 
                                   columnsInCorrelatio = columnsInDesign, maxInterchanges = maxInterchanges, 
                                   objectives = list(Objective(numberOfBlocks = 3, blocks = list(Block(1, 
                                                                                                       columnsInDesign, 1), Block(rowsInDesign, 1, 1), 
                                                                                                 Block()), corr = c1))))
  if (targetGroup == 1) 
    setAType(dfull, 1, "A11")
  if (targetGroup == 2) 
    setAType(dfull, 1, "A22")
  setSeeds(dfull, rngSeeds[1], rngSeeds[2])
  if (runSearch) {
    if (is.null(splitOpt)) {
      dfull <- run(dfull)
    }
    else {
      if (length(splitOpt) != 2) {
        print("#######################################")
        print("# Dimension pair expected in splitOpt #")
        print("# Final search has not yet been run   #")
        print("#######################################")
      }
      else {
        if (rowsInDesign%%splitOpt[1] != 0 | columnsInDesign%%splitOpt[2] != 
            0) {
          print("#######################################")
          print("# splitOpt blocks must tile evenly    #")
          print("# Final search has not yet been run   #")
          print("#######################################")
        }
        else {
          dfull <- subOpt(dfull, splitOpt[1], splitOpt[2])
        }
      }
    }
  }
  else {
    print("#####################################")
    print("# Final search has not yet been run #")
    print("#####################################")
  }
  dfull
}


############################################

# ------------------------------------------------------------------------------------------
#' @name ibDiGGer  
#' @aliases ibDiGGer  
#' @title ibDiGGer  
#' 
#' @description ibDiGGer modifies a DiGGer object to search for an incomplete block design.
#' @param numberOfTreatments The number of treatments in the design.
#' @param rowsInDesign The number of rows in the design.
#' @param columnsInDesign The number of columns in the design.
#' @param rowsInReplicate The number or rows in the template replicate block.
#' @param columnsInReplicate The number of columns in the template replicate block.
#' @param rowsInBlock The number of rows in each block.
#' @param columnsInBlock The number of columns in each block.
#' @param fixedBlocks Logical value whether blocks are treated as fixed.
#' @param blockGamma Variance component assigned to blocks.
#' @param aType A DiGGer A-measure type: "A++", "Agg", "A22", "A11", "A1+", "Aa2" or "Aaa". The default is "A++" for equally replicated treatments, "Agg" for unequal replication.
#' @param targetAValue The search stops in each search phase when the A value is below this value.
#' @param maxInterchanges Number of treatment interchanges to test
#' @param searchIntensity Percentage of possible interchanges to consider for non-improving interchanges.
#' @param runSearch Logical value, whether to run the search immediately.
#' @param rngSeeds Seeds c(s1,s2) to control the DiGGer search. s1 must be in the range [0,31328], s2 must be in the range [0,30081].
#' @param rngState Current state of the random number generator.
#' @param treatName Vector of treatment names to be associated with the design numbers.
#' @param treatNumber Vector of treatment numbers to be associated with the design numbers.
#' @param treatRepPerRep Vector of replication levels for each treatment within each replicate template block.
#' @param treatGroup Vector of group codes (up to 200 distinct values) associated with treatments. Group codes may be used to modify the A-efficiency measure in the optimisation.
#' @param initialDesign A matrix rowsInDesign by columnsInDesign giving design numbers in the initial design.
#' @param initialSwap A dimension pair or a matrix rowsInDesign by columnsInDesign of swap codes. Only plots with the same swap code may have treatment interchanges during the DiGGer search.
#' @source DiGGer 1.0.5
# ------------------------------------------------------------------------------------------

ibDiGGer <- function (numberOfTreatments, rowsInDesign, columnsInDesign, 
                      rowsInReplicate = NULL, columnsInReplicate = NULL, rowsInBlock, 
                      columnsInBlock, fixedBlocks = FALSE, blockGamma = 1, aType = NULL, 
                      targetAValue = 0, maxInterchanges = 1e+05, searchIntensity = 100, 
                      runSearch = TRUE, rngSeeds = NULL, rngState = NULL, treatName = NULL, 
                      treatNumber = NULL, treatRepPerRep = NULL, treatGroup = NULL, 
                      initialDesign = NULL, initialSwap = NULL) 
{
  if (fixedBlocks) 
    blockGamma <- 0
  if (is.null(treatRepPerRep)) {
    treatRepPerRep <- rep(1, numberOfTreatments)
  }
  nrbv <- rowsInBlock
  ncbv <- columnsInBlock
  nmissing <- 0
  if (!is.null(initialDesign)) {
    ttab <- table(initialDesign)
    if (names(ttab)[1] == "0") {
      nmissing <- ttab[1]
      nrbv <- rowsInDesign
      ncbv <- columnsInDesign
    }
  }
  if (is.null(aType)) {
    ttab <- table(treatRepPerRep)
    if (length(ttab) == 1) 
      aType <- "A++"
    else aType <- "Agg"
  }
  dib <- DiGGer(numberOfTreatments = numberOfTreatments, rowsInDesign = rowsInDesign, 
                columnsInDesign = columnsInDesign, rowsInReplicate = rowsInReplicate, 
                columnsInReplicate = columnsInReplicate, treatName = treatName, 
                treatNumber = treatNumber, treatRepPerRep = treatRepPerRep, 
                treatGroup = treatGroup, initialDesign = initialDesign, 
                rngSeeds = rngSeeds, rngState = rngState)
  addPhase(dib, newphase = Phase(rowsInSwapBlock = dib$.nrbr, 
                                 columnsInSwapBlock = dib$.ncbr, rowsInCorrelationBlock = nrbv, 
                                 columnsInCorrelationBlock = ncbv, maxInterchanges = maxInterchanges, 
                                 searchIntensity = searchIntensity, aType = aType, targetAValue = targetAValue, 
                                 objectives = list(Objective(numberOfBlocks = 1, blocks = list(Block(rowsInBlock, 
                                                                                                     columnsInBlock, blockGamma)), corr = Correlation("ID", 
                                                                                                                                                      0, "ID", 0, 1)))))
  if (runSearch) 
    dib <- run(dib)
  dib
}

############################################

# ------------------------------------------------------------------------------------------
#' @name rcDiGGer  
#' @aliases rcDiGGer  
#' @title rcDiGGer  
#' 
#' @description rcDiGGer modifies a DiGGer object to search for a row-column design.
#' @param numberOfTreatments The number of treatments in the design.
#' @param rowsInDesign The number of rows in the design.
#' @param columnsInDesign The number of columns in the design.
#' @param rowsInReplicate The number or rows in the template replicate block.
#' @param columnsInReplicate The number of columns in the template replicate block.
#' @param twoPhase If TRUE the bigger block is optimised in the first phase of the search and the smaller in the second phase. The order of the optimisation can be specified as "rowThenCol" or "colThenRow". If FALSE both blocks are optimised in one search phase.
#' @param nested If TRUE blocks do not extend past replicates. If a pair of numbers given, the first is the number of rows in a column-block, the second is the number of columns in a row-block. If nested is FALSE, row and column blocks extend across the whole design.
#' @param fixedBlocks Logical value whether blocks are treated as fixed.
#' @param blockGammas Variance component assigned to random row blocks, rgamma and random column blocks, cgamma.
#' @param aType A DiGGer A-measure type: "A++", "Agg", "A22", "A11", "A1+", "Aa2" or "Aaa". The default is "A++" for equally replicated treatments, "Agg" for unequal replication.
#' @param targetAValue The search stops in each search phase when the A value is below this value.
#' @param maxInterchanges number of treatment interchanges to test in each search phase. If twoPhase = TRUE two numbers may be given.
#' @param searchIntensity Percentage of possible interchanges to consider for non-improving interchanges. Two figures may be given for two phase searches.
#' @param runSearch Logical value, whether to run the search immediately.
#' @param rngSeeds Seeds c(s1,s2) to control the DiGGer search. s1 must be in the range [0,31328], s2 must be in the range [0,30081].
#' @param rngState Current state of the random number generator.
#' @param treatName Vector of treatment names to be associated with the design numbers.
#' @param treatNumber Vector of treatment numbers to be associated with the design numbers.
#' @param treatRepPerRep Vector of replication levels for each treatment within each replicate template block.
#' @param treatGroup Vector of group codes (up to 200 distinct values) associated with treatments. Group codes may be used to modify the A-efficiency measure in the optimisation.
#' @param initialDesign A matrix rowsInDesign by columnsInDesign giving design numbers in the initial design.
#' @param initialSwap A dimension pair or a matrix rowsInDesign by columnsInDesign of swap codes. Only plots with the same swap code may have treatment interchanges during the DiGGer search.
#' @source DiGGer 1.0.5
# ------------------------------------------------------------------------------------------

rcDiGGer <- function (numberOfTreatments, rowsInDesign, columnsInDesign, 
                      rowsInReplicate = NULL, columnsInReplicate = NULL, twoPhase = TRUE, 
                      nested = TRUE, fixedBlocks = FALSE, blockGammas = list(rgamma = 1, 
                                                                             cgamma = 1), aType = NULL, targetAValue = 0, maxInterchanges = 1e+05, 
                      searchIntensity = c(0, 100), runSearch = TRUE, rngSeeds = NULL, 
                      rngState = NULL, treatName = NULL, treatNumber = NULL, treatRepPerRep = NULL, 
                      treatGroup = NULL, initialDesign = NULL, initialSwap = NULL) 
{
  if (length(maxInterchanges) == 2) {
    maxI1 <- maxInterchanges[1]
    maxI2 <- maxInterchanges[2]
  }
  else {
    maxI1 <- maxInterchanges
    maxI2 <- maxInterchanges
  }
  if (length(searchIntensity) == 2) {
    srchI1 <- searchIntensity[1]
    srchI2 <- searchIntensity[2]
  }
  else {
    srchI1 <- searchIntensity
    srchI2 <- searchIntensity
  }
  if (length(targetAValue) == 1) {
    targetAValue <- c(targetAValue, targetAValue)
  }
  drc <- DiGGer(numberOfTreatments = numberOfTreatments, rowsInDesign = rowsInDesign, 
                columnsInDesign = columnsInDesign, rowsInReplicate = rowsInReplicate, 
                columnsInReplicate = columnsInReplicate, treatName = treatName, 
                treatNumber = treatNumber, treatRepPerRep = treatRepPerRep, 
                treatGroup = treatGroup, initialDesign = initialDesign, 
                rngSeeds = rngSeeds, rngState = rngState)
  s0 <- getSeeds(drc)
  if (length(nested) == 2) {
    rowsInBlock <- nested[1]
    columnsInBlock <- nested[2]
  }
  else {
    if (!is.logical(nested)) {
      stop("TRUE/FALSE required for single value of nested.")
    }
    if (nested) {
      if (drc$.nrbr == 1 | drc$.ncbr == 1) {
        stop("rows and column blocks must have more than 1 unit")
      }
      rowsInBlock <- drc$.nrbr
      columnsInBlock <- drc$.ncbr
    }
    else {
      rowsInBlock <- rowsInDesign
      columnsInBlock <- columnsInDesign
    }
  }
  nrbv <- rowsInBlock
  ncbv <- columnsInBlock
  nmissing <- 0
  if (!is.null(initialDesign)) {
    ttab <- table(initialDesign)
    if (names(ttab)[1] == "0") {
      nmissing <- ttab[1]
      nrbv <- rowsInDesign
      ncbv <- columnsInDesign
    }
  }
  if (fixedBlocks) {
    cblock <- Block(rowsInBlock, 1, 0)
    rblock <- Block(1, columnsInBlock, 0)
  }
  else {
    cblock <- Block(rowsInBlock, 1, blockGammas$cgamma)
    rblock <- Block(1, columnsInBlock, blockGammas$rgamma)
  }
  corr0 <- Correlation(rowModel = "ID", rowParameter = 0, 
                       columnModel = "ID", columnParameter = 0, spatialVarianceRatio = 1)
  if (is.logical(twoPhase)) {
    if (twoPhase) {
      if (rowsInBlock < columnsInBlock) {
        twoPhase <- "colThenRow"
      }
      else {
        twoPhase <- "rowThenCol"
      }
    }
  }
  if (is.null(aType)) {
    ttab <- table(treatRepPerRep)
    if (length(ttab) == 1) {
      aType <- "A++"
    }
    else {
      aType <- "Agg"
    }
  }
  if (twoPhase == "rowThenCol") {
    nrbv1 <- 1
    ncbv1 <- columnsInBlock
    if (nmissing > 0) {
      nrbv1 <- nrbv
      ncbv1 <- ncbv
    }
    addPhase(drc, newphase = Phase(rowsInSwapBlock = drc$.nrbr, 
                                   columnsInSwapBlock = drc$.ncbr, rowsInCorrelationBlock = nrbv1, 
                                   columnsInCorrelationBlock = ncbv1, maxInterchanges = maxI1, 
                                   searchIntensity = srchI1, objectives = list(Objective(numberOfBlocks = 1, 
                                                                                         blocks = list(rblock), corr = corr0))))
    nrbv2 <- rowsInBlock
    ncbv2 <- columnsInBlock
    if (nmissing > 0) {
      nrbv2 <- rowsInDesign
      ncbv2 <- columnsInDesign
    }
    addPhase(drc, newphase = Phase(rowsInSwapBlock = 1, 
                                   columnsInSwapBlock = columnsInBlock, rowsInCorrelationBlock = nrbv2, 
                                   columnsInCorrelationBlock = ncbv2, maxInterchanges = maxI2, 
                                   searchIntensity = srchI2, objectives = list(Objective(numberOfBlocks = 2, 
                                                                                         blocks = list(rblock, cblock), corr = corr0))))
  }
  else if (twoPhase == "colThenRow") {
    nrbv1 <- rowsInBlock
    ncbv1 <- 1
    if (nmissing > 0) {
      nrbv1 <- rowsInDesign
      ncbv1 <- columnsInDesign
    }
    addPhase(drc, newphase = Phase(rowsInSwapBlock = drc$.nrbr, 
                                   columnsInSwapBlock = drc$.ncbr, rowsInCorrelationBlock = nrbv1, 
                                   columnsInCorrelationBlock = ncbv1, maxInt = maxI1, 
                                   searchInt = srchI1, objectives = list(Objective(numberOfBlocks = 1, 
                                                                                   blocks = list(cblock), corr = corr0))))
    nrbv2 <- rowsInBlock
    ncbv2 <- columnsInBlock
    if (nmissing > 0) {
      nrbv2 <- rowsInDesign
      ncbv2 <- columnsInDesign
    }
    addPhase(drc, newphase = Phase(rowsInSwapBlock = rowsInBlock, 
                                   columnsInSwapBlock = 1, rowsInCorrelationBlock = nrbv2, 
                                   columnsInCorrelationBlock = ncbv2, maxInterchanges = maxI2, 
                                   searchIntensity = srchI2, objectives = list(Objective(numberOfBlocks = 2, 
                                                                                         blocks = list(rblock, cblock), corr = corr0))))
  }
  else {
    nrbv2 <- rowsInBlock
    ncbv2 <- columnsInBlock
    if (nmissing > 0) {
      nrbv2 <- rowsInDesign
      ncbv2 <- columnsInDesign
    }
    addPhase(drc, newphase = Phase(rowsInSwapBlock = drc$.nrbr, 
                                   columnsInSwapBlock = drc$.ncbr, rowsInCorrelationBlock = nrbv2, 
                                   columnsInCorrelationBlock = ncbv2, maxInterchanges = maxI1, 
                                   searchIntensity = srchI1, objectives = list(Objective(numberOfBlocks = 2, 
                                                                                         blocks = list(rblock, cblock), corr = corr0))))
  }
  nphases <- getNumberOfSearchPhases(drc)
  for (i in seq(nphases)) {
    setTargetAValue(drc, i, targetAValue[i])
  }
  if (runSearch) {
    drc <- run(drc)
  }
  drc
}

