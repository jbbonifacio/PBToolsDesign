# -------------------------------------------------------------------------------------
# Name             : designPRep 
# Description      : Generate randomization for p-rep design
# R Version        : 4.0.1 
# -------------------------------------------------------------------------------------
# Author           : Alaine A. Gulles 
# Author Email     : a.gulles@irri.org
# Script Version   : 3.3.0
# Modified by      : Justine B. Bonifacio
# Date Modified    : 2023.05.11
# Changes made     : remove the rule for replication of 20% test entries, remove blocks
# -------------------------------------------------------------------------------------

#' @title P-Rep Design
#' @aliases designPrep
#  @aliases designPrep.default
#'  
#' @description Generate randomization and layout for P-Rep Design.
#'
#' @param numTrmtPerTrmtGrp number of treatments per replicate group
#' @param numRepPerTrmtGrp number of replicates per treatment group
#' @param numTrmtPerGrp number of treatments per group
#' @param numTrial number of trials (randomization set-ups with the complete set of entries)
#' @param genLayout logical, whether a layout of the design will be generated or not
#' @param numFieldRow number of field rows
#' @param treatName NULL or vector of names associated with the treatment
#' @param listName column name associated with the list provided (treatName)
#' @param serpentine logical, whether plot number will be arranged as serpentine order
#' @param topToBottom logical, whether plot number will be written from to to bottom
#' @param display logical, whether randomization parameters will be displayed
#' 
#' @details 
#' If \code{genLayout} is \code{TRUE}, then parameters \code{numFieldRow} and \code{numFieldCol} must be specified.  
#' Values of \code{numFieldRow} and \code{numFieldCol} should be a factor of the total number of plots.
#' 
#' @return A list containing the following components:
#' \item{fieldbook}{a data frame}
#' \item{plan}{a data frame, if \code{genLayout} is \code{FALSE} or a list containing the following components 
#' if \code{genLayout} is \code{TRUE}}
#' \item{TrmtLayout}{a list whose length is equal to the 
#'      number of trials containing the treatment layout for each trial} 
#' \item{PlotNumLayout}{a matrix containing the plot number layout of the experiment}
#' 
#' @examples
#' ## Generate randomization for an experiment with 150 treatment levels in P-rep design.
#' ## There will be 2 groups (test and check)
#' ## Test group 1 - 175 entries replicated 1 time
#' ## Check group 1 - 5 entries replicated 5 times
#' prep1 <- designPrep(numRepPerTrmtGrp = c(1,5), numTrmtPerTrmtGrp = c(175,5), numTrmtPerGrp = c(175,5), numTrial = 1,
#'                    genLayout = FALSE, treatName = NULL, listName = "Entry",
#'                    display = TRUE)
#'                    
#' ## The experiment will be arranged in a 20 x 10 field in Left to Right Serpentine order.
#' prep2 <- designPrep(numRepPerTrmtGrp = c(1,2,5), numTrmtPerTrmtGrp = c(115,30,5), numTrmtPerGrp = c(145,5), numTrial = 1,
#'                    genLayout = TRUE, numFieldRow = 20, treatName = NULL, listName = "Entry",
#'                    serpentine = TRUE, topToBottom = FALSE, display = TRUE)
# -------------------------------------------------------------------------------------  

# install.packages("dplyr")
# install.packages("stringr")
# library("dplyr")
# library("stringr")

designPrep <- function(numTrmtPerTrmtGrp, numRepPerTrmtGrp, numTrmtPerGrp, numTrial = 1, 
                       genLayout = FALSE, numFieldRow = 1, 
                       treatName = NULL, listName = NULL, serpentine = FALSE, 
                       topToBottom = TRUE, display = TRUE){ #UseMethod("designPrep") }

#designPrep.default <- function(numTrmtPerTrmtGrp, numRepPerTrmtGrp, numTrmtPerGrp, numTrial = 1, 
#                       genLayout = FALSE, numFieldRow = 1, 
#                       treatName = NULL, listName = NULL, serpentine = FALSE, 
#                       topToBottom = TRUE, display = TRUE) {

    # --- check input from user:
    
    # default values when genLayout is TRUE or FALSE
    treatRepPerRep <- rep(numRepPerTrmtGrp,numTrmtPerTrmtGrp)
    treatGroup <- rep(c(1,2),numTrmtPerGrp)
    numberOfTreatments <- sum(numTrmtPerTrmtGrp)
    
    if(!genLayout){
        numFieldRow <- sum(treatRepPerRep)
        numFieldCol <- 1
        tempblk <- getMidFloorFactor(numFieldRow)
        blockSequence <- list(c(tempblk,1))
    } else {
      if(sum(treatRepPerRep)%%numFieldRow == 0){
        numFieldCol <- sum(treatRepPerRep)/numFieldRow
      } else {
        numFieldCol <- ceiling(sum(treatRepPerRep)/numFieldRow)
      }
      
      if(numFieldCol <= numFieldRow){
        tempblk <- getMidFloorFactor(numFieldRow)
        blockSequence <- list(c(tempblk,numFieldCol))
      } else{
        tempblk <- getMidFloorFactor(numFieldCol)
        blockSequence <- list(c(numFieldRow,tempblk))
      }
    }
    
    #assign group names
    cumTreatPerRep <- cumsum(numTrmtPerTrmtGrp)
    j <- 0
    for(i in 1:length(numTrmtPerTrmtGrp)){
      if(cumTreatPerRep[i] <= numTrmtPerGrp[1]){
        names(numTrmtPerTrmtGrp)[i] <- paste("Test Group",i,sep=" ")
        j <- j+1
      } else { names(numTrmtPerTrmtGrp)[i] <- paste("Check Group",i-j,sep=" ")}
    }
    
    treatList <- TRUE
    # default value when treatName is NULL
    if(is.null(treatName)){
      treatName <- c(paste("Test",(1:numTrmtPerGrp[1]),sep=""),paste("Check",1:numTrmtPerGrp[2],sep=""))
      treatList <- FALSE
    }
  
    #check if rep more than 1 in check groups
    grpchck <- as.data.frame(cbind(treatRepPerRep, treatGroup))
    if(any(grpchck$treatGroup == 2 & grpchck$treatRepPerRep < 2)){
      stop("ERROR: The number of replicates for check groups should be at least 2.")}
    
    # #check if rep more than 1 if only 1 test group
    # if(length(grpchck %>% filter(treatGroup == 1) %>% unique() %>% .$treatRepPerRep) == 1){
    #   if((grpchck %>% filter(treatGroup == 1) %>% unique() %>% .$treatRepPerRep) < 2){
    #     stop("ERROR: If there is only 1 test group, the number of replicates should be at least 2.")}
    # }
    
    #check if there is a duplicate in the number of replicates per group
    if(any(duplicated(numRepPerTrmtGrp[1:which(cumsum(numTrmtPerTrmtGrp) == numTrmtPerGrp[1])]) == TRUE)){
      stop("ERROR: The number of replicates for test groups should not be duplicated.")}
    
    if(any(duplicated(numRepPerTrmtGrp[(which(cumsum(numTrmtPerTrmtGrp) == numTrmtPerGrp[1])+1):length(numTrmtPerTrmtGrp)]) == TRUE)){
      stop("ERROR: The number of replicates for check groups should not be duplicated.")}
    
    # #check if at least 20% of test entries replicated
    # replicated <- nrow(grpchck %>% filter(treatGroup == 1 & treatRepPerRep > 1))
    # totalTest <- nrow(grpchck %>% filter(treatGroup == 1))
    # percRep <- (replicated / totalTest) * 100
    # if (percRep < 20) { 
    #   stop("ERROR: At least 20% of the total number of treatment levels of the test group should be replicated.") }
    
    # #check if total eu is not a prime number
    # flag <- primeNumber(sum(treatRepPerRep))
    # if(flag == 1) {
    #   stop("ERROR: The total number of plots should not be a prime number.")}
  
    #instantiate values
    randomize <- NULL
    plan <- list()
    plotNum <- NULL
    trmtLayout <- NULL
    
    #assign colors per group in the plot layout
    #if(length(numTrmtPerTrmtGrp) <= 8) { colorCode <- c(8:(8-length(numTrmtPerTrmtGrp) + 1))     
    #} else { colorCode <- rev(c(1:length(numTrmtPerTrmtGrp))) } 
    
    for (i in (1:numTrial)) {
      
      if(sum(treatRepPerRep) > 1500 | primeNumber(sum(treatRepPerRep)) == 1){
        if(sum(treatRepPerRep)%%numFieldRow != 0){
          if(topToBottom){
            # initial matrix
            initialMat <- matrix(data = 1, nrow = numFieldRow, ncol = numFieldCol)
            emptyCol <- ((numFieldRow*numFieldCol)-sum(treatRepPerRep))%/%numFieldRow
            notEmptyCol <- numFieldCol - emptyCol
            
            if(emptyCol > 0){
              initialMat[1:numFieldRow, (numFieldCol-emptyCol+1):numFieldCol] <- 0
            }
            
            if((serpentine & notEmptyCol%%2 != 0) | !serpentine){
              initialMat[(numFieldRow - (((numFieldRow*numFieldCol)-sum(treatRepPerRep))%%numFieldRow) + 1):numFieldRow, numFieldCol-emptyCol] <- 0
            } else{
              initialMat[1:(((numFieldRow*numFieldCol)-sum(treatRepPerRep))%%numFieldRow), numFieldCol-emptyCol] <- 0
            }  
          } else{
            # initial matrix
            initialMat <- matrix(data = 1, nrow = numFieldRow, ncol = numFieldCol)
            emptyRow <- ((numFieldRow*numFieldCol)-sum(treatRepPerRep))%/%numFieldCol
            notEmptyRow <- numFieldRow - emptyRow
            
            if(emptyRow > 0){
              initialMat[(numFieldRow-emptyRow+1):numFieldRow, 1:numFieldCol] <- 0
            }
            
            if((serpentine & notEmptyRow%%2 != 0) | !serpentine){
              initialMat[numFieldRow-emptyRow, (numFieldCol - (((numFieldRow*numFieldCol)-sum(treatRepPerRep))%%numFieldCol) + 1):numFieldCol] <- 0
            } else{
              initialMat[numFieldRow-emptyRow, 1:(((numFieldRow*numFieldCol)-sum(treatRepPerRep))%%numFieldCol)] <- 0
            }
          }
          
          # initial design
          choices <- rep(treatName,treatRepPerRep)
          iDesign <- initialMat
          iDesign[iDesign == 1] <- sample(choices, length(choices))
        } else {
          choices <- rep(treatName,treatRepPerRep)
          iDesign <- sample(choices, length(choices))
        }
        
        trmt <- data.frame("ID" = treatName, "ENTRY" = 1:numberOfTreatments)
        result <- data.frame("UNIT" = 1:(numFieldRow*numFieldCol), "ID" = as.vector(iDesign))
        result$ENTRY <- trmt$ENTRY[match(result$ID,trmt$ID)]
        
        prep_mat <- matrix(result$ENTRY, nrow = numFieldRow, ncol = numFieldCol, byrow = FALSE)
        prep <- data.frame(result,
                           "ROW" = rep(1:numFieldRow, numFieldCol),
                           "RANGE" = rep(1:numFieldCol, each = numFieldRow),
                           "REP" = rep(1,(numFieldRow*numFieldCol)),
                           "TRT" = result$ID)
        
        plan[[i]] <- matrix(prep$ID, nrow(prep_mat), ncol(prep_mat))
        
      } else {
        if(sum(treatRepPerRep)%%numFieldRow != 0){
          if(topToBottom){
            # initial matrix
            initialMat <- matrix(data = 1, nrow = numFieldRow, ncol = numFieldCol)
            emptyCol <- ((numFieldRow*numFieldCol)-sum(treatRepPerRep))%/%numFieldRow
            notEmptyCol <- numFieldCol - emptyCol
            
            if(emptyCol > 0){
              initialMat[1:numFieldRow, (numFieldCol-emptyCol+1):numFieldCol] <- 0
            }
            
            if((serpentine & notEmptyCol%%2 != 0) | !serpentine){
              initialMat[(numFieldRow - (((numFieldRow*numFieldCol)-sum(treatRepPerRep))%%numFieldRow) + 1):numFieldRow, numFieldCol-emptyCol] <- 0
            } else{
              initialMat[1:(((numFieldRow*numFieldCol)-sum(treatRepPerRep))%%numFieldRow), numFieldCol-emptyCol] <- 0
            }  
          } else{
            # initial matrix
            initialMat <- matrix(data = 1, nrow = numFieldRow, ncol = numFieldCol)
            emptyRow <- ((numFieldRow*numFieldCol)-sum(treatRepPerRep))%/%numFieldCol
            notEmptyRow <- numFieldRow - emptyRow
            
            if(emptyRow > 0){
              initialMat[(numFieldRow-emptyRow+1):numFieldRow, 1:numFieldCol] <- 0
            }
            
            if((serpentine & notEmptyRow%%2 != 0) | !serpentine){
              initialMat[numFieldRow-emptyRow, (numFieldCol - (((numFieldRow*numFieldCol)-sum(treatRepPerRep))%%numFieldCol) + 1):numFieldCol] <- 0
            } else{
              initialMat[numFieldRow-emptyRow, 1:(((numFieldRow*numFieldCol)-sum(treatRepPerRep))%%numFieldCol)] <- 0
            }
          }
          
          # initial design
          choices <- rep(treatName,treatRepPerRep)
          iDesign <- initialMat
          iDesign[iDesign == 1] <- sample(choices, length(choices))
          
          trmt <- data.frame("ID" = treatName, "ENTRY" = 1:numberOfTreatments)
          result <- data.frame("UNIT" = 1:(numFieldRow*numFieldCol), "ID" = as.vector(iDesign))
          result$ENTRY <- trmt$ENTRY[match(result$ID,trmt$ID)]
          
          prep_mat <- matrix(result$ENTRY, nrow = numFieldRow, ncol = numFieldCol, byrow = FALSE)
          prep <- data.frame(result,
                             "ROW" = rep(1:numFieldRow, numFieldCol),
                             "RANGE" = rep(1:numFieldCol, each = numFieldRow),
                             "REP" = rep(1,(numFieldRow*numFieldCol)),
                             "TRT" = result$ID)
          
          plan[[i]] <- matrix(prep$ID, nrow(prep_mat), ncol(prep_mat))
        } else{
          sink(tempfile())
          # Using prDiGGer function
          result <- try(prep <- DiGGer::prDiGGer(numberOfTreatments = numberOfTreatments,
                                                 treatRepPerRep = treatRepPerRep,
                                                 treatGroup = treatGroup,
                                                 blockSequence = blockSequence,
                                                 rowsInDesign = numFieldRow,
                                                 columnsInDesign = numFieldCol,
                                                 treatName = treatName,
                                                 runSearch = TRUE),
                        silent = TRUE)
          
          # Using nurseryDiGGer
          # result <- try(prep <- nurseryDiGGer(numberOfTreatments = numberOfTreatments,
          #                                     treatRepPerRep = treatRepPerRep,
          #                                     treatGroup = treatGroup,
          #                                     blockSequence = blockSequence,
          #                                     rowsInDesign = numFieldRow,
          #                                     columnsInDesign = numFieldCol,
          #                                     treatName = treatName,
          #                                     checkGroup = 2),
          #               silent = TRUE)
          sink()
          if(all(class(result) == "try-error")) {
            msg <- str_trim(strsplit(result, ":")[[1]])
            msg <- str_trim(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
            stop(paste("Error in DiGGer:", msg, sep = ""))
          }
          
          capture.output(prep)
          if(all(class(result) == "try-error")) {
            msg <- str_trim(strsplit(result, ":")[[1]])
            msg <- str_trim(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
            stop(paste("Error in DiGGer:", msg, sep = ""))
          }
          
          prep_mat <- getDesign(prep)
          plan[[i]] <- matrix(print(prep, option = "list")$ID, nrow(prep_mat), ncol(prep_mat))
        }
      }
      
      if (!topToBottom){
        if (i == 1){
          plotNum <- matrix(1:(numFieldRow*numFieldCol), nrow(plan[[i]]), ncol(plan[[i]]), byrow = TRUE)
          if (serpentine & numFieldRow > 1) { 
            for(k in seq(2, numFieldRow, by = 2)) { 
              plotNum[k,] <- rev(plotNum[k,]) }}
        } ## end stmt -- without blk if (i == 1)
        tempFieldOrder <- as.data.frame.table(plotNum)
      } else{
        if (i == 1){
          plotNum <- matrix(1:(numFieldRow*numFieldCol), nrow(plan[[i]]), ncol(plan[[i]]), byrow = FALSE)
          if (serpentine & numFieldCol > 1) { 
            for(k in seq(2, numFieldCol, by = 2)) { 
              plotNum[,k] <- rev(plotNum[,k]) }}
        } ## end stmt -- without blk if (i == 1)
        tempFieldOrder <- as.data.frame.table(plotNum)
      }
      
      tempFieldOrder[,"Var1"] <- as.numeric(tempFieldOrder[,"Var1"])
      tempFieldOrder[,"Var2"] <- as.numeric(tempFieldOrder[,"Var2"])
      
      if(sum(treatRepPerRep) > 1500 | primeNumber(sum(treatRepPerRep)) == 1 | (sum(treatRepPerRep)%%numFieldRow) != 0){
        randomize <- rbind(randomize, cbind(TRIAL = i, merge(prep, tempFieldOrder, by.x = c("ROW", "RANGE"), by.y = c("Var1", "Var2"))))
      } else{
        randomize <- rbind(randomize, cbind(TRIAL = i, merge(print(prep, option = "list"), tempFieldOrder, by.x = c("ROW", "RANGE"), by.y = c("Var1", "Var2"))))
      }
      
      if(genLayout){
        trmtLayout[[i]] <- plan[[i]]
        dimnames(trmtLayout[[i]]) <- list(paste("FieldRow", 1:nrow(plan[[i]]), sep = ""), paste("FieldColumn", 1:ncol(plan[[i]]), sep = ""))
      } else{
        trmtLayout[[i]] <- data.frame(Trial = i, plan[[i]])
        dimnames(trmtLayout[[i]]) <- list(1:nrow(plan[[i]]),c("Trial","EntryName"))
      }
      #dimnames(plan1[[i]]) <- dimnames(plan[[i]])
    } ## end stmt -- for (i in (1:numTrial))
    
    randomize <- mutate(randomize, TYPE = if_else(as.numeric(randomize$ENTRY) <= numTrmtPerGrp[1], "Test", "Check"))
    
    #generating fieldbook
    if (!genLayout){
      #dimnames(plotNum) <- list(1:nrow(plan[[1]]),"PlotNumber")
      names(trmtLayout) <- paste("Trial", 1:numTrial, sep = "")
      #names(plan1) <- names(plan)
      
      if(treatList){
        randomize <- randomize[,c("TRIAL", "ID", "TYPE", "Freq")]
        colnames(randomize)[match("TRIAL", colnames(randomize))] <- "Trial"
        colnames(randomize)[match("Freq", colnames(randomize))] <- "PlotNumber"
        if(!is.null(listName)){
          colnames(randomize)[match("ID", colnames(randomize))] <- listName
          colnames(randomize)[match("TYPE", colnames(randomize))] <- paste(listName,"Type",sep="")
        } else {
          colnames(randomize)[match("ID", colnames(randomize))] <- "Entry"
          colnames(randomize)[match("TYPE", colnames(randomize))] <- "EntryType"
        }
      } else {
        randomize <- randomize[,c("TRIAL", "ID", "Freq")]
        colnames(randomize)[match("TRIAL", colnames(randomize))] <- "Trial"
        colnames(randomize)[match("Freq", colnames(randomize))] <- "PlotNumber"
        if(!is.null(listName)){
          colnames(randomize)[match("ID", colnames(randomize))] <- listName
        } else {
          colnames(randomize)[match("ID", colnames(randomize))] <- "Entry"
        }
      }
    } else{
      dimnames(plotNum) <- dimnames(trmtLayout[[1]])
      names(trmtLayout) <- paste("Trial", 1:numTrial, sep = "")
      #names(plan1) <- names(plan)
      
      if(treatList){
        randomize <- randomize[,c("TRIAL", "ID", "TYPE", "Freq", "ROW", "RANGE")]
        colnames(randomize)[match("TRIAL", colnames(randomize))] <- "Trial"
        colnames(randomize)[match("Freq", colnames(randomize))] <- "PlotNumber"
        colnames(randomize)[match("ROW", colnames(randomize))] <- "FieldRow"
        colnames(randomize)[match("RANGE", colnames(randomize))] <- "FieldColumn"
        if(!is.null(listName)) { 
          colnames(randomize)[match("ID", colnames(randomize))] <- listName
          colnames(randomize)[match("TYPE", colnames(randomize))] <- paste(listName,"Type",sep="")
        } else {
          colnames(randomize)[match("ID", colnames(randomize))] <- "Entry"
          colnames(randomize)[match("TYPE", colnames(randomize))] <- "EntryType"
        }
      } else {
        randomize <- randomize[,c("TRIAL", "ID", "Freq", "ROW", "RANGE")]
        colnames(randomize)[match("TRIAL", colnames(randomize))] <- "Trial"
        colnames(randomize)[match("Freq", colnames(randomize))] <- "PlotNumber"
        colnames(randomize)[match("ROW", colnames(randomize))] <- "FieldRow"
        colnames(randomize)[match("RANGE", colnames(randomize))] <- "FieldColumn"
        if(!is.null(listName)) { 
          colnames(randomize)[match("ID", colnames(randomize))] <- listName
        } else {
          colnames(randomize)[match("ID", colnames(randomize))] <- "Entry"
        }
      }
    }
    
    #combining appropriate plan layouts
    if(genLayout){
      exptLayout <- list()
      exptLayout[[1]] <- trmtLayout
      exptLayout[[2]] <- plotNum
      names(exptLayout) <- c("TrmtLayout", "PlotNumLayout")
    } else {
      exptLayout <- data.frame()
      exptLayout <- randomize[order(randomize$Trial, randomize$PlotNumber), ]
      rownames(exptLayout) <- 1:nrow(exptLayout)
      exptLayout <- exptLayout[1:2]
    }
    
    # re-arrange the fieldbook
    randomize <- randomize[order(randomize$Trial, randomize$PlotNumber), ]
    rownames(randomize) <- 1:nrow(randomize)
  
    randomize <- as.data.frame(randomize %>% 
                                 group_by_at(vars(names(randomize[2]),names(randomize[1]))) %>% 
                                 mutate(Rep=1:n(),.after = 1))
    
    # print of design parameters:
    if(display){
      cat(toupper("Design Properties:"), "\n", sep = "")
      cat("\t", "Incomplete Block Design", "\n", sep = "")
      cat("\t", "P-Rep Design", "\n\n", sep = "")
      cat(toupper("Design Parameters"), "\n", sep = "")
      cat("\t", "Number of Trials = ", numTrial, "\n", sep = "")
      cat("\t", "Number of Groupings = ", length(numTrmtPerTrmtGrp), "\n", sep = "")
      
      for (i in (1:length(numTrmtPerTrmtGrp))) {
        cat("\t", "Group ", i, ": ", names(numTrmtPerTrmtGrp)[i], "\n", sep = "")
        cat("\t", "Number of Replicate for ", names(numTrmtPerTrmtGrp)[i], " = ",  numRepPerTrmtGrp[i],"\n", sep = "")
        cat("\t", "Levels of ", names(numTrmtPerTrmtGrp)[i]," = ", sep = "")
        
        #if (treatList) {
          if (i == 1) {
            start <- 1 
            end <- numTrmtPerTrmtGrp[[1]]
          } else {
            start <- sum(numTrmtPerTrmtGrp[1:(i-1)])+1
            end <- start + numTrmtPerTrmtGrp[i]-1
          }
          if (numTrmtPerTrmtGrp[[i]] <= 5) { 
            cat(paste(treatName[start:end], collapse = ", ", sep = ""), "\n", sep = "") 
          } else {
            cat(paste(treatName[start:(start+2)], collapse = ", ", sep = ""), sep = "") 
            cat(", ..., ", paste(treatName[end], sep = ""), "\n", sep = "") 
          }
        #} else {
        #  if (numTrmtPerTrmtGrp[[i]] <= 5) { cat(paste(names(numTrmtPerTrmtGrp)[i], paste(1:numTrmtPerTrmtGrp[[i]]), collapse = ", ", sep = ""), "\n", sep = "") 
        #  } else {
        #    cat(paste(names(numTrmtPerTrmtGrp)[i], paste(1:3), collapse = ", ", sep = ""), sep = "") 
        #    cat(", ..., ", paste(names(numTrmtPerTrmtGrp)[i], numTrmtPerTrmtGrp[[i]], sep = ""), "\n", sep = "") 
        #  }
        #}
        
        cat("\n")
      }
      
      if(genLayout){
        cat("\t", "Number of Field Rows = ", numFieldRow,"\n", sep = "")
        cat("\t", "Number of Field Columns = ", numFieldCol, "\n\n", sep = "")}
      
      #cat("Results of Randomization:\n")
      #print(randomize)
    }
    
    return(invisible(list(fieldbook = randomize, plan = exptLayout)))
    
  } ## end stmt -- designPRep function
