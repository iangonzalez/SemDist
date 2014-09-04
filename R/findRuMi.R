# Title: findRuMi.R
# Author: Ian Gonzalez, 2/2014, gonzalez.isv@gmail.com
# Package: SemDist
# Packages required: GO.db
# Contains functions: loadIA, getIA, findRUMI, RUMIcurve, getPredicitons, getTrues

# This file contains code to access existing information accretion data about the ontology
# and compare the true terms to a given set of predicted terms for a protein using the 
# remaining uncertainty and misinformation metrics (information content analogs to 
# type 1 and type 2 errors).

##----------------------------------------------------------------------------------------##

## function to load information accretion data into the environment.
## Accepts single organism or a char vector of them
loadIA <- function(organism, ont) {
  termcnt <- NULL
  parentcnt <- NULL
  if (length(organism) == 1) {
    fname <- paste("Info_Accretion",
                   organism,
                   ont,
                   sep="_")
    tryCatch(data(list=fname, envir=environment()))
    IA <- get("IAccr")
    org.ont.IA <- paste(organism,
                        ont,
                        "IA",
                        sep="")
    assign(eval(org.ont.IA),
           IA,
           envir=IAEnv)
    rm(IA)
    
  } else {
    for (i in 1:length(organism)) {
      fname1 <- paste("Term_Count",
                      organism[i],
                      ont,
                      sep="_")
      fname2 <- paste("Parent_Count",
                      organism[i],
                      ont,
                      sep="_")
      tryCatch(data(list=fname1, envir=environment()))
      tryCatch(data(list=fname2, envir=environment()))
      if (i == 1) {
        pcount <- parentcnt
        tcount <- termcnt
      } else {
        pcount <- pcount + parentcnt
        tcount <- tcount + termcnt
      }
    }
    
    pcond <- tcount/pcount
    IA <- -log2(pcond)
    org.ont.IA <- paste(paste(organism, collapse=""),
                        ont,
                        "IA",
                        sep="")
    assign(eval(org.ont.IA),
           IA,
           envir=IAEnv)
    rm(IA)
  }
}

## function to get information accretion data from the environment
## Accepts single organism or a char vector of them
getIA <- function(organism, ont) {  
  org.ont.IA <- paste(paste(organism, collapse=""),
                      ont,
                      "IA",
                      sep="")
  
  if(!exists(org.ont.IA, envir=IAEnv)) {
    loadIA(organism, ont)
  }
  IA <- get(org.ont.IA, envir=IAEnv)
  IA
}

## 2 functions to read in the predicted and true data from files.

getPredictions <- function(filename) {
  predictions           <- read.table(filename, colClasses = "character")
  colnames(predictions) <- c("seqids","terms","scores")
  predictions$scores <- as.numeric(predictions$scores)    
  predictions <- predictions[predictions$scores != 0,]   #remove predictions with score 0
  predictions <- as.data.frame(predictions)
  predictions
}

getTrues <- function(filename) {
  trues           <- read.table(filename, colClasses = "character")
  trues <- trues[,1:2]
  colnames(trues) <- c("seqids", "terms")
  trues <- as.data.frame(trues)
  trues
}



## findRUMI will take a file of predicted terms (1 column each for sequences, terms, and scores from 0-1),
## a threshold value, the relevant ontological info (ont/organism) and return a data frame containing
## RU and MI for each sequence whose terms were predicted. Alternatively, if fromfile is set to false,
## the function will take the predIDs and trueIDs as data frams with seqids, terms, and scores data columns.
## If seqtrueIAs is set to a vector of values >= 0, the function will treat these as the sums of the true
## terms for each sequence (otherwise it will calculate this on its own).

findRUMI <- function(ont, organism, threshold = 0.05, 
                     truefile="", predfile="", IAccr=NULL) {
  trueIDs <- getTrues(truefile)       ##Read in data from given files
  predIDs <- getPredictions(predfile)
  
  seqs    <- unique(append(predIDs$seqids,trueIDs$seqids))   ##Get list of sequences whose annotations have been predicted
  predIDs <- predIDs[predIDs$scores > threshold,]   ## Remove predictions below the threshold
  if (is.null(IAccr)) {
    IA <- getIA(organism, ont)
  } else {
    IA <- IAccr
  }
  
  ## Get all GO term ancestors for later propagation step:
  Ancestor <- as.list(.getAncestors(ont))
  Ancestor <- Ancestor[!is.na(Ancestor)]
  
  ## For each sequence, find its true and predicted terms, find their intersection,
  ## and find the sum of information accretion over each of these three domains.
  ## RU is calculated by subtracting the intersection from the true IA and MI is
  ## calculated by subtracting the intersection from the pred IA.
  
  seqtrues <- split(trueIDs$terms, f = factor(trueIDs$seqids, seqs))
  
  ## Propagate all the true terms with all ancestors up to the root:
  seqtrues <- lapply(seqtrues,function(terms){
    unique(c(terms, unlist(Ancestor[terms], use.names=FALSE)))
  })
  seqtrues <- lapply(seqtrues, function(x) x[x != "all"])
  
  message("Getting true IAs\n")
  # Helper function to sum IA over a sequence's term set: 
  IA_sum <- function(idx) sum(na.omit(IA[idx]))
  trueIA <- sapply(seqtrues, IA_sum)
  names(trueIA) <- seqs
  
  message("Getting sequence predicted terms.\n")
  seqpreds <- split(predIDs$terms, f = factor(predIDs$seqids, seqs))

  message("Getting IA values for predicted terms.\n")
  predIA <- sapply(seqpreds, IA_sum)
  
  message("Doing the same for the intersect.\n")
  crossover <- Map(intersect, seqtrues, seqpreds)
  crossoverIA <- sapply(crossover, IA_sum)

  message("Calculating RU, MI\n")
  RU <- trueIA - crossoverIA
  MI <- predIA - crossoverIA
  answers <- data.frame(RU=RU,MI=MI)
  
  ## These values are returned in a data frame with RU in first col, MI second.
  answers
}

## RUMIcurve is function that takes in a function predictors predictions, the true annotations,
## and information about which ontology to use and outputs a set of remaining uncertainty/misinformation
## data points that can be used to plot changes in RU/MI as the decision threshold changes.
## Thresholds are tried from 0 to 1 in steps indicated by increment. The add options allow for more
## indicators (weighted RUMI, prec/recall) to be output if needed.

RUMIcurve <- function(ont, organism, increment = 0.05, truefile, predfiles,
                      IAccr = NULL, add.weighted = FALSE, 
                      add.prec.rec = FALSE) {
  thresholds <- seq(1-increment, increment, -1*increment)    ## Create the sequence of thresholds to loop over
  trueIDs <- getTrues(truefile)                           ## Read in data from given files if from file
  
  #Get IA data, from R and getIA function if IAfile is null,
  #from specified .rda file otherwise. Must be called "IA" or "IAccr"
  if (is.null(IAccr)) {
    IA <- getIA(organism, ont)
  } else {
    IA <- IAccr
  }

  output <- list()
  
  ## Get all GO term ancestors for later propagation step:
  Ancestor <- as.list(.getAncestors(ont))
  Ancestor <- Ancestor[!is.na(Ancestor)]
  
  ## For each file given in predfiles:
  ## Get the predicted IDs from the file and generate a list of the IA sum for the true annotations
  ## for each relevant sequence (since this only needs to be computed once).
  for (file in predfiles) {                               
    message("Working on data for file: ", file, "\n")
    predIDs  <- getPredictions(file)
    predIDs  <- predIDs[predIDs$scores != 0,]
    seqs     <- unique(append(predIDs$seqids,trueIDs$seqids))
    seqs     <- sort(seqs)
    
    message("Getting true terms\n")
    seqtrues <- split(trueIDs$terms, f = factor(trueIDs$seqids, seqs))
    
    ## Propagate all the true terms with all ancestors up to the root:
    seqtrues <- lapply(seqtrues,function(terms){
      unique(c(terms, unlist(Ancestor[terms], use.names=FALSE)))
    })
    seqtrues <- lapply(seqtrues, function(x) x[x != "all"])

    # Helper function to sum IA for a sequence's terms
    IA_sum <- function(idx) sum(na.omit(IA[idx]))
    message("Getting true IAs\n")
    trueIA <- sapply(seqtrues, IA_sum)
    names(trueIA) <- seqs
    totalI <- sum(trueIA)
    
    ## Initialize data strctures to hold predicted terms for each sequence and output frame
    seqpreds <- as.list(rep("",length(seqs)))
    seqpreds <- lapply(seqpreds, function(x) x[x != ""])
    names(seqpreds) <- seqs
    
    if (!add.weighted) { 
      answers <- data.frame(threshold = thresholds,
                            RU = rep(0,length(thresholds)),
                            MI = rep(0,length(thresholds)),
                            SS = rep(0,length(thresholds)))
    } else if (!add.prec.rec) {
      answers <- data.frame(threshold = thresholds,
                            RU = rep(0,length(thresholds)),
                            MI = rep(0,length(thresholds)),
                            SS = rep(0,length(thresholds)),
                            WRU = rep(0,length(thresholds)),
                            WMI = rep(0,length(thresholds)),
                            WSS = rep(0,length(thresholds)))
    } else {
      answers <- data.frame(threshold = thresholds,
                            RU = rep(0,length(thresholds)),
                            MI = rep(0,length(thresholds)),
                            SS = rep(0,length(thresholds)),
                            WRU = rep(0,length(thresholds)),
                            WMI = rep(0,length(thresholds)),
                            WSS = rep(0,length(thresholds)),
                            precision = rep(0,length(thresholds)),
                            recall = rep(0,length(thresholds)),
                            specificity = rep(0,length(thresholds)),
                            Wprecision = rep(0,length(thresholds)),
                            Wrecall = rep(0,length(thresholds)))
    }
    
    ## For each threshold, add the predicted terms to the list corresponding to 
    ## their sequence, calc IA, RU, MI, (and other requested values) and put in output frame
    for (thresh in thresholds) {    
      message("Now working on threshold: ", thresh, "\n")
      message("Getting sequence predicted terms.\n")
      
      # Get ONLY the new predicted IDs for this threshold range:
      newpreds <- predIDs[(predIDs$scores > thresh &
                             predIDs$scores <= thresh + increment),]
      
      # if there aren't any, just return the last RUMI values that were calculated
      # or 0 if this is the first thresh being looked at
      if (length(newpreds$seqids) < 1) {
        answers$RU[answers$threshold == thresh] <- mean(RU)
        answers$MI[answers$threshold == thresh] <- mean(MI)
        next
      }
      
      # Then add them to the seqpreds list of lists.
      # This step also propagates all the new terms to the root:
      newpreds <- split(newpreds$terms, f = factor(newpreds$seqids, seqs))
      seqpreds <- Map(union, seqpreds, newpreds)
      seqpreds <- Map(union, seqpreds, lapply(newpreds, function(x) unlist(Ancestor[x], use.names=FALSE)))
      seqpreds <- lapply(seqpreds, function(x) x[x != "all"])
      
      message("Getting IA values for predicted terms.\n")
      predIA <- sapply(seqpreds, IA_sum) 
      
      message("Doing the same for the intersect.\n")
      crossover <- Map(intersect, seqtrues, seqpreds)
      crossoverIA <- sapply(crossover, IA_sum)
      
      RU <- trueIA - crossoverIA
      MI <- predIA - crossoverIA
      SS <- crossoverIA
      message("RU: ", mean(RU), ", MI: ", mean(MI), "\n")
      
      answers$RU[answers$threshold == thresh] <- mean(RU)
      answers$MI[answers$threshold == thresh] <- mean(MI)
      answers$SS[answers$thresholds == thresh] <- mean(SS)
      if (add.weighted) {
        WRU <- sum(RU * trueIA)/totalI ## for wwprecision replace RU
        WMI <- sum(MI * trueIA)/totalI ## with Wprecision
        WSS <- sum(MI * trueIA)/totalI
        
        answers$WRU[answers$threshold == thresh] <- WRU
        answers$WMI[answers$threshold == thresh] <- WMI
        answers$WSS[answers$thresholds == thresh] <- WSS
      }
      if (add.prec.rec) {
        precision <- sum(sapply(crossover, length)) / sum(sapply(seqpreds, length))
        recall <- sum(sapply(crossover, length)) / sum(sapply(seqtrues, length))
        Wprecision <- SS / predIA
        Wrecall <- SS / trueIA
        TN <- sum(sapply(seqtrues, function(x) length(seqs) - length(x)))
        specificity <- TN/(TN + sum(sapply(seqs, function(x) length(seqpreds[[x]]) - length(crossover[[x]]))))
        
        answers$precision[answers$threshold == thresh] <- precision
        answers$recall[answers$threshold == thresh] <- recall
        answers$specificity[answers$thresholds == thresh] <- specificity
        answers$Wprecision[answers$threshold == thresh] <- Wprecision
        answers$Wrecall[answers$thresholds == thresh] <- Wrecall
      }
    }
    output <- append(output, list(answers))
  }
  names(output) <- predfiles
  output
}