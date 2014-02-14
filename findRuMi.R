# Title: findRuMi.R
# Author: Ian Gonzalez, 2/2014
# Package: SemDist
# This file contains code to access existing information accretion data about the ontology
# and compare the true terms to a given set of predicted terms for a protein using the 
# remaining uncertainty and misinformation metrics (information content analogs to 
# type 1 and type 2 errors).

##----------------------------------------------------------------------------------------##

## function to load information accretion data into the environment.

loadIA <- function(organism,ont) { ##taken from GOSemSim, dont get it yet
	if(!exists("ICEnv")) .initial() 
	fname <- paste("Info_Accretion",
                   organism,
                   ont,
                   sep="_")
    if (ont == "DO") {
        tryCatch(utils::data(list=fname,
                             package="DOSE"))
    } else { ## This method is outdated -- see "data()" help page
        tryCatch(utils::data(list=fname, envir=environment())) 
        ##what env should I load it into? -IG
                             ##package="GOSemSim"))
    }
    IA <- get("IAccr")
    org.ont.IA <- paste(organism,
                        ont,
                        "IA",
                        sep="")
    assign(eval(org.ont.IA),
           IA,
           envir=ICEnv) ## At some point change the name of this Env (maybe IAEnv)
    rm(IA)
}

## function to get information accretion data from the environment

getIA <- function(organism, ont) { ##also lifted from GOSemSim
    if(!exists("ICEnv")) {
        .initial()
    }

    org.ont.IA <- paste(organism,
                        ont,
                        "IA",
                        sep="")

    if(!exists(org.ont.IA, envir=ICEnv)) {
        loadIA(organism, ont)
    }
    IA <- get(org.ont.IA, envir=ICEnv)
    return(IA)
}

## Given a set of predicted IDs and a set of true IDs, computeRU returns the 
## remaining uncertainty for the given ontology using information accretion data. 

computeRU <- function(trueIDs, predIDs, ont, organism, method = "Clark") {
    IA      <- getIA(organism, ont)
    misses  <- setdiff(trueIDs, predIDs)
    RU      <- sum(IA[misses][!is.na(IA[misses])])

    return(RU)
}

## Given a set of predicted IDs and a set of true IDs, computeMI returns the 
## misinformation for the given ontology using information accretion data. 

computeMI <- function(trueIDs, predIDs, ont, organism, method = "Clark") {
    IA      <- getIA(organism, ont)
    wrongs  <- setdiff(predIDs, trueIDs)
    MI      <- sum(IA[wrongs][!is.na(IA[wrongs])])
    
    return(MI)
}

## Adding the threshold values:
## add a set of "scores" from 0:1 that the function predictor must give for each
## term it considered. We step through the possible threshold acceptance values
## from 0 to 1 and graph it to produce the relevant plot.

## Next: 2 functions to read in the predicted and true data from files.

getPredictions <- function(filename) {
    predictions           <- read.table(filename, colClasses = "character")
    colnames(predictions) <- c("seqids","terms","scores")
    predictions$scores <- as.numeric(predictions$scores)    
    predictions <- predictions[predictions$scores != 0,]   #remove predictions with score 0
    return(predictions)
}

getTrues <- function(filename) {
    trues           <- read.table(filename, colClasses = "character")
    colnames(trues) <- c("seqids", "terms")
    return(trues)
}

## Next: A function that the user can access that will return all the desired information.
## find.RU.MI will take a file of predicted terms (1 column each for sequences, terms, and scores from 0-1),
## a threshold value, the relevant ontological info (ont/organism) and return a data frame containing
## RU and MI for each sequence whose terms were predicted.

## In progess: theshcurve=TRUE will create ru-mi curves by stepping through threshold values 0-1.

find.RU.MI <- function(predfile, truefile, threshold = 0.0, ont, organism, threshcurve = FALSE) {   ## NEEDS A BETTER NAME
     
    trueIDs <- getTrues(truefile)       ##Read in data from given files
    predIDs <- getPredictions(predfile)
    seqs    <- unique(predIDs$seqids)   ##Get list of sequences whose annotations have been predicted
    answers <- data.frame("MI" = rep(NA, length(seqs)),     ##Initialize data frame for RU/MI values
                          "RU" = rep(NA, length(seqs)), 
                          row.names = seqs)

    ## If no threshold curve desired, calculate the RU and MI for predicted terms above threshold value
    ## provided and put into the answers data frame. 

    ##SELF NOTE: LOOK INTO WAYS TO MAKE THIS FASTER. MAYBE COMBINE SAPPLYs SOMEHOW? -IG
    if (!threshcurve) {
        #firsterr <- TRUE

        ## Find misinformation for each sequence by calling computeMI on the sets of 
        ## terms related to each sequence:
        answers$MI <- sapply(seqs, function(seq){
            ## If the sequence isn't in the true terms data, return an NA.
            if (!(seq %in% trueIDs$seqids)) {
                if (firsterr) {
                    #cat("WARNING: One or more predicted sequences not found in true term file.\n")
                    firsterr = FALSE
                }
                NA
            } else {
               computeMI(trueIDs$terms[trueIDs$seqids == seq],
                          predIDs$terms[predIDs$seqids == seq & predIDs$scores > threshold],
                          ont,
                          organism)
            }
        })
        
        ## Find remaining uncertainty for each sequence by calling computeRU on the sets of 
        ## terms related to each sequence:
        answers$RU <- sapply(seqs, function(seq){
            ## Same as above.
            if (!(seq %in% trueIDs$seqids)) {
                    if (firsterr) {
                        #cat("WARNING: One or more predicted sequences not found in true term file.\n")
                        firsterr = FALSE
                    }
                    NA
            } else {
                computeRU(trueIDs$terms[trueIDs$seqids == seq],
                      predIDs$terms[predIDs$seqids == seq & predIDs$scores > threshold],
                      ont,
                      organism)
            }
        })          
        return(answers)
    }
    
    else {
        ## code to produce ru-mi curve (unless we want it in another function)
    }
}