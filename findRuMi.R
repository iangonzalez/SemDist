# Title: findRuMi.R
# Author: Ian Gonzalez, 2/2014
# Package: SemDist
# Pakcages required: GO.db
# This file contains code to access existing information accretion data about the ontology
# and compare the true terms to a given set of predicted terms for a protein using the 
# remaining uncertainty and misinformation metrics (information content analogs to 
# type 1 and type 2 errors).

##----------------------------------------------------------------------------------------##

## function to load information accretion data into the environment.
setwd("~/SemDist")
source('~/SemDist/gene2GO.R')
source('~/SemDist/utilities.R')



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

## Next: 2 functions to read in the predicted and true data from files.

getPredictions <- function(filename) {
    predictions           <- read.table(filename, colClasses = "character")
    colnames(predictions) <- c("seqids","terms","scores")
    predictions$scores <- as.numeric(predictions$scores)    
    predictions <- predictions[predictions$scores != 0,]   #remove predictions with score 0
    predictions <- as.data.frame(predictions)
    return(predictions)
}
  
getTrues <- function(filename) {
    trues           <- read.table(filename, colClasses = "character")
    trues <- trues[,1:2]
    colnames(trues) <- c("seqids", "terms")
    trues <- as.data.frame(trues)
    return(trues)
}

## Next: A function that the user can access that will return all the desired information.
## find.RU.MI will take a file of predicted terms (1 column each for sequences, terms, and scores from 0-1),
## a threshold value, the relevant ontological info (ont/organism) and return a data frame containing
## RU and MI for each sequence whose terms were predicted. Alternatively, if fromfile is set to false,
## the function will take the predIDs and trueIDs as data frams with seqids, terms, and scores data columns.
## If seqtrueIAs is set to a vector of values >= 0, the function will treat these as the sums of the true
## terms for each sequence (otherwise it will calculate this on its own).

clarkIA <- read.table("MFO_IA.txt",colClasses="character")
clarkIA2 <- as.numeric(clarkIA[,2])
names(clarkIA2) <- clarkIA[,1]
clarkIA <- clarkIA2


find.RU.MI <- function(predIDs="", trueIDs="", ont, organism, 
                        threshold = 0.05, fromfile = TRUE, truefile="", 
                        predfile="", seqtrueIAs = -1,seqtrues = -1) {
    if (fromfile) {
        trueIDs <- getTrues(truefile)       ##Read in data from given files if from file
        predIDs <- getPredictions(predfile)
    }
    seqs    <- unique(append(predIDs$seqids,trueIDs$seqids))   ##Get list of sequences whose annotations have been predicted
    predIDs <- predIDs[predIDs$scores > threshold,]   ## Remove predictions below the threshold
    IA <- clarkIA
    #IA      <- getIA(organism, ont)

    ## For each sequence, find its true and predicted terms, find their intersection,
    ## and find the sum of information accretion over each of these three domains.
    ## RU is calculated by subtracting the intersection from the true IA and MI is
    ## calculated by subtracting the intersection from the pred IA.
    cat("Getting sequence predicted terms.\n")
    seqpreds <- as.list(rep("",length(seqs)))
    names(seqpreds) <- seqs
    for (i in 1:length(predIDs$seqids)) {
      end <- length(seqpreds[[ predIDs$seqids[i] ]])
      seqpreds[[ predIDs$seqids[i] ]][end + 1] <- predIDs$terms[i]
    }

    seqpreds <- lapply(seqpreds, function(x) x[x!=""])
    cat("Getting IA values for predicted terms.\n")
    predIA <- sapply(seqpreds, function(preds) sum(IA[preds][!is.na(IA[preds])]))
    cat("Doing the same for the intersect.\n")
    crossover <- sapply(seqs, function(seq) {
      intersect(seqtrues[[seq]],seqpreds[[seq]])
    })
    cat(length(crossover),"\n")
    crossoverIA <- sapply(crossover,function(int) sum(IA[int][!is.na(IA[int])]))
    cat(length(crossoverIA),"\n")
    cat("Calculating RU, MI\n")
    RU <- seqtrueIAs - crossoverIA
    MI <- predIA - crossoverIA
    answers <- data.frame(RU=RU,MI=MI)

    ## These values are returned in a data frame with RU in first col, MI second.
    return(answers)
}

## RUMIcurve is function that takes in a function predictors predictions, the true annotations,
## and information about which ontology to use and plots a (base) scatterplot that shows the RU/MI
## curve based on incrementing the threshold by the chosen value.

##IDEA: ADD ... SO THEY CAN CHOOSE PLOT OPTION
truefile <- "MFO_LABELS.txt"
predfiles <- "MFO_BLAST.txt"
ont <- "MF"
organism <- "human"
RUMIcurve <- function(predfiles, truefile, ont, organism, increment = 0.05,...) {
    thresholds <- seq(1-increment, increment, -1*increment)    ## Create the sequence of thresholds to loop over
    trueIDs <- getTrues(truefile)                           ## Read in data from given files if from file
    colors <- c("blue","red","green","orange","purple","brown")
    #plot(0, 0, xlab="Remaining Uncertainty",
    #     ylab="Misinformation",xlim=c(0,10), ylim=c(0,15), type="n")  ## Initialize the plotting space
    i <- 1
    IA <- clarkIA
    #IA <- getIA(organism, ont)

    output <- list()


    ## For each file given in predfiles:
    ## Get the predicted IDs from the file and generate a list of the IA sum for the true annotations
    ## for each relevant sequence (since this only needs to be computed once).
    ## Then step through the threshold values, calculating all the RU/MI information each time
    ## and returning a 2 column data frame containing the average values obtained at each threshold.
    ## Finally, plot the data using the base plot system.
    for (file in predfiles) {                               
        cat("Plotting data for file: ", file, "\n")
        predIDs  <- getPredictions(file)
        predIDs  <- predIDs[predIDs$scores != 0,]
        seqs     <- unique(append(predIDs$seqids,trueIDs$seqids))

        cat("Getting true terms\n")
        seqtrues <- as.list(rep("",length(seqs)))
        names(seqtrues) <- seqs
        for (i in 1:length(trueIDs$seqids)) {
          seqtrues[[ trueIDs$seqids[i] ]] <- append(seqtrues[[ trueIDs$seqids[i] ]],
                                                    trueIDs$terms[i]
        }
        seqtrues <- lapply(seqtrues, function(x) x[x!=""])
        cat("Getting true IAs\n")
        trueIA <- sapply(seqtrues, function(trues){
          sum(IA[trues][!is.na(IA[trues])])
        })
        names(trueIA) <- seqs
        ## Get the RU/MI data by looping through thresholds and calculating the mean RU and MI obtained for each value
## TESTING a new way: ------------------------------------------------------------------------------------------------
        seqpreds <- as.list(rep("",length(seqs)))
        seqpreds <- lapply(seqpreds, function(x) x[x != ""])
        names(seqpreds) <- seqs
        predIA <- rep(0,length(seqs))
        names(predIA) <- seqs
        
        for (thresh in thresholds) {
            answers <- data.frame(threshold = thresholds,   #initialize frame to hold answers
                                   RU = rep(0,length(thresholds)),
                                   MI = rep(0,length(thresholds)))
            cat("Now working on threshold: ",thresh,"\n")
            cat("Getting sequence predicted terms.\n")
            # Get ONLY the new predicted IDs for this threshold range:
            newpreds <- predIDs[(predIDs$scores > thresh &
                                predIDs$scores <= thresh + increment),] 
            #Then add them to the seqpreds list of lists:
            for (i in 1:length(newpreds$seqids)) {
              seqpreds[[ newpreds$seqids[i] ]] <- append(seqpreds[[ newpreds$seqids[i] ]],
                                                         newpreds$terms[i])
            }
            
            cat("Seqpreds len: ",length(seqpreds),"\n")
            cat("Getting IA values for predicted terms.\n")
            for (i in 1:length(seqpreds)) {
              predIA[i] <- sum(IA[ seqpreds[[i]] ][!is.na(seqpreds[[i]])])
            }

            cat("Doing the same for the intersect.\n")
            crossover <- sapply(seqs, function(seq) {
              intersect(seqtrues[[seq]],seqpreds[[seq]])
            })

            crossoverIA <- sapply(crossover, function(int) sum(IA[int][!is.na(IA[int])]))

            cat("Calculating RU, MI\n")
            RU <- trueIA - crossoverIA
            MI <- predIA - crossoverIA
            cat("RU: ", mean(RU),"\n")
            cat("MI: ", mean(MI), "\n")
            answers$RU[answers$threshold == thresh] <- mean(RU)
            answers$MI[answers$threshold == thresh] <- mean(MI)
            output <- append(output, answers)
        }
##END testing a new way--------------------------------------------------------------------------------##

        # data <- data.frame(as.numeric(data[1,]),as.numeric(data[2,]))
        # colnames(data) <- c("RU","MI")
        #points(data, type="l", col=colors[i])  ## plot the data
        i <- i + 1
    }
    #legend(0, 6, legend = predfiles, fill = colors)
    output
}
attempt <- RUMIcurve(predfiles,truefile,ont,organism)
save(attempt, file="attempt2.rda")
