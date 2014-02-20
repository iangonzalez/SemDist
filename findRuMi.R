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
## RU and MI for each sequence whose terms were predicted. Alternatively, if fromfile is set to false,
## the function will take the predIDs and trueIDs in the correct format if they've already been read in.

find.RU.MI <- function(predIDs="", trueIDs="", ont, organism, 
                        threshold = 0.0, fromfile = TRUE, truefile="", predfile="") {
    if (fromfile) {
        trueIDs <- getTrues(truefile)       ##Read in data from given files if from file
        predIDs <- getPredictions(predfile)
    }
    seqs    <- unique(predIDs$seqids)   ##Get list of sequences whose annotations have been predicted
    answers <- data.frame("MI" = rep(NA, length(seqs)),     ##Initialize data frame for RU/MI values
                          "RU" = rep(NA, length(seqs)), 
                          row.names = seqs)
    predIDs <- predIDs[predIDs$scores > threshold,]

    ## If no threshold curve desired, calculate the RU and MI for predicted terms above threshold value
    ## provided and put into the answers data frame. 

    ##SELF NOTE: LOOK INTO WAYS TO MAKE THIS FASTER. MAYBE COMBINE SAPPLYs SOMEHOW? -IG
        #firsterr <- TRUE

    ## Find misinformation for each sequence by calling computeMI on the sets of 
    ## terms related to each sequence:
    answers$MI <- sapply(seqs, function(seq){
        ## If the sequence isn't in the true terms data, return an NA.
        if (!(seq %in% trueIDs$seqids)) {
            #if (firsterr) {
            #cat("WARNING: One or more predicted sequences not found in true term file.\n")
            #    firsterr = FALSE
            #}
            NA
        } else {
            computeMI(trueIDs$terms[trueIDs$seqids == seq],
                       predIDs$terms[predIDs$seqids == seq],
                       ont,
                       organism)
        }
    })
        
    ## Find remaining uncertainty for each sequence by calling computeRU on the sets of 
    ## terms related to each sequence:
    answers$RU <- sapply(seqs, function(seq){
        ## Same as above.
        if (!(seq %in% trueIDs$seqids)) {
            #if (firsterr) {
            #    #cat("WARNING: One or more predicted sequences not found in true term file.\n")
            #    firsterr = FALSE
            #}
            NA
        } else {
            computeRU(trueIDs$terms[trueIDs$seqids == seq],
                      predIDs$terms[predIDs$seqids == seq],
                      ont,
                      organism)
        }
    })          
    return(answers)
}

## RUMIcurve is function that takes in a function predictors predictions, the true annotations,
## and information about which ontology to use and plots a (base) scatterplot that shows the RU/MI
## curve based on incrementing the threshold by the chosen value.

RUMIcurve <- function(predfiles, truefile, ont, organism, increment = 0.05) {
    thresholds <- seq(increment, 1-increment, increment)    ## Create the sequence of thresholds to loop over
    trueIDs <- getTrues(truefile)                           ## Read in data from given files if from file
    colors <- c("blue","red","green","orange","purple","brown")
    plot(0, 0, xlim=c(0,10), ylim=c(0,15), type="n")
    i <- 1
    for (file in predfiles) {                               ## For each file given in predfiles:
        cat("Plotting data for file: ", file, "\n")
        predIDs <- getPredictions(file)
        predIDs <- predIDs[predIDs$scores != 0,]
        ## Get the RU/MI data by looping through thresholds and calculating the mean RU and MI obtained for each value
        data    <- sapply(thresholds, function(thresh) {
            threshdata <- find.RU.MI(predIDs, trueIDs, ont, organism, threshold = thresh, fromfile = FALSE)
            c(mean(threshdata$RU[!is.na(threshdata$RU)]), mean(threshdata$MI[!is.na(threshdata$MI)]))
            })

        ## Manipulate the data into ggplot-readable form and add its points and fit line onto the plot:
        data <- data.frame(as.numeric(data[1,]),as.numeric(data[2,]))
        colnames(data) <- c("RU","MI")
        points(data, type="l", col=colors[i])
        i <- i + 1
    }
    legend(0,4,legend = predfiles, fill = colors)
}
