# Title: computeIA.R
# Author: Ian Gonzalez, 1/2014
# Package: SemDist
# Contains functions: readOntology, computeIA

# A method to compute the information accretion of all the GO terms 
# for a given ontology.
# Currently, this is a test modification of GOSemSim's computeIC function intended as a framework
# within which to test my method. All changes from the original marked and explained with comments.

##----------------------------------------------------------------------------##

#The readOntology function reads in an ontology file (all parent-child links specified
#in 2 tab delineated columns) and outputs the links, the terms, a term to ancestors object
# and a term to parents object.
readOntology <- function(ontology="mfo_ontology.txt") {
  myGOlinks <- read.csv(ontology, sep="|", colClasses = "character", header=FALSE)
  myGOterms <- unique(append(myGOlinks[,1], myGOlinks[,3]))
  names(myGOlinks) <- c("parents", "pterm", "children", "cterm")
  myGO <- myGOlinks #this saves links so they can be returned
  for (ID in myGOlinks$parents) {
    if (!(ID %in% myGOlinks$children)) {
      root <- ID
      break
    }
  }
  
  Ancestor <- as.list(rep("",length(myGOterms)))
  names(Ancestor) <- myGOterms
  Ancestor <- lapply(Ancestor, function(x) x[x!=""])
  
  Parent <- as.list(rep("",length(myGOterms)))
  names(Parent) <- myGOterms
  Parent <- lapply(Parent, function(x) x[x!=""])
  
  queue <- root
  done <- 0
  
  while (length(queue) > 0) {
    cat("Links done: ", done, "\n")
    todo <- length(queue)
    for (i in 1:todo) {
      if (is.na(queue[1])) {
        cat("missing val detected\n")
        cat("i = ", i, "\n")
        cat("Qlen = ", length(queue), "\n")
        cat("Expected: ", todo, "\n")
      }
      if (queue[1] %in% Ancestor[[queue[1]]]) {
        cat("Error: Term is an ancestor of itself\n")
        break
      }
      termchildren <- myGOlinks$children[myGOlinks$parents == queue[1]]
      myGOlinks <- myGOlinks[myGOlinks$parents != queue[1],]
      done <- done + length(termchildren)
      for (child in termchildren) {
        Parent[[child]] <- append(Parent[[child]], queue[1])
        Ancestor[[child]] <- append(Ancestor[[child]], queue[1])
        Ancestor[[child]] <- append(Ancestor[[child]], Ancestor[[ queue[1] ]])
      }
      queue <- append(queue, termchildren)
      queue <- queue[-1]
    }
  }
  
  if(length(myGOlinks$children) != 0) {
    warning("Some non-root terms did not have parents and were not included.\n")
  }
  return(list(myGO, myGOterms, Ancestor, Parent))
}


#The compute IA function calculates IA values for the specified ontology (obtained
#from built in R data).
computeIA <- function(ont, organism, evcodes="", specify.ont=FALSE, myont=NULL,
                      specify.annotations=FALSE, annotfile=NULL) {
  if (specify.annotations && !specify.ont) {
    cat("Error: Must specify ontology if annotations are being speicified.\n")
    return(0)
  }
## Read in the user's ontology if that option has been specified.
  if (specify.ont) { 
    GOinfo <- readOntology(myont)
    myGO <- GOinfo[[1]]
    myTerms <- GOinfo[[2]]
  }

  ##if the annotations are given, omit all the steps of reading in term/sequence data
  ##and propogating that data
  if (specify.annotations) {  
    myannot <- read.table(annotfile, colClasses = "character")
    names(myannot) <- c("sequences", "GOids")
    goids <- myTerms
    term2seq  <- as.list(rep(character(1),length(goids)))
    names(term2seq) <- goids
    for (ID in goids) {
      term2seq[ID] <- append(term2seq[ID],
                              myannot$sequences[myannot$GOids == ID])
    }



  } else {
  ## Otherwise, get the data from R packages:
  ## Get all the sequence-goid data from GO.db:

    loadGOMap(organism) ##method from gene2GO.R (GOSemSim)
    gomap        <- get("gomap", envir=GOSemSimEnv)
    mapped_genes <- mappedkeys(gomap)
    gomap        <- AnnotationDbi::as.list(gomap[mapped_genes])
    gomap        <- lapply(gomap, function(x){x[!(x$Evidence %in% evcodes)]})
    gomap        <- sapply(gomap, function(x) sapply(x, function(y) y$Ontology))
    
  ## Manipulate the gomap data into a more useful form:

    seq2terms <- sapply(gomap, function(x) {names(x[x==ont])})
    ##seq2terms <- seq2terms[length(seq2terms)==0]
    ##seq2terms <- sapply(gomap, function(x) {names(x[x==ont])})  ##READ THIS: THESE LINES ARE LIKELY REMOVABLE
    
    #If ontology has been specified remove any terms that don't appear in it
    if (specify.ont) {
      seq2terms <- sapply(seq2terms, function(x) { x[x %in% myTerms] })
    }
    
    seq2terms <- seq2terms[sapply(seq2terms, function(x) {if (length(x)==0) {FALSE} else {TRUE}})]
    
    

  ## This ends up being a list that maps sequences to their GO terms in this ont

  ## Now we just need to get ancestors of every GO term to make sure everything is
  ## propagated when we make our table. Either get it from user or from R.
    if (specify.ont) {
      Ancestor <- GOinfo[[3]]
    } else { 
      Ancestor.name <- switch(ont,
                              MF = "GOMFANCESTOR",
                              BP = "GOBPANCESTOR",
                              CC = "GOCCANCESTOR"
      )
      Ancestor <- AnnotationDbi::as.list(get(Ancestor.name,envir=GOSemSimEnv))
      Ancestor <- Ancestor[!is.na(Ancestor)]
    }
    
  ##Propagate annotations:
  ##For each term set in the seq2term object, we append all the ancestors of each term in
  ##the set and apply unique() to remove repeats.
    seq2terms <- lapply(seq2terms,function(terms){
      temp <- unique(c(terms, unlist(sapply(terms, function(term){Ancestor[[term]]}))))
      temp[temp != "all"]
    })

  ## next a list of lists that maps each term to the sequences annotated with it
  ## is created:

  ## Initializing the data structure:
    ## require(GO.db)
    if ( !exists("ALLGOID", envir=GOSemSimEnv) ) {
      assign("ALLGOID", toTable(GOTERM), envir=GOSemSimEnv )
    }
    # Get all the possible terms:
    goids     <- get("ALLGOID", envir=GOSemSimEnv)
    # Make sure they're unique so we don't double-count any:
    goids     <- unique(goids[goids[,"Ontology"] == ont, "go_id"])
    # Initialize the empty list:
    term2seq  <- as.list(rep(character(1),length(goids)))

    names(term2seq) <- goids
    
  ## Loop through sequences and update cell table with each match:    
    for (i in names(seq2terms)) {      
      #Add the sequence to each id in the list that it's annotated with:  
      term2seq[ seq2terms[[i]] ] <- lapply(term2seq[ seq2terms[[i]] ], 
                                            function(x) append(x, i))
    }

    for (i in 1:length(term2seq)) {
      term2seq[[i]] <- term2seq[[i]][ term2seq[[i]] != "" ]
    }
  }

## Now that we have this object, the next step is to calculate IA for each
## term by taking -log of parent count/term count for each term:

  #First, Get all the parent terms for each term.
  if (specify.ont) {
    Parents <- GOinfo[[4]]
  } else { 
    Parents.name <- switch(ont,
                        MF = "GOMFPARENTS",
                        BP = "GOBPPARENTS",
                        CC = "GOCCPARENTS")
    Parents <- AnnotationDbi::as.list(get(Parents.name,envir=GOSemSimEnv))
    Parents <- Parents[!is.na(Parents)]
  }
  
  #Next, a parent count list is created that counts the number of times a
  #term's parent set annotates sequences

  parentcnt <- sapply(names(term2seq), 
                      function (x){ 
                        seqs <- term2seq[[ Parents[[x]][1] ]]
                        for (i in Parents[[x]]) {
                          seqs <- intersect(seqs, term2seq[[i]])
                        }
                        length(seqs)
                      })
  parentcnt <- parentcnt + 1 #add a pseudocount of 1 to each count to prevent 0s
    
  #Now that parent count has been computed, we need only divide this by the 
  #annotation count for each term and take -log2 to get information accretion. 
  
  termcnt <- sapply(term2seq, length)
  #apply a pseudocount of 1 to all terms here as well:
  termcnt <- termcnt + 1
  #sets parent count of root term to #terms so its probability will be 1:
  parentcnt[termcnt==max(termcnt)] <- max(termcnt) 
    
  #Calculate conditional probability, IA for all terms:
  pcond <- termcnt/parentcnt 
  IAccr   <- -log2(pcond)
  save(IAccr, 
       file=paste(paste("Info_Accretion", organism, ont, sep="_"), ".rda", sep=""), 
       compress="xz")
  save(termcnt,
       file=paste(paste("Term_Count", organism, ont, sep="_"), ".rda", sep=""), 
       compress="xz")
  save(parentcnt,
      file=paste(paste("Parent_Count", organism, ont, sep="_"), ".rda", sep=""), 
      compress="xz")
  return(1)

}

##---------------------------------------------------------------### -IG


## General notes on the code:
## Perhaps there is a way to add a functionality that would allow for multiple IA
## data structures to be combined if the user wanted to use multiple (from the same
## section of the ontology).
## Update: this can be done. save the count structures as well as the IA values so that
## multiple species' data can be combined. Then in findRuMi.R, add code to combine a
## specified vector of organisms.
## Also: Add function to read in sequence data set if they want to calculate IA over that.
## Also add user option to read in sequence-term pairs if they so desire
## Possible change: consider combining creation of ancestor and parents object.