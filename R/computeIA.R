# Title: computeIA.R
# Author: Ian Gonzalez, 1/2014, gonzalez.isv@gmail.com
# Package: SemDist
# Contains functions: readOntology, computeIA

# A method to compute the information accretion of all the GO terms 
# for a given ontology.

##----------------------------------------------------------------------------##

#The readOntology function reads in an ontology file (all parent-child links specified
#in 2 tab delineated columns) and outputs the links, the terms, a term to ancestors object
# and a term to parents object.
readOntology <- function(ontology="mfo_ontology.txt") {
  myGOlinks <- read.table(ontology, colClasses = c("character","character"), header=FALSE)
  myGOterms <- unique(append(myGOlinks[,1], myGOlinks[,2]))
  names(myGOlinks) <- c("parents", "children")
  myGO <- myGOlinks #this saves links so they can be returned

  # Create parent and ancestor objects:
  Parents <- split(myGOlinks$parents, f = myGOlinks$children)

  #Ancestor includes all terms above a term in ontology:
  ## Build a list of each parent generation from the previous generation,
  ## append to ancestor object.
  map <- Ancestor <- Parents
  while (length(map)) {
      child <- rep(names(map), sapply(map, length))
      map <- Parents[unlist(map, use.names=FALSE)]
      child <- rep(child, sapply(map, length))
      if (all(sapply(map, is.null))) {
        break
      }
      map <- split(unlist(map, use.names=FALSE), child)
      idx <- names(map)
      Ancestor[idx] <- Map(append, Ancestor[idx], map)
  }
  Ancestor <- lapply(Ancestor, unique)

  list(myGO, myGOterms, Ancestor, Parents)
}

#The compute IA function calculates IA values for the specified ontology (obtained
#from built in R data).
computeIA <- function(ont, organism, evcodes=NULL, specify.ont=FALSE, myont=NULL,
                      specify.annotations=FALSE, annotfile=NULL) {
  
  if (specify.annotations && !specify.ont) {
    stop("Error: Must specify ontology if annotations are being specified.\n")
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
    myannot <- read.table(annotfile, colClasses = c("character","character"))
    names(myannot) <- c("sequences", "GOids")
    goids <- myTerms

    term2seq <- split(myannot$sequences, f = myannot$GOids)
    
  } else {
    ## Otherwise, get the data from R packages:
    ## Get all the sequence-goid data from GO.db:
    loadGOMap(organism) ##method from gene2GO.R (GOSemSim)
    db.name <- getDb(organism)
    annoDb <- get(db.name, envir=SemDistEnv)
    
    # Get the protein IDs appropriate to each package:
    if (grepl(".eg.", db.name, fixed=TRUE)) {
      identifier <- "ENTREZID"
    } else if (grepl(".sgd.", db.name, fixed=TRUE)) {
      identifier <- "SGD"
    } else if (grepl(".plasmo.", db.name, fixed=TRUE)) {
      identifier <- "GENENAME"
    } else if (grepl(".tair.", db.name, fixed=TRUE)) {
      identifier <- "TAIR"
    } else stop("Unsupported annotation package: ", sQuote(db.name))
    
    seq2terms <- suppressWarnings(select(annoDb, keys(annoDb, keytype=identifier), "GO", identifier))
    seq2terms <- seq2terms[seq2terms$ONTOLOGY == ont,]
    seq2terms <- seq2terms[!is.na(seq2terms$GO),]

    # If the evcodes arg was specified, any annotations not in those
    # evidence codes are removed from the list
    if (!is.null(evcodes)) {
      if (length(seq2terms$EVIDENCE[seq2terms$EVIDENCE %in% evcodes]) == 0) {
        stop("No annotations match the evidence code specifications for organism ", sQuote(organism))
      } else {
        seq2terms <- seq2terms[seq2terms$EVIDENCE %in% evcodes,]
      }
    }

    seq2terms <- seq2terms[, 1:2]
    
    #If ontology has been specified remove any terms that don't appear in it
    if (specify.ont) {
      seq2terms <- seq2terms[seq2terms$GO %in% myTerms,]
    }

    ## This ends up being a data frame that maps sequences to their GO terms in this ont
    
    ## Now we just need to get ancestors of every GO term to make sure everything is
    ## propagated when we make our table. Either get it from user or from R.
    if (specify.ont) {
      Ancestor <- GOinfo[[3]]
    } else { 
      Ancestor <- as.list(.getAncestors(ont))
      Ancestor <- Ancestor[!is.na(Ancestor)]
    }
    
    ##Propagate annotations:
    ##For each term set in the seq2term object, we add entries for all the ancestors of 
    ##each term in the set.

    seq2terms <- split(seq2terms$GO, f = seq2terms[,1])
    seq2terms <- lapply(seq2terms, function(terms){
      temp <- unique(c(terms, unlist(sapply(terms, function(term){Ancestor[[term]]}))))
      temp[temp != "all"]
    })
    new_id_vec <- unlist(lapply(names(seq2terms), function(x) rep(x, length(seq2terms[[x]]))))
    seq2terms <- data.frame(IDS=new_id_vec, GO=unlist(seq2terms), stringsAsFactors=FALSE)

    ## next a list of lists that maps each term to the sequences annotated with it
    ## is created:
    
    term2seq <- split(seq2terms$IDS, f = seq2terms$GO)
  }
  
  ## Now that we have this object, the next step is to calculate IA for each
  ## term by taking -log of parent count/term count for each term:
  
  #First, Get all the parent terms for each term.
  if (specify.ont) {
    Parents <- GOinfo[[4]]
  } else { 
    Parents <- as.list(.getParents(ont))
    Parents <- Parents[!is.na(Parents)]
  }
  
  #Next, a parent count list is created that counts the number of times a
  #term's parent set annotates sequences
  
  parentcnt <- sapply(names(term2seq), 
                      function (x){ 
                        if (!(x %in% names(Parents))) {
                          return(0)
                        }
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
  parentcnt[termcnt == max(termcnt)] <- max(termcnt) 
  
  #Calculate conditional probability, IA for all terms:
  pcond <- termcnt/parentcnt 
  IAccr   <- -log2(pcond)
  IAccr  
}

##---------------------------------------------------------------### -IG