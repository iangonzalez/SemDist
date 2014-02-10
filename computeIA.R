# Author: Ian Gonzalez
# Package: SemDist
# A method to compute the information accretion of all the GO terms 
# for a given ontology.
# Currently, this is a test modification of GOSemSim's computeIC function intended as a framework
# within which to test my method. All changes from the original marked and explained with comments.

##----------------------------------------------------------------------------##

computeIA <- function(organism,ont) {
	  loadGOMap(organism)
    gomap   <- get("gomap", envir=GOSemSimEnv)
    #Understanding this is crucial.Here the code grabs all the gene
    #information for the organism then converts it into a new format. -IG
    mapped_genes <- mappedkeys(gomap)
    gomap <- AnnotationDbi::as.list(gomap[mapped_genes])
    ## GOMAP IS THE SEQUENCES
    gomap <- sapply(gomap, function(x) sapply(x, function(y) y$Ontology))
    
## Manipulate the gomap data into a more useful form:
##
	  test<-sapply(gomap, function(x) {names(x[x==ont])})
	  test<-test[length(test)==0]
	  test<-sapply(gomap, function(x) {names(x[x==ont])})
	  test<-test[sapply(test, function(x) {if (length(x)==0) {FALSE} else {TRUE}})]
## This ends up being a list that maps sequences to their GO terms in this ont

## Now we just need to get the Ancestor list to make sure everything is
## propagated when we make our table.
    Ancestor.name <- switch(ont,
                        MF = "GOMFANCESTOR",
                        BP = "GOBPANCESTOR",
                        CC = "GOCCANCESTOR"
    )
    Ancestor <- AnnotationDbi::as.list(get(Ancestor.name,envir=GOSemSimEnv))
    Ancestor <- Ancestor[!is.na(Ancestor)]
## This Ancestor object contains a mapping of each GO term to ALL its
## ancestors in the given ontology.

## Finally, we step through each sequence in the gomap (test object for now),
## step through each ID for the sequence and all of its ancestors, and add the
## sequence to the list mapped to each ID in a 2D array:

## Initializing the data structure:
    ## require(GO.db)
    if ( !exists("ALLGOID", envir=GOSemSimEnv) ) {
      assign("ALLGOID", toTable(GOTERM), envir=GOSemSimEnv )
    }
    # Get all the possible terms:
    goids   <- get("ALLGOID", envir=GOSemSimEnv)
    # Make sure they're unique so we don't double-count any:
    goids   <- unique(goids[goids[,"Ontology"] == ont, "go_id"])
    # Initialize the empty list: (CONSIDER OTHER IMPLEMENTATIONS OF THE DATA)
    term2seq <- as.list(rep(character(1),length(goids)))
    names(term2seq)<-goids

##Propagate annotations:
    for (i in names(test)) {
      idlist<-test[[i]] #get the ids pertaining to each sequence
      
      #Add the ancestor terms as well (propagates all annotations):
      idlist<-c(idlist, 
                unlist(sapply(test[[i]], function(x){Ancestor[[x]]})))
      #Remove repeats and the meaningless "all" tag:
      idlist<-unique(idlist)
      idlist<-idlist[idlist != "all"]
      test[[i]] <- idlist
    }
    
## Loop through sequences and update cell table with each match:
## (this part will be computationally intensive)
  
    for (i in names(test)) {      
      #Add the sequence to each id in the list:  
      term2seq[ test[[i]] ] <- lapply(term2seq[ test[[i]] ], 
                                      function(x) append(x, i))
    }

  #this removes the empty string from each element. It's kind of a hack, 
  #the initialization of the object should be modified so this isn't needed.
    for (i in 1:length(term2seq)) {
      term2seq[[i]] <- term2seq[[i]][ term2seq[[i]]!="" ]
    }

## Now that we have this object, the next step is to calculate IA for each
## term. To get the parent ocurrence count for each term, look at its parent
## set and find the size of the intersection between their annotated seqs.
## Then compute IA by taking -log of parent count/term count for each term.

  #First, Get all the parent terms for each term.
    Parents.name <- switch(ont,
                        MF = "GOMFPARENTS",
                        BP = "GOBPPARENTS",
                        CC = "GOCCPARENTS")

    Parents <- AnnotationDbi::as.list(get(Parents.name,envir=GOSemSimEnv))
    Parents <- Parents[!is.na(Parents)]
  
  #Next, create the parent count list with a sapply that steps through
  #term2seq's name object and finds the number of times the parents of each
  #term appear together by taking the set intersect of the sets of sequences
  #that each parent annotates.

    parentcnt <- sapply(names(term2seq), 
                        #maybe make this a named function? this code is ugly
                        function (x){ 
                          seqs <- term2seq[[ Parents[[x]][1] ]]
                          for (i in Parents[[x]]) {
                            seqs <- intersect(seqs, term2seq[[i]])
                          }
                          length(seqs)
                        })

  #Now that parent count has been computed, we need only divide this by the 
  #count for each term and take -log2 to get information accretion. 
  
    termcnt <- sapply(term2seq, function (x) {length(x)})
    #apply a pseudocount of 1 to terms that never appear:
    termcnt <- sapply(termcnt, function (x) {if (x == 0) {1} else {x}})

    pcond   <- parentcnt/termcnt
    IAccr   <- -log2(pcond)

    save(IAccr, 
         file=paste(paste("Info_Accretion", organism, ont, sep="_"), ".rda", sep=""), 
         compress="xz")
}

##---------------------------------------------------------------### -IG


## General notes on the code:

## Use set functions: Intersect of A and B should equal length of A
## CREATE A LABEL MATRIX USING GO DATA
## SEE HIS MATLAB CODE
## Creating a label matrix:
## goids = the relevant "sequences"
## cnt properly gives how many times a term appears
## What we need to do is:
## 1. Make a 2D matrix with all goterms as the column terms
## 2. Go through each goid, append a new row with goid as rowname and 1s for each term describing it (0s elsewhere).

## BETTER WAY:
## For each sequence, find annotations, find ancestors of those annotations,
## go through list of all terms, add the sequence name to the list of
## sequences corresponding to that term.
## Maybe use R hash package. Consider just listing the indexes of sequences
## rather than the sequence names themselves