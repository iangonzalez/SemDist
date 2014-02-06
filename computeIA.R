# Author: Ian Gonzalez
# Package: SemDist
# A method to compute the information accretion of all the GO terms 
# for a given ontology.
# Currently, this is a test modification of GOSemSim's computeIC function intended as a framework
# within which to test my method. All changes from the original marked and explained with comments.

## IG: .getParents is a modified .getOffsprings that obtains parent data instead. 
.getAncestor <- function(ont="MF") {
    if(!exists("GOSemSimEnv")) .initial()
    wh_Ancestor <- switch(ont,
                            MF = "MFANCESTOR",
                            BP = "BPANCESTOR",
                            CC = "CCANCESTOR"
                            )
    Ancestor <- switch(ont,
                         MF = AnnotationDbi::as.list(GOMFANCESTOR) ,
                         BP = AnnotationDbi::as.list(GOBPANCESTOR) ,
                         CC = AnnotationDbi::as.list(GOCCANCESTOR)
                         )
    assign(eval(wh_Ancestor), Ancestor, envir=GOSemSimEnv)#environment will be changed before package is complete.
}
##----------------------------------------------------------------------------##

computeIA <- function(organism,ont){
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
    term2seq<-as.list(rep(character(1),length(goids)))
    names(term2seq)<-goids
    
## Loop through sequences and update cell table with each match:
## (this part will be computationally intensive)
    for (i in names(test)) {
      idlist<-test[[i]] #get the ids pertaining to each sequence
      
      #Add the ancestor terms as well (propagates all annotations):
      idlist<-c(idlist, 
                unlist(sapply(test[[i]], function(x){Ancestor[[x]]})))
      #Remove repeats and the meaningless "all" tag:
      idlist<-unique(idlist)
      idlist<-idlist[idlist != "all"]
      
      #Add the sequence to each id in the list:
      for (j in idlist) {
        term2seq[[j]]<-append(term2seq[[j]],i)
      }
    }


##---------------------------------------------------------------### -IG


    ## all GO terms appearing in an given ontology ###########
    goterms <- unlist(sapply(gomap, function(x) names(x[x == ont])), use.names=FALSE)

    ## require(GO.db)
    if ( !exists("ALLGOID", envir=GOSemSimEnv) ) {
        assign("ALLGOID", toTable(GOTERM), envir=GOSemSimEnv )
    }
    #Get calls an R object using a character string.-IG
    goids   <- get("ALLGOID", envir=GOSemSimEnv)
    ##goids <- toTable(GOTERM)

    ## all go terms which belong to the corresponding ontology..
    goids   <- unique(goids[goids[,"Ontology"] == ont, "go_id"])
    gocount <- table(goterms)
    goname  <- names(gocount) #goid of specific organism and selected category.

    ## ensure goterms not appearing in the specific annotation have 0 frequency..
    go.diff        <- setdiff(goids, goname)
    m              <- double(length(go.diff))
    names(m)       <- go.diff
    gocount        <- as.vector(gocount)
    names(gocount) <- goname
    gocount        <- c(gocount, m)
    ##^^^^????? -IG

    Offsprings.name <- switch(ont,
                              MF = "MFOffsprings",
                              BP = "BPOffsprings",
                              CC = "CCOffsprings"
                              )
    if (!exists(Offsprings.name, envir=GOSemSimEnv)) {
        .getOffsprings(ont)
    }
    Offsprings <- get(Offsprings.name, envir=GOSemSimEnv)
    cnt        <- sapply(goids,function(x){ n=gocount[ Offsprings[[x]] ]; gocount[x]+sum(n[!is.na(n)])})
    names(cnt) <- goids
  	## the probabilities of occurrence of GO terms in a specific corpus.
  	p          <- cnt/sum(gocount)
  	## IC of GO terms was quantified as the negative log likelihood.
  	IC         <- -log(p)

##--------------------------------------------------------------------------------------------------##
    ## IG: Added this section to read in Parent data.
    Ancestor.name <- switch(ont,
                              MF = "MFANCESTOR",
                              BP = "BPANCESTOR",
                              CC = "CCANCESTOR"
                              )
    if (!exists(Ancestor.name, envir=GOSemSimEnv)) {#environment will be changed when package complete
        .getAncestor(ont="BP")
    }
    
    squeeze<-character(length=length(Parents))
    for (i in 1:length(Parents)){##collapse all parents into one string for easier use
        squeeze[i]<-paste(Parents[[i]],sep="",collapse="")
    }
    #Calculate how many times each parent set appears (WORK IN PROGRESS)
    #Address the case in which parent sets are equivalent but differently ordered.
    parentcnt<-sapply(goids,function(x){length(squeeze[squeeze==squeeze[x]])}) 
    pcond   <- parentscnt/cnt
    IA      <- -log(pcond)
##---------------------------------------------------------------------------------------------------##

    save(IC, file=paste(paste("Info_Contents", organism, ont, sep="_"), ".rda", sep=""), compress="xz")
}

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