##From GOSemSim package

<<<<<<< HEAD
SemDistEnv <- new.env(parent=emptyenv())
IAEnv <- new.env(parent=emptyenv())

.onLoad <- function(libname, pkgname) {
    SemDistEnv[["SupportedSpecies"]] <- c("anopheles",
                                         "arabidopsis",
                                         "bovine",
                                         "canine",
                                         "chicken",
                                         "chimp",
                                         "ecolik12",
                                         "ecsakai",
                                         "fly",
                                         "human",
                                         "malaria",
                                         "mouse",
                                         "pig",
                                         "rat",
                                         "rhesus",
                                         "worm",
                                         "xenopus",
                                         "yeast",
                                         "zebrafish")
    
=======
.initial <- function() {
    #evironment names changed from GOSemSimEnv, ICEnv
    assign("SemDistEnv", new.env(),.GlobalEnv)
    assign("SemSimCache", new.env(), .GlobalEnv)
    assign("IAEnv", new.env(), .GlobalEnv)
    assign("SupportedSpecies", c("anopheles",
                                 "arabidopsis",
                                 "bovine",
                                 "canine",
                                 "chicken",
                                 "chimp",
                                 "ecolik12",
                                 "ecsakai",
                                 "fly",
                                 "human",
                                 "malaria",
                                 "mouse",
                                 "pig",
                                 "rat",
                                 "rhesus",
                                 "worm",
                                 "xenopus",
                                 "yeast",
                                 "zebrafish"),
           envir=SemDistEnv)
    ## remove "coelicolor" as it is not supported by Bioconductor
>>>>>>> fce8d5b26e6ecfb454a764e7aa3c606d22b60638
}

##' @importFrom GO.db GOMFPARENTS
##' @importFrom GO.db GOBPPARENTS
##' @importFrom GO.db GOCCPARENTS
.getParents <- function(ont) {
<<<<<<< HEAD
    Parents.name <- switch(ont,
                      MF = "GOMFPARENTS",
                      BP = "GOBPPARENTS",
                      CC = "GOCCPARENTS",
                      stop("Unsupported ontology: ", sQuote(ont))
                      )
    ontEnv <- suppressPackageStartupMessages(loadNamespace("GO.db"))
    get(Parents.name, envir=ontEnv)
}

.getAncestors <- function(ont) {
    Ancestor.name <- switch(ont,
                      MF = "GOMFANCESTOR",
                      BP = "GOBPANCESTOR",
                      CC = "GOCCANCESTOR",
                      stop("Unsupported ontology: ", sQuote(ont))
                      )
    ontEnv <- suppressPackageStartupMessages(loadNamespace("GO.db"))
    get(Ancestor.name, envir=ontEnv)
}
=======
    Parents <- switch(ont,
                      MF = "GOMFPARENTS",
                      BP = "GOBPPARENTS",
                      CC = "GOCCPARENTS",
                      DO = "DOPARENTS"
                      )
    if (ont == "DO") {
        db <- "DO.db"
        require(db, character.only=TRUE)
    }
    Parents <- eval(parse(text=Parents))
    return(Parents)
}

>>>>>>> fce8d5b26e6ecfb454a764e7aa3c606d22b60638
