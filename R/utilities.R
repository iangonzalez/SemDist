##From GOSemSim package

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
    
}

##' @importFrom GO.db GOMFPARENTS
##' @importFrom GO.db GOBPPARENTS
##' @importFrom GO.db GOCCPARENTS
.getParents <- function(ont) {
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