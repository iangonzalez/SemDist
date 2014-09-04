##From GOSemSim package

SemDistEnv <- new.env(parent=emptyenv())
IAEnv <- new.env(parent=emptyenv())

##' @author Yu Guangchuang
##' @modified by Ian Gonzalez
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

##' @author Yu Guangchuang
##' @modified by Ian Gonzalez
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

##' @author Yu Guangchuang
##' @modified by Ian Gonzalez
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