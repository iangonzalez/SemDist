##From GOSemSim package

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
}

##' @importFrom GO.db GOMFPARENTS
##' @importFrom GO.db GOBPPARENTS
##' @importFrom GO.db GOCCPARENTS
.getParents <- function(ont) {
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

