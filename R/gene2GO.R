##From GOSemSim package

getDb <- function(organism) {
    annoDb <- switch(organism,
                     anopheles   = "org.Ag.eg.db",
                     arabidopsis = "org.At.tair.db",
                     bovine      = "org.Bt.eg.db",
                     canine      = "org.Cf.eg.db",
                     chicken     = "org.Gg.eg.db",
                     chimp       = "org.Pt.eg.db",
                     ecolik12    = "org.EcK12.eg.db",
                     ecsakai     = "org.EcSakai.eg.db",
                     fly         = "org.Dm.eg.db",
                     human       = "org.Hs.eg.db",
                     malaria     = "org.Pf.plasmo.db",
                     mouse       = "org.Mm.eg.db",
                     pig         = "org.Ss.eg.db",
                     rat         = "org.Rn.eg.db",
                     rhesus      = "org.Mmu.eg.db",
                     worm        = "org.Ce.eg.db",
                     xenopus     = "org.Xl.eg.db",
                     yeast       = "org.Sc.sgd.db",
                     zebrafish   = "org.Dr.eg.db",
                     stop("Unsupported organism: ", sQuote(organism))
                     )
    annoDb
}

loadGOMap_internal <- function(organism){
    annoDb <- getDb(organism)

    ## loading annotation package
    annoEnv <- suppressPackageStartupMessages(loadNamespace(annoDb))

    gomap <- switch(organism,
                    anopheles    = "org.Ag.egGO",
                    arabidopsis  = "org.At.tairGO",
                    bovine       = "org.Bt.egGO",
                    canine       = "org.Cf.egGO",
                    chicken      = "org.Gg.egGO",
                    chimp        = "org.Pt.egGO",
                    ecolik12     = "org.EcK12.egGO",
                    ecsakai      = "org.EcSakai.egGO",
                    fly          = "org.Dm.egGO",
                    human        = "org.Hs.egGO",
                    malaria      = "org.Pf.plasmoGO",
                    mouse        = "org.Mm.egGO",
                    pig          = "org.Ss.egGO",
                    rat          = "org.Rn.egGO",
                    rhesus       = "org.Mmu.egGO",
                    worm         = "org.Ce.egGO",
                    xenopus      = "org.Xl.egGO",
                    yeast        = "org.Sc.sgdGO",
                    zebrafish    = "org.Dr.egGO",
                    stop("Unsupported organism: ", sQuote(organism))
                    )
    SemDistEnv[[annoDb]] <- annoEnv[[annoDb]]
    SemDistEnv[["gomap"]] <- get(gomap, annoEnv)
    SemDistEnv[["gomap.flag"]] <- organism
}

##' @importMethodsFrom AnnotationDbi exists
##' @importMethodsFrom AnnotationDbi get
loadGOMap <- function(organism) {
    annoDb <- getDb(organism)
    if (!exists(annoDb, envir=SemDistEnv) || !exists("gomap", envir=SemDistEnv)) {
        loadGOMap_internal(organism)
    } else {
        flag <- get("gomap.flag", envir=SemDistEnv)
        if (flag != organism)
            loadGOMap_internal(organism)
    }
}

##' @importMethodsFrom AnnotationDbi get
gene2GO <- function(gene, organism, ont, dropCodes) {
    gene <- as.character(gene)
    loadGOMap(organism)
    gomap <- get("gomap", envir=SemDistEnv)
    go <- gomap[[gene]]

    if (all(is.na(go)))
        return (NA)

    goid <- sapply(go, function(i) i$GOID)
    evidence <- sapply(go, function(i) i$Evidence)
    ontology <- sapply(go, function(i) i$Ontology)

    idx <- (evidence %in% dropCodes) & (ontology == ont)
    if (sum(idx) == 0)
        NA
    else
        unname(goid[idx])
}
