\name{computeIA}
\alias{computeIA}
\title{
Compute information accretion for an ontology
}
\description{
Calculates information accretion for each term in the specified ontology using either user
-specified data or the sequence annotations for the organisms specified (note that organism-specific pacakges must be downloaded separately. See "note" section).
}
\usage{
computeIA(ont, organism, evcodes = NULL, specify.ont = FALSE, 
          myont = NULL, specify.annotations = FALSE, 
          annotfile = NULL)
}

\arguments{
   \item{ont}{Character representation of ontology version to use. One of "CC", "MF", or "BP"
              , corresponding to Cellular Component, Molecular Function, and Biological
              Process.}
  \item{organism}{A character vector indicating which organism's annotation data to use.}
<<<<<<< HEAD
  \item{evcodes}{A character vector specifying which evidence codes to use in the ontology
                data. Default NULL value causes all codes to be used.}
=======
  \item{evcodes}{A character vector specifying which evidence codes to ignore in the ontology
                data. Defaults to none.}
>>>>>>> fce8d5b26e6ecfb454a764e7aa3c606d22b60638
  \item{specify.ont}{
  A boolean indicating whether the user wants to specify their own version
                      of the ontology.}
  \item{myont}{
  Character object indicating what file to read in the specified ontology from. 
<<<<<<< HEAD
  The ontology should be specified as a tab-delimited file with 2 columns (no header). 
  Each row in the file should indicate a parent-child relationship between two GO accessions (e.g. "GO:0003674   GO:0004000")}
=======
  The ontology should be specified as a pipe-delimited (| symbol) file with 4
               columns. The first column contains all parent term accessions and the third column
               contains all child term accessions. Each row should represent a parent-child
               link. Every link in the ontology should be specified. The 2nd and 4th columns
               should contain the english terms corresponding to the accessions in the 1st and
               3rd columns.}
>>>>>>> fce8d5b26e6ecfb454a764e7aa3c606d22b60638
  \item{specify.annotations}{Boolean indicating whether the user wants to specify sequence
                          annotations from a file. Should only be TRUE if specify.ont is TRUE.}
  \item{annotfile}{Character object indicating which file to read sequence annotations from.
                    Should be a tab-delimited file with 2 columns. The first column is a
                    list of sequences, the second is a list of GO accessions in the same rows
                    as the sequences they annotate.}
}
\value{
Does not return a specific value. Saves the information accretion values for each term in the ontology
in a .rda file that specifies the organism and the ont version. Parent count and term count
objects are also saved in similarly formatted files so that IA calculations from multiple
organisms can be combined.
}

\note{
In order to compute IA for an organism, the specific annotation data set for that organism must be installed by the user. Here is a list of supported organisms (names in the format that should be passed to computeIA) and the correspoinding packages needed:

anopheles   = org.Ag.eg.db

arabidopsis = org.At.tair.db

bovine      = org.Bt.eg.db

canine      = org.Cf.eg.db

chicken     = org.Gg.eg.db

chimp       = org.Pt.eg.db

ecolik12    = org.EcK12.eg.db

fly         = org.Dm.eg.db

human       = org.Hs.eg.db

malaria     = org.Pf.plasmo.db

mouse       = org.Mm.eg.db

pig         = org.Ss.eg.db

rat         = org.Rn.eg.db

rhesus      = org.Mmu.eg.db

worm        = org.Ce.eg.db

xenopus     = org.Xl.eg.db

yeast       = org.Sc.sgd.db

zebrafish   = org.Dr.eg.db
}

\author{
Ian Gonzalez and Wyatt Clark
}

\seealso{
\code{\link{RUMIcurve}}
\code{\link{findRUMI}}
}

\examples{
<<<<<<< HEAD
          
# Calculate IA, specify ontology and annotations
ontfile <- system.file("extdata", "mfo_ontology.txt", package="SemDist")
=======
\dontrun{
# Calculate IA, specify ontology
ontfile <- system.file("extdata", "mfo_ontology.txt", package="SemDist")
computeIA("MF", "human", specify.ont=TRUE,
         myont=ontfile, specify.annotations=FALSE)
          
# Calculate IA, specify ontology and annotations
>>>>>>> fce8d5b26e6ecfb454a764e7aa3c606d22b60638
annotations <- system.file("extdata", "MFO_LABELS_TEST.txt", package="SemDist")
computeIA("my", "values", specify.ont=TRUE, 
          myont=ontfile, specify.annotations=TRUE, 
          annotfile=annotations)
<<<<<<< HEAD
=======
  }
>>>>>>> fce8d5b26e6ecfb454a764e7aa3c606d22b60638
}
