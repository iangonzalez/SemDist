\name{computeIA}
\alias{computeIA}
\title{
Compute information accretion for an ontology
}
\description{
Calculates information accretion for each term in the specified ontology using either user
-specified data or the sequence annotations for the organisms specified (data obtained from
annotate package).
}
\usage{
computeIA(ont, organism, evcodes = "", specify.ont = FALSE, 
          myont = NULL, specify.annotations = FALSE, 
          annotfile = NULL)
}

\arguments{
   \item{ont}{Character representation of ontology version to use. One of "CC", "MF", or "BP"
              , corresponding to Cellular Component, Molecular Function, and Biological
              Process.}
  \item{organism}{A character vector indicating which organism's annotation data to use.}
  \item{evcodes}{A character vector specifying which evidence codes to ignore in the ontology
                data. Defaults to none.}
  \item{specify.ont}{A boolean indicating whether the user wants to specify their own version
                      of the ontology.}
  \item{myont}{Character object indicating what file to read in the specified ontology from.
               The ontology should be specified as a pipe-delimited (| symbol) file with 4
               columns. The first column contains all parent term accessions and the third column
               contains all child term accessions. Each row should represent a parent-child
               link. Every link in the ontology should be specified. The 2nd and 4th columns
               should contain the english terms corresponding to the accessions in the 1st and
               3rd columns.}
  \item{specify.annotations}{Boolean indicating whether the user wants to specify sequence
                          annotations from a file. Should only be TRUE if specify.ont is TRUE.}
  \item{annotfile}{Character object indicating which file to read sequence annotations from.
                    Should be a tab-delimited file with 2 columns. The first column is a
                    list of sequences, the second is a list of GO accessions in the same rows
                    as the sequences they annotate.}
}
\value{
Returns a 1 if successful. Saves the information accretion values for each term in the ontology
in a .rda file that specifies the organism and the ont version. Parent count and term count
objects are also saved in similarly formatted files so that IA calculations from multiple
organisms can be combined.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ian Gonzalez and Wyatt Clark
}

\seealso{
\code{\link{RUMIcurve}}
\code{\link{findRUMI}}
\code{\link{SemDist}}
}
