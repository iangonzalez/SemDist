\name{findRUMI}
\alias{findRUMI}
\title{
Information accretion based predictor assessment
}

\description{
  Reads in a file containing the true terms annotating a set of sequences and a file
  containing the predicted terms and scores for a set of sequences and outputs a data
  frame containing the remaining uncertainty and misinformation values for the predictions
  made for each sequence.
}

\usage{
findRUMI(ont, organism, threshold = 0.05, truefile="",
          predfile = "", IAccr = NULL)
}

\arguments{
  \item{ont}{Character representation of ontology version to use. One of "CC", "MF", or "BP"
              , corresponding to Cellular Component, Molecular Function, and Biological
              Process.}
  \item{organism}{A character vector indicating which organism(s) annotation data to use.}
  \item{threshold}{Score above which a predicted annotation should be included in the 
                    calculation. Must be a numeric value between 0 and 1, or else findRUMI
                    throws an error.}
  \item{truefile}{Text file from which to read true annotations of sequences. Should be
                  a tab-delineated file with 2 columns: Sequences and GO terms (accessions).}
  \item{predfile}{Text file from which to read predicted annotations of sequences. Should be
                  a tab-delineated file with 3 columns: Sequences, GO terms (accessions), and
                  probability score from 0 to 1 for each prediction.}
    \item{IAccr}{A variable containing a named numeric vector of IA values for all
  the GO terms being used that will be used for calculations instead of R
  packages. This argument is optional.}
}

\value{
  A data frame containing the RU and MI values for each sequence in the file.
}

\author{
Ian Gonzalez and Wyatt Clark
}

\seealso{
\code{\link{computeIA}}
\code{\link{RUMIcurve}}
}

\examples{
# Using test data sets from SemDist, calculate RU and MI:
truefile <- system.file("extdata", "MFO_LABELS_TEST.txt", package="SemDist")
predfile <- system.file("extdata", "MFO_PREDS_TEST.txt", package="SemDist")
rumiTable <- findRUMI("MF", "human", 0.75, truefile, predfile)
avgRU <- mean(rumiTable$RU)
avgMI <- mean(rumiTable$MI)
}
