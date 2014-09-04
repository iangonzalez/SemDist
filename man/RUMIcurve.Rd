\name{RUMIcurve}
\alias{RUMIcurve}
\title{
Information accretion based predictor assessment (across many thresholds)
}
\description{
Reads in a (tab-delimited) file containing the true annotations for a set of sequences, a
(tab-delimited) file containing the predicted annotations and corresponding scores for the
same sequences. Calculates and outputs the average remaining uncertainty, misinformation,
and semantic similarity at a series of user-specified thresholds.
}
\usage{
RUMIcurve(ont, organism, increment = 0.05, truefile, predfiles, 
<<<<<<< HEAD
          IAccr = NULL, add.weighted = FALSE, 
          add.prec.rec = FALSE)
=======
          IAfile = NULL, outfile = "rumicurve.rda", 
          add.weighted = FALSE, add.prec.rec = FALSE)
>>>>>>> fce8d5b26e6ecfb454a764e7aa3c606d22b60638
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ont}{Character representation of ontology version to use. One of "CC", "MF", or "BP"
              , corresponding to Cellular Component, Molecular Function, and Biological
              Process.}
  \item{organism}{A character vector indicating which organism(s) annotation data to use.}
  \item{increment}{A numeric value between 0 and 1 indicating the distance between each
                    threshold that should be calculated. Note that the iteration starts from
                    a threshold of 1, so an increment value of 0.08 will result in the
                    thresholds 0.92, 0.84, 0.76 ... being used.}
  \item{truefile}{A character vector indicating the file from which to read the true
                  annotations for the given sequences. Should be tab-delimited, with the
                  first column containing the sequence ids and the second containing GO
                  accessions.}
  \item{predfiles}{A character vector containing which files to read in as the predicted
                    annotations. Should be tab-delimited, with the first column containing
                    sequences, the second column containing GO accessions, and the third
                    column containing the predictors 0-1 score for that prediction.}
<<<<<<< HEAD
  \item{IAccr}{A variable containing a named numeric vector of IA values for all
  the GO terms being used that will be used for calculations instead of R
  packages. This argument is optional.}
=======
  \item{IAfile}{A .rda filename (character) containing a named numeric vector of IA values for all
  the GO terms being used that will be used for calculations instead of R
  packages. The vector saved in the file must be called "IA" or "IAccr" (see computeIA). This argument is optional.}
  \item{outfile}{A file in which to save the function output as an R object. Defaults to
                  "rumicurve.rda". Should be a character vector.}
>>>>>>> fce8d5b26e6ecfb454a764e7aa3c606d22b60638
  \item{add.weighted}{A boolean indicating whether or not to add calculation of information
                    content weighted versions of RU, MI, and SS to the output.}
  \item{add.prec.rec}{A boolean indicating whether or not to calculate precision, recall
                      and specificity values for the prediction at each threshold and add
                      to the output.}
}

\value{
Returns a named list with the same number of elements as the input "predfiles".
Each element is a data frame containing all of the user-requested values for the data at
each threshold.
}

\author{
Ian Gonzalez and Wyatt Clark
}

\seealso{
\code{\link{computeIA}}
\code{\link{findRUMI}}
}

\examples{
# Using test data sets from SemDist, plot a RUMI curve:
truefile <- system.file("extdata", "MFO_LABELS_TEST.txt", package="SemDist")
predfile <- system.file("extdata", "MFO_PREDS_TEST.txt", package="SemDist")
avgRUMIvals <- RUMIcurve("MF", "human", 0.05, truefile, predfile)
firstset <- avgRUMIvals[[1]]
plot(firstset$RU, firstset$MI)
}