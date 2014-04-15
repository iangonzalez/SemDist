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
RUMIcurve(ont, organism, increment = 0.05, predfiles, truefile, 
          outfile = "rumicurve.rda", add.weighted = FALSE, 
          add.prec.rec = FALSE)
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
  \item{predfiles}{A character vector containing which files to read in as the predicted
                    annotations. Should be tab-delimited, with the first column containing
                    sequences, the second column containing GO accessions, and the third
                    column containing the predictors 0-1 score for that prediction.}
  \item{truefile}{A character vector indicating the file from which to read the true
                  annotations for the given sequences. Should be tab-delimited, with the
                  first column containing the sequence ids and the second containing GO
                  accessions.}
  \item{outfile}{A file in which to save the function output as an R object. Defaults to
                  "rumicurve.rda". Should be a character vector.}
  \item{add.weighted}{A boolean indicating whether or not to add calculation of information
                    content weighted versions of RU, MI, and SS to the output.}
  \item{add.prec.rec}{A boolean indicating whether or not to calculate precision, recall
                      and specificity values for the prediction at each threshold and add
                      to the output.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a named list with the same number of elements as the input "predfiles".
Each element is a data frame containing all of the user-requested values for the data at
each threshold.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ian Gonzalez and Wyatt Clark
}
\note{
%%  ~~further notes~~
}

\seealso{
\code{\link{computeIA}}
\code{\link{findRUMI}}
\code{\link{SemDist}}
}