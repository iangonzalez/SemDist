\name{findRUMI}
\alias{findRUMI}
\title{findRUMI}

\description{
  Reads in a file containing the true terms annotating a set of sequences and a file
  containing the predicted terms and scores for a set of sequences and outputs a data
  frame containing the remaining uncertainty and misinformation values for the predictions
  made for each sequence.
}

\usage{
findRUMI(ont, organism, threshold = 0.05, truefile="trues.txt",
          predfile="preds.txt", outfile="rumi.rda")
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
  \item{outfile}{Filename to save the results of the function call as an R object.}
}

\value{
  A data frame containing the RU and MI values for each sequence in the file.
}
