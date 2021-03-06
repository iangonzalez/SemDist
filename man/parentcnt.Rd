\name{parentcnt}
\alias{parentcnt}
\docType{data}
\title{
Parent Count Data  
}
\description{
This data set contains the parent count values for each term in the requested ontology/species (the number of times that the term's parents annotate a protein). This can be used along with the term count to calculate information accretion.
}
\usage{parentcnt}
\format{
  A named numeric vector with one value corresponding to each GO accession in the ontology.
}
\examples{
data("Parent_Count_mouse_CC", package = "SemDist")
str(parentcnt)
}
\source{
The gene ontology data was obtained from the GO.db package and the annotation data was obtained from the following packages for each organism: 

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
\keyword{datasets}