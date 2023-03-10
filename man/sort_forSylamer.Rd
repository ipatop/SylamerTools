% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sylamer_Tools.R
\name{sort_forSylamer}
\alias{sort_forSylamer}
\title{Sort differential gene expression resutls}
\usage{
sort_forSylamer(
  DEres = "res.txt",
  sortBy = "log2FoldChange",
  outName = "ForSylamer.txt",
  return = F
)
}
\arguments{
\item{DEres}{table with differential expression, one column must be gene names and must be named "gene", there has to be at least another column that will be used to sort genes}

\item{sortBy}{name of column used to sort, default is "log2FoldChange"}

\item{outName}{name of the file to write the sorted list to. Deafult is ="ForSylamer.txt"}

\item{return}{Boolean, if True, then return the sorted list. Default is False.}
}
\value{
return the sorted list of genes list to run in sylamer
}
\description{
This function sort results from differential gene expression by any parameter.
}
\examples{

#sort_forSylamer(DEres="res.txt",sortBy="log2FoldChange",outName="ForSylamer.txt")



}
