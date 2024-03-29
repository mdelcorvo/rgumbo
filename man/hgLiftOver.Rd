% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hgLiftOver.R
\name{hgLiftOver}
\alias{\%liftover\%}
\alias{hgLiftOver}
\title{hgLiftOver}
\usage{
hgLiftOver(df,from="hg38",to="hg19")
}
\arguments{
\item{df}{data.frame in BED format with chromosome number, start and end position or chromosome number and a single genomic position respectively}
\item{from}{starting assembly}
\item{to}{final assembly}
}
\value{
}
\description{
Function to converts genome coordinates between assemblies
}
\details{
}
\examples{
hgLiftOver(df = bed_file,from="hg38",to="hg19")
}
