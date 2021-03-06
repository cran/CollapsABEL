% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_bedinfo.R
\name{fidIid}
\alias{fidIid}
\title{FID and IID columns from fam file}
\usage{
fidIid(pl_info)
}
\arguments{
\item{pl_info}{PlInfoC object}
}
\value{
data.frame of two columns "FID" and "IID"
}
\description{
FID and IID columns from fam file
}
\examples{
\dontrun{
pl_info = plInfo(bedstem = "mmp13", db_setup = TRUE)
fidiid = fidIid(pl_info)
fam = readFam("mmp13.fam", c("FID", "IID"))
all(fam == fidiid)
}

}
\author{
Kaiyin Zhong, Fan Liu
}

