% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0_readFactory.R
\name{readPlinkOut}
\alias{readPlinkOut}
\title{Read plink output files}
\usage{
readPlinkOut(filename, ...)
}
\arguments{
\item{filename}{Filenames of plink output files, see \code{collenv$.plink_out_ext}}

\item{...}{passed to one of \code{readAssoc, readQassoc, readLinear, readLogistic}}
}
\value{
data.frame
}
\description{
Read plink output files
}
\examples{
\dontrun{
dat1 = readPlinkOut("assoc/mmp13.assoc")
dat2 =    readAssoc("assoc/mmp13.assoc")
all(na.omit(dat1 == dat2))
dat1 = readPlinkOut("assoc/mmp13.assoc", c("CHR", "SNP", "P", "OR"))
dat2 =    readAssoc("assoc/mmp13.assoc", c("CHR", "SNP", "P", "OR"))
all(na.omit(dat1 == dat2))
dat1 = readPlinkOut("assoc/mmp13.qassoc")
dat2 =   readQassoc("assoc/mmp13.qassoc")
all(na.omit(dat1 == dat2))
dat1 = readPlinkOut("assoc/mmp13.qassoc", c("CHR", "SNP", "P", "R2"))
dat2 =   readQassoc("assoc/mmp13.qassoc", c("CHR", "SNP", "P", "R2"))
all(na.omit(dat1 == dat2))
}

}
\author{
Kaiyin Zhong, Fan Liu
}

