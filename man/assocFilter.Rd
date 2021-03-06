% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_rbedInfo.R
\name{assocFilter}
\alias{assocFilter}
\title{Filter a PlGwasC object by the results of a \code{plink --assoc} run}
\usage{
assocFilter(pl_gwas, plink_out_stem = NULL, p_threshold = 0.1,
  db_setup = FALSE, force = TRUE)
}
\arguments{
\item{pl_gwas}{PlGwasC object}

\item{plink_out_stem}{character. Output plink file stem (without .bed extension). The default is to add a "_filtered_{RANDOM_ID}" suffix to the original.}

\item{p_threshold}{numeric. P-value threshold.}

\item{db_setup}{logical. Whether to setup the PlGwasC object.}

\item{force}{logical. Overwrite existing PLINK files.}
}
\value{
a new PlGwasC object.
}
\description{
This is meant for reduction in computational burden. The \code{plink --assoc} does not
accept covariates makes some assumptions accordingly, and thus runs faster than \code{--linear} and
\code{--logistic}. SNPs that does not produce a p-value more significant than a user-set threshold will
be filtered out. A new PLINK file is made and a corresponding new PlGwasC object is returned.
}
\examples{
\dontrun{
rbed_info = rbedInfo(bedstem = "mmp13", db_setup = FALSE)
pl_gwas = plGwas(rbed_info, 
		pheno = "mmp13.phe",
		pheno_name = "Page", 
		gwas_tag = "mmp13_page_sex_age")
runGwas(pl_gwas)
x = readGwasOut(pl_gwas, c("SNP", "P"), rmGwasOut = FALSE)
pl_gwas1 = assocFilter(pl_gwas, p_threshold = 0.001)
runGwas(pl_gwas1)
x1 = readGwasOut(pl_gwas1, c("SNP", "P"), rmGwasOut = FALSE)
y = dplyr::inner_join(x, x1, by = "SNP")
all(y$P.x == y$P.y)
all(y$P.y < 0.001)
}

}
\author{
Kaiyin Zhong, Fan Liu
}

