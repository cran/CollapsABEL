#' Retrive SNP positions from UCSU database
#' @param snps A vector of SNP names
#' @param rm_underscore Remove irregular chromosome names
#' @param ref Either "hg18" or "hg19"
#' @param snpdb Either "snp138" or "snp137"
#' @return A data frame containing positions of given SNPs
#' 
#' @author kaiyin
#' @export
snpPos = function(snps, rm_underscore = TRUE, ref = c("hg18", "hg19"), snpdb = c("snp138", "snp137")) {
	cmdStart = sprintf("mysql --user=genome --host=genome-mysql.cse.ucsc.edu -NA     -e \"select chrom, chromStart, chromEnd, name from %s.%s where name in ", ref, snpdb)
	tmpfile = tempfile()
	cmdEnd = sprintf("\" %s > %s", ref, tmpfile)
	snplist = strVectorRepr(snps, start_with_c = FALSE)
	cmd = sprintf("%s %s %s", cmdStart, snplist, cmdEnd)
	system(cmd)
	res = read.table(tmpfile, header = FALSE)
	colnames(res) = c("chrom", "chromStart", "chromEnd", "SNP")
	hasUnderscore = stringr::str_detect(res$chrom, "_")
	res[!hasUnderscore, ]
}