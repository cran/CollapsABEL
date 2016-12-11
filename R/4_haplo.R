# Translates a genotype vec into a format that haplo.stats understands:
# 0 -> 1, 1
# 1 -> 1, 2 
# 2 -> 2, 2
haplo.g = function(g) {
	a1 = (g > 1) + 1
	a2 = (g > 0) + 1
	data.frame(a1, a2)
}



running.bool = function(true_id, len) {
	stopifnot(all(true_id >= 1 & true_id <= len))
	res = rep(FALSE, len)
	res[true_id] = TRUE 
	res
}

# given a vector of maxid and a vector of lengths,
# return a bool vector
# for example:
# maxid = c(1 2 1 3)
# l     = c(2 3 4 5), then
# ret = c(T, F;          # 2 bools, 1st chosen
#         F, T, F;       # 3 bools, 2nd chosen
#         T, F, F, F;    # 4 bools, 1st chosen
#         F, F, T, F, F) # 5 bools, 3rd chosen
running.bools = function(maxid, len) {
	stopifnot(length(maxid) == length(len))
	res = lapply(1:length(len), function(i) running.bool(maxid[i], len[i]))
	do.call(c, res)
}

haploFormat = function(geno) {
	haplo = data.frame(matrix(NA, nrow(geno), 0))
	geno_ncol = ncol(geno)
	for(i in 1:geno_ncol) {
		# each SNP gets two alleles
		haplo = cbind(haplo, haplo.g(geno[, i]))
	}
	old_cnames = rep(colnames(geno), each = 2)
	cnames = paste0(old_cnames, ".a", rep(1:2, geno_ncol))
	colnames(haplo) = cnames
	haplo
}

#' Infer haplotypes for a pair of SNPs
#' 
#' @param geno Genotype data frame. Must have 4 columns, the first two being "FID" and "IID", the last two being the genotypes.
#' @param format_idx Column indices used for formatting haplotype string.
#' @return A data frame of haplotypes
#' @importFrom haplo.stats haplo.em
#' @importFrom dplyr group_by summarise id 
#' 
#' @author kaiyin
getHaplo = function(geno, format_idx = NULL) {
	stopifnot(all(c("FID", "IID") %in% colnames(geno)))
	geno = geno[complete.cases(geno), ]
	nsnp = ncol(geno)
	labels = colnames(geno)[3:nsnp]
	geno$id = 1:nrow(geno)
	save.em = haplo.em(haploFormat(geno[, 3:nsnp]), labels)
	hapdat = data.frame(id = save.em$subj.id, h1 = save.em$hap1code, h2 = save.em$hap2code, p = save.em$post)
	id_groups = group_by(hapdat, id)
	idx =  summarise(id_groups, maxid = which.max(p))
	len   =  summarise(id_groups, len = length(p))
	select = running.bools(idx$maxid, len$len)
	hapdat = hapdat[select, ]
	haptab = save.em$haplotype
	haptab$h = 1:nrow(haptab)
	hapdat = dplyr::left_join(hapdat, haptab, by = c("h1" = "h"))
#	geno.hapdat = renameLoci(geno.hapdat, 1)
	hapdat = dplyr::left_join(hapdat, haptab, by = c("h2" = "h"))
	h = haploString(hapdat, snp_idx = 5:ncol(hapdat), format_idx)
	dplyr::left_join(geno[, c("FID", "IID", "id")], h[, c("id", "diplo")], by = "id")
}

haploString = function(hapdat, snp_idx = NULL, format_idx = NULL, generic_diplo_name = TRUE) {
	if(is.null(snp_idx)) {
		snp_idx = 1:nc
	}
	hapdat1 = hapdat[, snp_idx]
	nc = ncol(hapdat1)
	stopifnot(nc %% 2 == 0)
	half = nc / 2
	stopifnot(length(format_idx) <= half)
	stopifnot(all(format_idx > 0))

	format1 = function(x, idx1, idx2)  {
		s1 = paste(x[idx1], collapse = "-")
		s2 = paste(x[idx2], collapse = "-")
		s1 = paste0("-", s1, "-")
		s2 = paste0("-", s2, "-")
		paste0(s1, "+", s2)
	}
	
	cnames = gsub(x = colnames(hapdat1), pattern = "\\.[xy]$", replacement = "")
	newcolname = function(idx) {
		paste(cnames[idx], collapse = "_")
	}
	
	if(is.null(format_idx)) {
		idx1 = 1:half
		idx2 = idx1 + half
		f = function(x) {
			format1(x, idx1, idx2)
		}
		if(generic_diplo_name) {
			hapdat[, "diplo"] = apply(hapdat1, 1, f)
		} else {
			hapdat[, newcolname(idx2)] = apply(hapdat1, 1, f)
		}
		hapdat
	} else {
		stopifnot(max(format_idx) <= half)
		idx2 = format_idx + half
		f = function(x) {
			format1(x, format_idx, idx2)
		}
		if(generic_diplo_name) {
			hapdat[, "diplo"] = apply(hapdat1, 1, f)
		} else {
			hapdat[, newcolname(idx2)] = apply(hapdat1, 1, f)
		}
		hapdat
	}
}

#' Inferring haplotypes from two genotype data frames, and join with phenotypes
#' @param g1 First genotype data frame
#' @param g2 Second genotype data frame, must be of the same dimension as the first. The first two column must be FID and IID.
#' @param phe Phenotype data frame, the first two columns must be FID and IID
#' @param pool A genotype data frame, assumed to be different from g1 and g2, used for pooling.
#' @return A data frame containing phenotype and haplotype for each individual.
#' 
#' @author kaiyin
#' @export
getHaplos = function(g1, g2, phe, pool = NULL) {
	stopifnot(all(dim(g1) == dim(g2)))
	stopifnot(all(g1$FID == g2$FID))
	stopifnot(all(g1$IID == g2$IID))
	nc = ncol(g1)
	stopifnot(nc > 2)
	stopifnot(ncol(phe) > 2)
	stopifnot(all(c("FID", "IID") %in% colnames(phe)))
	stopifnot(all(c("FID", "IID") %in% colnames(g1)))
	stopifnot(all(c("FID", "IID") %in% colnames(g2)))
	ret = list()
	for (i in 3:nc) {
		x = dplyr::left_join(g1[, c(1, 2, i)], g2[, c(1, 2, i)], by = c("FID", "IID"))
		if(!is.null(pool)) {
			ncpool = ncol(pool)
			stopifnot(ncpool > 4) # fid and iid should be there
			stopifnot(all(c("FID", "IID") %in% colnames(geno_pool)))
			idx = c(1, 2, sample(3:ncpool, 3))
			x = dplyr::left_join(x, pool[, idx], by = c("FID", "IID"))
		}
		h = getHaplo(x, 1:2)
		dat = dplyr::inner_join(phe, h[, c("FID", "IID", "diplo")], by = c("FID", "IID"))
		name = paste0(names(g1)[i], "_", names(g2)[i])
		ret[[name]] = dat
	}
	ret
}
