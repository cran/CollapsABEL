# TODO: Add comment
# 
# Author: kaiyin
###############################################################################

#' Find cytoband at a given position
#' @param chr integer or character. Chromosome number. If it's an integer it should be in range [1, 22]. If it's a string it's should be in the format as "chr1, chr2, ..., chr22, chrX, chrY"
#' @param pos integer. Position on chromosome.
#' @param ref character. Reference data. Should be either "hg18" or "hg19"
#' @return Vector of cytobands.
#' 
#' @author kaiyin
#' @export
cytoband = function(chr, pos, ref = "hg19") {
	if (is.numeric(chr)) {
		stopifnot(all(1 <= chr & chr <= 22))
		chr = paste0("chr", chr)
	} else {
		correctFormat = stringr::str_match(pattern="^chr(([0-9]{1,2})|([XY]))$", string = chr)
		incorrectFormat = which(is.na(correctFormat[, 3]) & is.na(correctFormat[, 4]))
		if(length(incorrectFormat) != 0) {
			message(sprintf("Chr not in the right format, index: %d", incorrectFormat))
			print(chr[incorrectFormat])
			stop()
		}
		chrInt = as.integer(na.omit(correctFormat[, 3]))
		stopifnot(
				all(1 <= chrInt && chrInt <= 22)
		)
	}
	stopifnot(ref == "hg18" || ref == "hg19")
	if(length(chr) == 1) {
		chr = rep(chr, length(pos))
	} else {
		stopifnot(length(pos) == length(chr))
	}
	ref = get(ref)
	idx = sapply(1:length(pos), function(j) {
				p = pos[j]
				chr = chr[j]
				i = which(ref$chr == chr & ref$start < p & p <= ref$end)	
				if(length(i) == 0) {
					.Machine$double.xmax # use an extremely large values so that the indexing will return NA
				} else {
					i
				}
			})
	ref$cyto[idx]
}

