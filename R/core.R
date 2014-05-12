

fromton <- function(x, scale, sample.period) {

	step <- floor(scale * sample.period) 
	
	# default from = 1
#	from <- 1
#	to <- object@seqlength
	from <- min(x) - (step/2)
	to   <- max(x) + (step/2)
	
	# n is the number of sample points
	n <- 2 + ( (to - from - 1) %/% step )
	
	# update 'to'
	to <- ((n-1) * step) + from
	
	list(from=from, to=to, n=n)
}


log.seq <- function(from, to, length.out) {
	exp(seq(log(from), log(to), length.out=length.out))
}

segments <- function(x) {
	# 'x' is a logical vector
	
	x <- c(FALSE, x, FALSE) # to make sure we capture the boundary segments as well
	
	trans_vec <- x[-length(x)] - x[-1]
	# at 0 --> 1 a region start
	# at 1 --> 0 a the region ends
	start_pos <- which(trans_vec == -1)
	end_pos <- which(trans_vec == 1) - 1
	
	list(start.pos=start_pos, end.pos=end_pos, n=length(start_pos))
}


whichLocalMaxima <- function(x, min.val=0.1) {
	####### a port of Jeroens matlab code #######
	diff_vec            <- x[-1] - x[-length(x)]
	# Make a tri-state vector indicating the transitions
	tristate <- sign(diff_vec)
	# Remove the flat parts
	non_flat_index      <- which(tristate != 0)
	bistate             <- tristate[tristate != 0]
	# Find the 1-->-1 transitions
	trans_vec           <- bistate[-length(bistate)] - bistate[-1]
	trans_index         <- which(trans_vec == 2);
	if (length(trans_index) == 0) {
		list(max_pos=c(), max_val=c())
	}
	max_pos             <- non_flat_index[trans_index+1];
	max_val             <- x[max_pos];
	#############################################
	
	max_pos[max_val > min.val]
}





##
## Associate a CIS with genes
##

.associateGenes <- function(cisStart, cisEnd, peakLoc, chr_genes, margin = 100e3) {
	
	# the genes within the region + margin
	gene.idx <- which(chr_genes$start_position < (cisEnd + margin) & chr_genes$end_position > (cisStart - margin))
	
	if (length(gene.idx) == 0) {
		list(
			associated = numeric(0),
			other      = numeric(0)
		)
	} else if (length(gene.idx) == 1) {
		list(
			associated = gene.idx,
			other      = numeric(0)
		)
	} else {

		# calculate the distances
		dists <- pmin(abs(chr_genes$start_position[gene.idx] - peakLoc), abs(chr_genes$end_position[gene.idx] - peakLoc))
		
		# genes which contain the peak get distance 0
		dists[chr_genes$start_position[gene.idx] <= peakLoc & chr_genes$end_position[gene.idx] >= peakLoc] <- 0


		gene.order <- order(dists, decreasing=FALSE)	
		
		min.ties <- length(which(dists == min(dists)))
		
		list(
			associated = gene.idx[gene.order][1:min.ties],
			other      = gene.idx[gene.order][-(1:min.ties)]
		)
	}
}

.formatScales <- function(scales) {

	as.character(as.integer(scales))

	#str <- as.character(scales)
	#k.idx <- scales %% 1e3 == 0
	#m.idx <- scales %% 1e6 == 0
	#str[k.idx] <- paste(scales[k.idx] / 1e3, 'k', sep='')
	#str[m.idx] <- paste(scales[k.idx] / 1e6, 'm', sep='')
	#str
}


getEnsemblGenes <- function(chromosomes, geneIdentifiers=c('ensembl_gene_id', 'external_gene_id'), mart=useMart("ensembl", dataset = "mmusculus_gene_ensembl")) {
  genes <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position'), filters='chromosome_name', values=substring(chromosomes, 4), mart=mart)
  
  # , 'unigene'
  # filter out unknown 
  if (!all(i <- geneIdentifiers %in% listAttributes(mart)[,1])) {
    stop(paste(geneIdentifiers[!i], ' is/are no valid biomaRt attributes for this mart.'))
  }
  
  if (length(geneIdentifiers) > 0) {
    if ( !(length(geneIdentifiers) == 1 & geneIdentifiers[1] == 'ensembl_gene_id' )) {
      annot <- getBM(attributes=unique(c('ensembl_gene_id', geneIdentifiers)), filters='ensembl_gene_id', values=genes$ensembl_gene_id, mart=mart)
      
      extra_annot <- sapply(genes$ensembl_gene_id, function(ens) {
        idx <- annot$ensembl_gene_id == ens
        sapply(2:ncol(annot), function(col) {
          vals <- annot[idx, col]
          vals[vals==''] <- NA
          paste(unique(na.omit(vals)), collapse='|')
        })
      })
      
      if (!is.null(dim(extra_annot)))
        extra_annot <- t(extra_annot)
      
      extra_annot <- data.frame(extra_annot, stringsAsFactors=FALSE)
      colnames(extra_annot) <- colnames(annot)[-1]
      genes <- cbind(genes, extra_annot)
    }
  }
  if (! 'ensembl_gene_id' %in% geneIdentifiers) {
    genes$ensembl_gene_id <- NULL
  }
  
  genes$chromosome_name <- paste('chr', genes$chromosome_name, sep='')
  
  attr(genes, 'geneIdentifiers') <- geneIdentifiers
  
  genes
}


block.convolve <- function(x, bg, width, verbose=FALSE) {
	if (verbose) {
		cat('x        =', x, '\n')
		cat('bg       =', bg, '\n')
		cat('width    =', width, '\n')
		cat('\n')
	}

	vec <- sort(c(x, bg+.1))
	x.indices  <- match(x, vec)
	vec.is.x <- vec %in% x
	bg.indices <- match(bg+.1, vec)
	vec <- sort(c(x, bg))

	n <- ceiling(width / (x[2] - x[1]) / 2)

	val <- rep(0, length(x))

	if (verbose) {
		cat('vec        =', vec, '\n')			
		cat('x.indices  =', x.indices, '\n')			
		cat('bg.indices =', bg.indices, '\n')			
		cat('vec.is.x   =', vec.is.x, '\n')
		cat('n          =', n, '\n')
		cat('\n')
	}


	for (i in 1:length(x)) {
#		if (i %% 100 == 0) {
#			cat(i/length(x), '...\r')
#		}

		vec.idx <- x.indices[i]
		loc <- vec[vec.idx]

		if ((i - n) < 1) {
			start.vec.idx <- 1
		} else {
			start.vec.idx <- x.indices[i - n]
		}

		if ((i + n) > length(x)) {
			end.vec.idx   <- length(vec)
		} else {
			end.vec.idx   <- x.indices[i + n]
		}

		if (verbose) {
			cat('i=', i, '\n')			
			cat(' vec.idx       =', vec.idx, '\n')			
			cat(' loc           =', loc, '\n')			
			cat(' start.vec.idx =', start.vec.idx, '\n')			
			cat(' end.vec.idx   =', end.vec.idx, '\n')			
		}


		val[i] <- sum(
			abs(vec[start.vec.idx:end.vec.idx] - loc) <= width/2 & 
			!vec.is.x[start.vec.idx:end.vec.idx]
		)
	}



	val
}




