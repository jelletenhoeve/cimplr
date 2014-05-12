

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

buildKSEDistributions <- function(scales, D=4, nBins=1000, kNormal=120, n.cores=1, output.dir) {

	kse.dists <- mclapply(scales, mc.cores=n.cores, FUN=function(scale) {
		cat('calculating KSEDistribution scale = ', scale, '...\n', sep='')
		new('KSEDistribution', scale=scale, nBins=nBins, kNormal=kNormal)
	})

	names(kse.dists) <- .formatScales(scales)



	kse.dists.object <- list(
		cumdens = lapply(x@histograms, '[[', 'cumdens'); names(val) <- 1:length(x); val),
		cumdens = lapply(x@histograms, '[[', 'cumdens'); names(val) <- 1:length(x); val)
	)

	#file <- paste(output.dir, '/KSEDistributions.rda', sep='')
	#cat('writing ', file, '...', sep='')
	#save(kse.dists, file=file)
	#cat('done!\n')


	kse.dists.object
}

buildReferences <- function(pattern, bsgenome, chromosomes, scales,
	stepsize=1000, D=4, n.cores=c(1,1),
	output.dir) {
	

	# Load 'pattern' sites
	cat('Matching pattern ', pattern, ' to ', paste(chromosomes, collapse=', '), '...', '\n', sep='')
	bgsites <- mclapply(chromosomes, mc.cores=n.cores[1], FUN=function(chr) {
		start(matchPattern(pattern, bsgenome[[chr]], max.mismatch=0))
	})
	names(bgsites) <- chromosomes

	info.file <- paste(output.dir, '/', providerVersion(bsgenome), '.', pattern, '.txt', sep='')

	info.df <- data.frame(
		chromosomes,
		seqlengths(bsgenome)[chromosomes],
		sapply(bgsites, length)
	)
	colnames(info.df) <- c('chromosome', 'length', 'TA count')

	write.table(info.df, file=info.file, sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

	mclapply(scales, mc.cores=n.cores[2], FUN=function(scale) {
		grList <- mclapply(chromosomes, mc.cores=n.cores[1], FUN=function(chr) {
			cat(' ', chr, '\n', sep='')
			sl <- seqlengths(bsgenome)[chr]

			ftn <- list(
				from = 1,
				to   = sl - (sl %% stepsize) + 1
			)
			ftn$n <- 1 + ((ftn$to - ftn$from) / stepsize)

			# background density
			cat('  Calculating background density: bgd(x) at scale=', scale, '...', '\n', sep='')

			densObj <- density(bgsites[[chr]], bw=scale, from=ftn$from, to=ftn$to, n=ftn$n, kernel='rectangular')
			x       <- densObj$x
			
			bgdens  <- densObj$y * sl / sum(as.numeric(seqlengths(bsgenome)[chromosomes]))
			bgcount <- block.convolve(x, bgsites[[chr]], width=scale * D * 2)

			grdens  <- GRanges(seqnames=chr, ranges=IRanges(x, x), score=bgdens, seqlengths=seqlengths(bsgenome))
			grcount <- GRanges(seqnames=chr, ranges=IRanges(x, x), score=bgcount, seqlengths=seqlengths(bsgenome))

			list(grdens=grdens, grcount=grcount)
		})

		# write to wig
		wig.count.file <- paste(output.dir, '/', providerVersion(bsgenome), '.count',
			pattern, '.scale', as.integer(round(scale)), '.D', D, '.stepsize', stepsize, '.wig', sep=''
		)
		wig.density.file <- paste(output.dir, '/', providerVersion(bsgenome), '.density',
			pattern, '.scale', as.integer(round(scale)), '.stepsize', stepsize, '.wig', sep=''
		)


		export.wig(con=wig.density.file, do.call('c', lapply(grList, '[[', 'grdens')))
		export.wig(con=wig.count.file, do.call('c', lapply(grList, '[[', 'grcount')))
	})
}




