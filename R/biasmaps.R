

loadBiasmapInfo <- function(biasmap) {
  read.delim(file=paste(biasmap.dir, '/info.txt', sep=''), row.names=1)
}

buildBiasmapFromBSgenome <- function(
  bsgenome,
  pattern,
  chromosomes,
  scales,
	stepsize=100, D=4,
  n.cores=c(chromosomes=1,scales=1),
  biasmaps.dir
  ) {
	

	# Find 'pattern' sites
	cat('Matching pattern ', pattern, ' to ', paste(chromosomes, collapse=', '), '...', '\n', sep='')
	bgsites <- mclapply(chromosomes, mc.cores=n.cores[1], FUN=function(chr) {
		start(matchPattern(pattern, bsgenome[[chr]], max.mismatch=0))
	})
	names(bgsites) <- chromosomes

  
  # Create the output dir, f.e. bsgenome-mm9-TA-S100-D4
  output.dir <- paste(biasmaps.dir, '/bsgenomes-', providerVersion(bsgenome), '-', pattern, '-stepsize', stepsize, '-D', D, sep='')
  dir.create(output.dir, showWarnings=FALSE, recursive=TRUE)
	
  
  # Writing info file
	info.file <- paste(output.dir, '/info.txt', sep='')

	info.df <- data.frame(
		chromosomes,
		seqlengths(bsgenome)[chromosomes],
		sapply(bgsites, length)
	)
	colnames(info.df) <- c('chromosome', 'length', 'TA count')

	write.table(info.df, file=info.file, sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

  
  # for each scale, generate the background counts and densities
  
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
			cat('  Calculating background density and counts at scale=', scale, '...', '\n', sep='')

			densObj <- density(bgsites[[chr]], bw=scale, from=ftn$from, to=ftn$to, n=ftn$n, kernel='rectangular')
			x       <- densObj$x
			
			bgdens  <- densObj$y * sl / sum(as.numeric(seqlengths(bsgenome)[chromosomes]))
			bgcount <- block.convolve(x, bgsites[[chr]], width=scale * D * 2)

			grdens  <- GRanges(seqnames=chr, ranges=IRanges(x, x), score=bgdens, seqlengths=seqlengths(bsgenome))
			grcount <- GRanges(seqnames=chr, ranges=IRanges(x, x), score=bgcount, seqlengths=seqlengths(bsgenome))

			list(grdens=grdens, grcount=grcount)
		})

		# write to wigs
		wig.count.file <- paste(output.dir, '/count-scale', as.integer(round(scale)), '.wig', sep='')
		wig.density.file <- paste(output.dir, '/density-scale', as.integer(round(scale)), '.wig', sep='')

		export.wig(con=wig.density.file, do.call('c', lapply(grList, '[[', 'grdens')))
		export.wig(con=wig.count.file, do.call('c', lapply(grList, '[[', 'grcount')))
	})
  
  # return the name
	paste(providerVersion(bsgenome), '_', pattern)
}




