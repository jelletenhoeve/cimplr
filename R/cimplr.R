# cimplr.R
#
# Description: CIMPLR main function, from insertions to CIS.
# Author: Jelle ten Hoeve
#
#
# Arguments:
#   insertions:                   a data.frame with 'chr' and 'location' columns.
#   reference:                    the reference genome to which the insertion are mapped --> buildReference
#   pattern:                      the background pattern, f.e. TA for Sleaping Beauty data
#   exclude.chromosomes:          these chromosomes are excluded, typically the donor chromosome(s).
#   scales:                       the set of scales for kernel convolutions, default: log.seq(5000, 500000, 100)
#   alpha.level:                  the alpha level at which CISes are called are multiple testing, default: 0.05
#   p.adjust.method:              the multiple testing correction method. 'n.bgsites' for 'n.peaks'
#   p.adjust.n.method:            the method to calculate the number of tests used for multiple testing
#   test.chromosomes.seperately:  test each chromosome seperately, or test the genome and a whole
#   genes:                        a data.frame of genes used for annotating the CIS. The data.frame should contain columns: 'start_position', 'end_position' and 'chromosome_name'. An attr(genes, 'geneIdentifiers') should point to the "geneIdentifiers" column(s).
#   output.dir:                   the directory to which the output is written. This directory should exist.
#   sample.period:                this is sample period at which KSE(x) is sampled. This is dependent on scale, f.e. at scale 30000, x is sampled every .1 * 30000 = 3000 bp.
#   n.cores:                      a vector of length 2. First element is used for loops over scales, second element for loops over chromosomes
#

cimplr <- function(insertions,
	reference='mm9', pattern='TA',

	exclude.chromosomes=c(),

	scales=log.seq(5000, 500000, 100),
	pkse.method=c('count', 'density'),
	alpha.level=.05,
	p.adjust.method=c("bonferroni", "fdr", "none"),
	p.adjust.n.method=c('n.x', 'n.bgsites', 'n.peaks', 'n.bgsites/scale'),
	mtf=1,
	test.chromosomes.seperately=FALSE,
	stepsize=1000,
	D=4,

	genes.file='data/genes.rda',
	reference.dir='./references',
	n.cores=c(1,1)

	) {





	#
	# argument checks and processing
	#
	if (missing(insertions)) {
		stop("Need to provide 'insertions'.")
	}

	if (!all(c("chr", "location") %in% colnames(insertions))) {
		stop("'insertions' must contain 'chr' and 'location' columns")
	}


	if (alpha.level <= 0 | alpha.level >= 1) {
		stop("Need to provide 0 < 'alpha.level' < 1.")
	}


	pkse.method <- match.arg(pkse.method)
	p.adjust.method <- match.arg(p.adjust.method)
	p.adjust.n.method <- match.arg(p.adjust.n.method)




	# determine the available chromosomes
	chromosomes <- unique(as.character(insertions$chr))
	chromosomes <- chromosomes[!chromosomes %in% exclude.chromosomes]
	#chromosomes <- c('chr2', 'chr3', 'chr15', 'chr18')
	#chromosomes <- c('chr2')
	#chromosomes <- c('chr2', 'chr3', 'chr4')



	# make 'insertions' a list of insertions
	insertions.org <- insertions
	insertions  <- insertions[insertions$chr %in% chromosomes, ]
	insertions  <- split(insertions$location, insertions$chr)[chromosomes]
	insertions  <- insertions[chromosomes]






	# load the reference info and make a chr.info
	chr.info <- read.delim(file=paste(reference.dir, '/', reference, '.', pattern, '.txt', sep=''), row.names=1)
	chr.info <- chr.info[chromosomes, ]
	chr.info$n.insertions=sapply(insertions, length)

	#chr.info.totals <- sapply(chr.info, sum)
	n.insertions.total <- sum(chr.info[[3]])
	n.bgsites.total <- sum(chr.info[[2]])
	n.bp.total <- sum(chr.info[[1]])

	input <- as.list(environment())




	#
	# cat a summary
	#
	cat('*** CIMPLR *** \n\n')
	str(input)

	cat('\n')

	list(
		input = input
	)

}



cimplr.replace_insertions <- function(cimplr.object, new.insertions) {

	message('cimplr.replace_insertions()')
	with(cimplr.object$input, {

		insertions <- new.insertions
		insertions.org <- insertions
		insertions  <- insertions[insertions$chr %in% chromosomes, ]
		insertions  <- split(insertions$location, insertions$chr)[chromosomes]
		insertions  <- insertions[chromosomes]






		# load the reference info and make a chr.info
		chr.info <- read.delim(file=paste(reference.dir, '/', reference, '.', pattern, '.txt', sep=''), row.names=1)
		chr.info <- chr.info[chromosomes, ]
		chr.info$n.insertions=sapply(insertions, length)

		#chr.info.totals <- sapply(chr.info, sum)
		n.insertions.total <- sum(chr.info[[3]])
		n.bgsites.total <- sum(chr.info[[2]])
		n.bp.total <- sum(chr.info[[1]])

		input <- as.list(environment())


		##

		## clear the fields in data

		scale.objects <- mclapply(cimplr.object$data, mc.cores=n.cores[1], FUN=function(scale.object) {
			mclapply(scale.object, mc.cores=n.cores[2], FUN=function(chr.object) {
				if ('th' %in% names(chr.object)) {
					chr.object[c("chr", "scale", "x", "bg", "n.insertions", "n.insertions.total", "n.bgsites", "n.bgsites.total", "n.bp", "n.bp.total", "th")]
				} else {
					chr.object[c("chr", "scale", "x", "bg", "n.insertions", "n.insertions.total", "n.bgsites", "n.bgsites.total", "n.bp", "n.bp.total")]
				}
			})
		}) # end scale space


		list(
			input = input,
			data  = scale.objects,
			kse.distributions = cimplr.object$kse.distributions
		)

	})
}


cimplr.reference <- function(cimplr.object) {
	message('cimplr.reference()')

	with(cimplr.object$input, {	


		#
		# 1. KSE convolution and each scale
		#
		scale.objects <- mclapply(scales, mc.cores=n.cores[1], FUN=function(scale) {
			
			message('cimplr.reference() - load reference scale = ', scale)

			if (pkse.method == 'count') {
				wig.count.file <- paste(reference.dir, '/', reference, '.count',
					pattern, '.scale', as.integer(round(scale)), '.D', D, '.stepsize', stepsize, '.wig', sep=''
				)

				bgranges <- import.wig(con=wig.count.file, asRangedData=FALSE)
			}

			if (pkse.method == 'density') {
				wig.density.file <- paste(reference.dir, '/', reference, '.density',
					pattern, '.scale', as.integer(round(scale)), '.stepsize', stepsize, '.wig', sep=''
				)
				bgranges <- import.wig(con=wig.density.file)
			}

			chr.objects <- mclapply(chromosomes, mc.cores=n.cores[2], FUN=function(chr) {
	
				x  <- start(bgranges)[as.logical(seqnames(bgranges) == chr)]
				bg <- score(bgranges)[as.logical(seqnames(bgranges) == chr)]


				list(
					chr = chr,
					scale = scale,
					x          = x,
					bg         = bg,
					n.insertions       = chr.info[chr, 3],
					n.insertions.total = n.insertions.total,
					n.bgsites          = chr.info[chr, 2],
					n.bgsites.total    = n.bgsites.total,
					n.bp               = chr.info[chr, 1],
					n.bp.total         = n.bp.total
				)
			})

			names(chr.objects) <- chromosomes
			

			chr.objects

		})

		names(scale.objects) <- .formatScales(scales)


		kse.distributions <- mclapply(scales, mc.cores=n.cores[1], FUN=function(scale) {
			message('cimplr.reference() - compute kse distribution scale = ', scale)
			new('KSEDistribution', scale=scale, D=D, kNormal=30, verbose=FALSE)
		})

		names(kse.distributions) <- .formatScales(scales)

		list(
			input    = cimplr.object$input,
			data     = scale.objects,
			kse.distributions = kse.distributions
		)
	}) # end 'with'

}


cimplr.convolve <- function(cimplr.object) {
	message('cimplr.convolve()')

	with(cimplr.object$input, {	

		#
		# 1. KSE convolution and each scale
		#
		scale.objects <- mclapply(cimplr.object$data, mc.cores=n.cores[1], FUN=function(scale.object) {
			
			scale <- scale.object[[1]]$scale
			kse_dist <- cimplr.object$kse.distributions[[.formatScales(scale)]]


			#
			# Calculation kse objects for each chromosome
			#
			scale.object <- mclapply(scale.object, mc.cores=n.cores[2], FUN=function(chr.object) {
				message('cimplr.convolve() - calculcate KSE (scale = ', scale, ', chr = ', chr.object$chr, ')')
				calc.kse(chr.object, kse_dist, insertions[[chr.object$chr]], pkse.method, n.insertions.total, n.bgsites.total, alpha.level, p.adjust.method, p.adjust.n.method, mtf, test.chromosomes.seperately)
			})

			
			
      
      
#			n.peaks.total <- sum(sapply(scale.object, '[[', 'n.peaks'))

#      
#			scale.object <- mclapply(scale.object, mc.cores=n.cores[2], FUN=function(chr.object) {
#				if (is.null(chr.object[['th']])) {
#					message('cimplr.convolve() - calculcate threshold (scale = ', scale, ', chr = ', chr.object$chr, ')')
#					#chr.object <- calc.threshold(chr.object, kse_dist, alpha=chr.object$adjusted.alpha.level)
#					chr.object <- calc.threshold(chr.object, kse_dist)
#				}
#				
#				message('cimplr.convolve() - calculcate ratio (scale = ', scale, ', chr = ', chr.object$chr, ')')
#				chr.object$ratio <- chr.object$inskse / chr.object$th
#				chr.object$significant <- chr.object$inskse > chr.object$th
#				
#				chr.object
#			})


			scale.object
		})
    
		
    list(
			input = cimplr.object$input,
			data  = scale.objects,
			kse.distributions = cimplr.object$kse.distributions
		)

	}) # end 'with'
}


cimplr.threshold <- function(cimplr.object, alpha) {
  message('cimplr.threshold()')

  # step 1: calculate the FDR rate
  fdr.pval <- calc.fdr.pval(cimplr.object)
  
  
  with(cimplr.object$input, {
    scale.objects <- mclapply(cimplr.object$data, mc.cores=n.cores[1], FUN=function(scale.object) {
      
      scale <- scale.object[[1]]$scale
      kse_dist <- cimplr.object$kse.distributions[[.formatScales(scale)]]
      
      
  		scale.object <- mclapply(scale.object, mc.cores=n.cores[2], FUN=function(chr.object) {
				if (is.null(chr.object[['th']])) {
					message('cimplr.convolve() - calculcate threshold (scale = ', scale, ', chr = ', chr.object$chr, ')')
					chr.object <- calc.threshold(chr.object, kse_dist, alpha=fdr.pval)
				}
				
				message('cimplr.convolve() - calculcate ratio (scale = ', scale, ', chr = ', chr.object$chr, ')')
				chr.object$ratio <- chr.object$inskse / chr.object$th
				chr.object$significant <- chr.object$inskse > chr.object$th
				
				chr.object
			})
    
    
      scale.object
    })
    
    
    list(
      input = cimplr.object$input,
      data  = scale.objects,
      kse.distributions = cimplr.object$kse.distributions,
      fdr.pval = fdr.pval
    )
    
  })
}


cimplr.call <- function(cimplr.object) {

	message('cimplr.call()')

	with(cimplr.object$input, {
    
    
    
    

		scale.objects <- mclapply(cimplr.object$data, mc.cores=n.cores[1], FUN=function(scale.object) {

			#
			# Calculate the CIS calling for each chromosome
			#
			scale.object <- mclapply(scale.object, mc.cores=n.cores[2], FUN=function(chr.object) {
				message('cimplr.call() - calculcate CIS (scale = ', chr.object$scale, ', chr = ', chr.object$chr, ')')
				calc.cis(chr.object)
			})



			scale.object
		}) # end scale space


		#
		# merge scale space, plopping
		#
		message('cimplr.call() - merge CIS list')
		all_cises <- do.call('rbind', lapply(scale.objects, function(cimplr) do.call('rbind', lapply(cimplr, '[[', 'cises'))))
		# fix rownames, somehow rbind screws them
		if (!is.null(all_cises)) {
			rownames(all_cises) <- paste('CIS', substring(all_cises$chromosome, 4), ':', all_cises$peak.location, '_', .formatScales(all_cises$scale), sep='')
		}



		collapsed_cises <- do.call('rbind', mclapply(chromosomes, mc.cores=n.cores[2], FUN=function(chr) {
			message('cimplr.call() - call CIS across scales (chr = ', chr, ')')
			collapse.cis(scale.objects, scales, chr)
		}))



		if (!is.null(collapsed_cises)) {
			rownames(collapsed_cises) <- paste('CIS', substring(collapsed_cises$chromosome, 4), ':', collapsed_cises$peak.location, '_', .formatScales(collapsed_cises$scale), sep='')
		}





		if ( !is.null(genes.file) ) {
			

			message('cimplr.call() - annotate CIS')
			load(genes.file) # contains 'genes' object.

			all_cises       <- annotate.cis(all_cises, genes)
			collapsed_cises <- annotate.cis(collapsed_cises, genes)


		}

		output <- list(
			all.cises = all_cises,
			collapsed.cises = collapsed_cises
		)

		list(
			input  = cimplr.object$input,
			data   = scale.objects,
			kse.distributions = cimplr.object$kse.distributions,
			output = output
		)
	}) # end 'with'
}


#
# Does the kernel convolution
#
calc.kse <- function(chr.object, kse_dist, locs, pkse.method, n.insertions.total, n.bgsites.total, alpha.level, p.adjust.method, p.adjust.n.method, mtf, test.chromosomes.seperately) {
	with(chr.object, {
		
		# insertion density
		ftn     <- list(from=min(x), to=max(x), n=length(x))
		insdens <- density(locs, bw=scale, from=ftn$from, to=ftn$to, n=ftn$n)
		inskse  <- insdens$y * length(locs) / dnorm(0, 0, sd=scale)

		# n.peaks
		peak.idx <- whichLocalMaxima(inskse)
		peaks    <- x[peak.idx]
		n.peaks  <- length(peaks)

		# MT correction
		if (test.chromosomes.seperately) {
			#'n.x', 'n.bgsites', 'n.peaks', 'n.bgsites/scale', 'none'

			p.adjust.n <- switch(p.adjust.n.method,
				n.x               = length(x), 
				n.bgsites         = n.bgsites,
				n.peaks           = n.peaks,
				'n.bgsites/scale' = n.bgsites / kse_dist@scale
			)
		} else {
			p.adjust.n <- switch(p.adjust.n.method,
				n.x               = length(x), 
				n.bgsites         = n.bgsites.total,
				n.peaks           = n.peaks.total,
				'n.bgsites/scale' = n.bgsites.total / kse_dist@scale
			)
		}

		adjusted.alpha.level <- alpha.level / p.adjust.n / mtf

		p  <- rep(NA, length(x))
		idx <- bg > 0
    
		#p[idx] <- 1- pkse(kse_dist, x=x[idx], n=bg[idx], p=n.insertions.total / n.bgsites.total)
		p[idx] <- 1 - pkse(kse_dist, x=inskse[idx], n=bg[idx], p=n.insertions.total / n.bgsites.total)
		
    
		c(
			chr.object,
			list(
				locs        = locs,
				n.insertions.total = n.insertions.total,
				inskse      = inskse,
        p           = p,
				peaks       = peaks,
				n.peaks     = n.peaks,
				pkse.method = pkse.method,
				p.adjust.n           = p.adjust.n,
				p.adjust.method      = p.adjust.method,
				p.adjust.n.method    = p.adjust.n.method,
				mtf                  = mtf,
				alpha.level          = alpha.level,
				adjusted.alpha.level = adjusted.alpha.level
			)
		)
	})
}









calc.fdr.pval <- function(cimplr.object) {
  message('calc.fdr.pval()')
  
  with(cimplr.object$input, {
    # contains all pvalues across all scales and all chromosomes
    all.pvals <- sapply(cimplr.object$data, function(scale.object) sapply(scale.object, '[[', 'p'))

    
    p.sorted <- sort(all.pvals, na.last=NA) # note NA's are removed. p is NA when the background is < 0
    
    m <- length(p.sorted)
    k <- 1:m
    
    th <- k * alpha.level / m
    
    fdr.pval.idx <- which(!p.sorted < th)[1]
    
    fdr.pval <- p.sorted[fdr.pval.idx]
        
    fdr.pval
  })
  
}


calc.threshold <- function(chr.object, kse_dist, alpha) {
	with(chr.object, {

		th <- rep(NA, length(x))

    idx <- bg > 0
	
		#if (p.adjust.method == 'bonferroni') {
#
#			th[idx] <- kseThreshold(kse_dist, n=bg[idx], p=n.insertions.total / n.bgsites.total, max.kse=ceiling(max(inskse[idx])) * 2, alpha=adjusted.alpha.level, max.k=120, n.unique.p=2000, n.x=1000, summed.weights.margin = 1e-15)

#		} else if (p.adjust.method == 'fdr') {

      th[idx] <- kseThreshold(kse_dist, n=bg[idx], p=n.insertions.total / n.bgsites.total, max.kse=ceiling(max(inskse[idx])) * 2, alpha=alpha, max.k=120, n.unique.p=2000, n.x=1000, summed.weights.margin = 1e-15)
#		}

		c(
			chr.object,
			list(
				th = th
			)
		)
	}) # end 'with'
}



#
# Calculates CIS for each scale
#
calc.cis <- function(chr.object) {
	with(chr.object, {



		# get CIS as segments

		segs <- segments(significant)

		if (segs$n == 0) {
			cises <- NULL
		} else {

			# peak is defined as the 
			peak.idx <- mapply(FUN=function(start, end) {
				which.min(ratio[start:end]) + start - 1
			}, segs$start.pos, segs$end.pos)


			peak.location   <- x[peak.idx]
			peak.kse.value  <- inskse[peak.idx]
			peak.bg         <- bg[peak.idx]
			peak.ratio      <- ratio[peak.idx]


			start <- x[segs$start.pos]
			#end   <- pmin(x[segs$end.pos], sl)
			end <- pmin(x[segs$end.pos])


			n.insertions.within.3.sigma <- mapply(FUN=function(start, end) {
				sum( locs >=  (start - 3 * scale) & locs <=  (end + 3 * scale) )
			}, start, end)




			cises <- data.frame(
				chromosome      = chr,
				start           = start,
				end             = end,
				width           = end - start,
				peak.location   = peak.location,
				peak.kse.value  = peak.kse.value,
				peak.bg         = peak.bg,
				peak.ratio      = peak.ratio,
				n.insertions.within.3.sigma = n.insertions.within.3.sigma,
				scale           = scale,
				alpha.level     = alpha.level,
				p.adjust.n      = p.adjust.n,
				p.adjust.method = p.adjust.method,
				p.adjust.n.method = p.adjust.n.method,
				mtf             = mtf,
				adjusted.alpha.level = adjusted.alpha.level,
				stringsAsFactors=FALSE
			)

			rownames(cises) <- paste('CIS', substring(chr, 4), ':', peak.location, '_', .formatScales(scale), sep='')

		}


		c(
			chr.object,
			list(
#				p.adjust.n           = p.adjust.n,
#				p.adjust.method      = p.adjust.method,
#				p.adjust.n.method    = p.adjust.n.method,
#				mtf                  = mtf,
#				adjusted.alpha.level = adjusted.alpha.level,
#				significant = significant,
#				th          = th,
				cises       = cises
			)
		)

	}) # end 'with'
}





#
# Does the cross-scale CIS calling given a chromosome.
#
collapse.cis <- function(scale.objects, scales, chr) {

	# Calculate the min.p.value CIS
	x <- scale.objects[[1]][[chr]]$x
	#locs <- scale.objects[[1]][[chr]]$locs

	#p.adjusted <- sapply(scale.objects, function(cimplr) cimplr[[chr]]$p.adjusted)
	
	# If one scale only
	#if (length(scale.objects) == 1) {
	#	p.adjusted <- matrix(p.adjusted, ncol=1)
	#}


	#p.adjusted[p.adjusted==0] <- min(p.adjusted[p.adjusted!=0])
	
	#min.p.adjusted <- apply(p.adjusted, 1, min)
	#significant <- min.p.adjusted < alpha.level

	# From now on we take the max kse/threshold ratio...

	ratios <- sapply(scale.objects, function(cimplr) cimplr[[chr]]$ratio)
	#suppressWarnings (max_ratios <- apply(ratios, 1, max, na.rm=TRUE) )
	
	#significant <- max_ratios > 1

	significant.mat <- sapply(scale.objects, function(cimplr) cimplr[[chr]]$significant)
	significant <- apply(significant.mat, 1, any, na.rm=TRUE)
	significant[is.na(significant)] <- FALSE

	# get CIS
	segs <- segments(significant)


	if (segs$n == 0) {
		cises <- NULL
	} else {


		cises <- do.call('rbind', mapply(SIMPLIFY=FALSE, FUN=function(seg.start, seg.end) {
			

			## Convert it to a raster object

			#sub.p.adjusted <- p.adjusted[seg.start:seg.end, length(scales):1, drop=FALSE]
			#r <- raster(t(sub.p.adjusted))

			# take the th / inske kse ratio
			#sub.ratio <- th[seg.start:seg.end, length(scales):1, drop=FALSE] / inskse[seg.start:seg.end, length(scales):1, drop=FALSE]

			sub.sign <- significant.mat[seg.start:seg.end, length(scales):1, drop=FALSE]
			r_significant <- raster(t(sub.sign))
			extent(r_significant) <- extent(c(seg.start-1, seg.end, 0, length(scales)) + .5)


			sub.ratios <- ratios[seg.start:seg.end, length(scales):1, drop=FALSE]
			sub.ratios[sub.ratios < 1] <- NA
			r <- raster(t(sub.ratios))

			#extent(r) <- extent(c(0, ncol(sub.p.adjusted), 0, nrow(sub.p.adjusted)) + 0.5)
			extent(r) <- extent(c(seg.start-1, seg.end, 0, length(scales)) + .5)

			## Find the min value within the 9x9-cell neighborhood of each cell
			max.func <- function(x) {
				if (all(is.na(x))) {
					-Inf
				} else if ( sum(x == max(x, na.rm=TRUE), na.rm=TRUE) > 1 ) {
					-Inf
				} else {
					max(x, na.rm=TRUE)
				}

			}

#			localmax <- focal(r, w = matrix(rep(1, 21*21), nrow=21, ncol=21), fun = max.func, pad=TRUE, padValue=1)
			n <- 3
			localmax <- focal(r, w = matrix(rep(1, n*n), nrow=n, ncol=n), fun = max.func, pad=TRUE, padValue=1)

			## Does each cell have the maximum value in its neighborhood?
			r2 <- r==localmax


			## Get x-y coordinates of those cells that are local minima
			maxXY <- xyFromCell(r2, Which(r2==1, cells=TRUE))

			if (nrow(maxXY) == 0) {
				stop('No local maxima found in scale space!')
			}


			# two or more local maxima found at the same x position
			# this probably caused by the irri
#			if (nrow(maxXY) > 1 & all(maxXY[, 'x'] == maxXY[1, 'x'])) {				
        # > more than 2, call only one CIS, the one with the highest ratio
#        browser()
        
        
#			} else {
        
        # 
        
#			}

			## For each maxXY (a maximimam in the ratio space), select the CIS
			cises <- unique(do.call('rbind', lapply(1:nrow(maxXY), function(i) {
				
				peak.location <- x[maxXY[i, 'x']]
				cises_at_selected_scale <- scale.objects[[maxXY[i, 'y']]][[chr]]$cises

				idx <- cises_at_selected_scale$start <= peak.location & cises_at_selected_scale$end >= peak.location
				if (sum(idx) != 1) {
					plot(r)
					plot(r2)
					plot(r_significant)
					dev.off()
					browser()
				}

				stopifnot( sum(idx) == 1 )


        
				cises_at_selected_scale[idx, ]
			})))


      # select the CIS at the global max of the cloud.
      cises[which.max(cises$peak.ratio), ]


		}, segs$start.pos, segs$end.pos))
	}
	cises
}




annotate.cis <- function(cises, genes) {

	geneIdentifiers <- attr(genes, 'geneIdentifiers')
	genes <- split(genes, genes$chromosome_name)

	n.cis <- nrow(cises)

	for(id in geneIdentifiers) {
		cises[paste('associated_', id, sep='')] <- rep('', n.cis)
		cises[paste('other_', id, sep='')] <- rep('', n.cis)
	}

	for (i in 1:n.cis) {
		ags <- .associateGenes(cises$start[i], cises$end[i], cises$peak.location[i], genes[[cises$chromosome[i]]])
		for(id in geneIdentifiers) {
			
			a <- genes[[cises$chromosome[i]]][ags$associated, id]
			a <- a[a!='']
			
			o <- genes[[cises$chromosome[i]]][ags$other, id]
			o <- o[o!='']

			cises[i, paste('associated_', id, sep='')] <- paste(a, collapse='|')
			cises[i, paste('other_', id, sep='')] <- paste(o, collapse='|')
		}
	}

	cises
}
























#
# Exports a cimplr object
#
export.cimplr <- function(cimplr.object, output.dir, include.scales=TRUE) {

	# create the output dir
	if (!file.exists(output.dir)) {
		#dir.create(output.dir, showWarnings=FALSE, recursive=TRUE)
		stop(paste("output.dir:", output.dir, "does not exist."))
	}

	with(c(cimplr.object$input, cimplr.object['data'], cimplr.object$output), {

		n.bgsites <- chr.info[[2]]
		names(n.bgsites) <- chromosomes

		if (include.scales) {
			cat('exporting scale profiles...')
			for (scale in scales) {
				dir.create(paste(output.dir, '/scales/', .formatScales(scale), sep=''), showWarnings=FALSE, recursive=TRUE)
				cis_file     <- paste(output.dir, '/scales/', .formatScales(scale), '/cises.txt', sep='')


				# Make figures...

				## Make plots
				cimplr <- data[[.formatScales(scale)]]

				n.peaks.total <- sum(sapply(cimplr, '[[', 'n.peaks'))
				mclapply(names(cimplr), mc.cores=n.cores[2], FUN=function(chr) {

					figures_file <- paste(output.dir, '/scales/', .formatScales(scale), '/figures.', chr, '.png', sep='')
					png(file=figures_file, width=2000, height=2500)
					layout(matrix(1:4, ncol=1))
				
					plot(cimplr[[chr]]$x, cimplr[[chr]]$inskse, type='l', main='insertion kse')
					lines(cimplr[[chr]]$x, cimplr[[chr]]$th, col='red')
					legend('topleft', legend=c(
						paste(chr, ', scale =', scale),
						paste('# insertions = ', length(insertions[[chr]]), '/', n.insertions.total, ' = ', signif(length(insertions[[chr]])/ n.insertions.total, 3), sep=''),
						paste('# ', pattern, ' = ', n.bgsites[chr], '/', n.bgsites.total, ' = ', signif(n.bgsites[chr]/ n.bgsites.total, 3), sep='')
					))
					plot(cimplr[[chr]]$x, cimplr[[chr]]$bg , type='l', main='background density')
					plot(cimplr[[chr]]$x, cimplr[[chr]]$bg * n.bgsites.total / dnorm(0, 0, sd=scale) , type='l', main='background kse count')

#					plot(cimplr[[chr]]$x, -log10(cimplr[[chr]]$p), type='l')
#					abline(h=-log10(alpha.level / n.peaks.total), lwd=2, col='blue')
#					abline(h=-log10(alpha.level / n.bgsites.total), lwd=2, col='red')

					plot(cimplr[[chr]]$x, cimplr[[chr]]$ratio, type='l')
					abline(h= 1, lwd=2, col='red')

#					hist(cimplr[[chr]]$p, breaks=100)
#					plot(cimplr[[chr]]$x, -log10(cimplr[[chr]]$p.adjusted), type='l', ylim=c(0, 10))
#					abline(h=-log10(alpha.level), lwd=2, col='red')
				#	abline(h=-log10(.05), col='red')
				#	abline(h=-log10(.01), col='red')
				#	abline(h=-log10(.001), col='red')
					dev.off()

				})


				cises <- do.call('rbind', lapply(cimplr, '[[', 'cises'))

				# fix rownames, somehow rbind screws them
				if (!is.null(cises)) {
					rownames(cises) <- paste('CIS', substring(cises$chromosome, 4), ':', cises$peak.location, '_', .formatScales(cises$scale), sep='')
				}
				
				write.table(cises, file=cis_file, row.names=FALSE, quote=FALSE, sep='\t')
			}
			cat('done!\n')
		}



		if (!is.null(collapsed.cises)) {
			collapsed_cises <- split(collapsed.cises, collapsed.cises$chromosome)
		} else {
			collapsed_cises <- NULL
		}

		cat('exporting cross-scale profiles: ')
		dir.create(paste(output.dir, '/scale_space', sep=''), showWarnings=FALSE, recursive=TRUE)

		mclapply(chromosomes, mc.cores=n.cores[2], FUN=function(chr) {
			cat(chr, '... ', sep='')
			scale_space_file <- paste(output.dir, '/scale_space/figures.', chr,'.png', sep='')
		
#			png(file=scale_space_file, width=2000, height=1000)

	#		x <- sapply(scale_space, function(cimplr) cimplr[[chr]]$x)
			x <- data[[1]][[chr]]$x
#			p <- sapply(data, function(cimplr) cimplr[[chr]]$p)


#			p.adjusted <- sapply(data, function(cimplr) cimplr[[chr]]$p.adjusted)
#			p.adjusted[p.adjusted==0] <- min(p.adjusted[p.adjusted!=0])
			
#			min.p.adjusted <- apply(p.adjusted, 1, min)

#			th <- sapply(data, function(cimplr) cimplr[[chr]]$th)
			inskse <- sapply(data, function(cimplr) cimplr[[chr]]$inskse)

			ratio <- sapply(data, function(cimplr) cimplr[[chr]]$ratio)
			#th.inskes.diff <- inskse - th
#			ratio <- inskse / th

			# overall
#			matplot(x, -log10(p.adjusted), lty=1, lwd=1, type='l', col=rainbow(ncol(p.adjusted)))
#			lines(x, -log10(min.p.adjusted), type='l', lty=2, lwd=2, col='darkblue')
#			legend('topright', names(data), lty=1, lwd=1, col=rainbow(length(data)))
#			abline(h=-log10(alpha.level), lty=2, col='red')
#			#abline(h=-log10(alpha.level / n.peaks.total), lty=2, lwd=2, col='blue')
#			abline(h=-log10(alpha.level / n.bgsites.total), lty=2, lwd=2, col='red')

#			dev.off()



#			p.adjusted[p.adjusted >= alpha.level] <- NA
#			p[p >= alpha.level] <- NA

			# for each cross scale cis
			if (!is.null(collapsed_cises[[chr]])) {
				
				# plot 1: kse's around CIS
				scale_space_cises_file <- paste(output.dir, '/scale_space/figures.', chr, '.cises.pdf', sep='')
				
				pdf(file=scale_space_cises_file, width=20, height=10)
				layout(matrix(1:10, c(2,5), byrow=FALSE))
				for (i in 1:nrow(collapsed_cises[[chr]])) {
					margin <- 200e3
					xlim <- c(collapsed_cises[[chr]][i, 'start'] - margin, collapsed_cises[[chr]][i, 'end']+margin)
					idx <- x >= xlim[1] & x <= xlim[2]
					#browser()
#					matplot(x[idx], -log10(p.adjusted[idx, ]), lty=1, lwd=1, type='l', col=rainbow(ncol(p.adjusted)), main='kse')
					#lines(x[idx], -log10(min.p.adjusted[idx]), type='l', lty=2, lwd=2, col='darkblue')
#					abline(h=-log10(alpha.level), lty=2, col='red')
					#abline(h=-log10(alpha.level / n.peaks.total), lty=2, lwd=2, col='blue')
#					abline(h=-log10(alpha.level / n.bgsites.total), lty=2, lwd=2, col='red')
					#legend('topright', names(data), lty=1, lwd=1, col=rainbow(length(data)))


					# ratio image
					matplot(x[idx], ratio[idx, ], lty=1, lwd=1, type='l', col=rainbow(length(scales)), main='kse / threshold')
					abline(h=1, lty=2, lwd=1, col='red')

					locs <- insertions[[chr]][insertions[[chr]] >= xlim[1] & insertions[[chr]] <= xlim[2]]

					#image(x=x[idx], y=log10(scales), z=-log10(p.adjusted[idx, ]), main='p.adjusted')


#					image(x=x[idx], y=log10(scales), z=-log10(p[idx, ]), main='p')


					# ratio
					image(x=x[idx], y=log10(scales), z=ratio[idx, ], main='kse / threshold')

					abline(v=locs, lty=2, col='blue', lwd=1)
					lines(
						x=c(collapsed_cises[[chr]][i, 'start'], collapsed_cises[[chr]][i, 'end']),
						y=c(log10(collapsed_cises[[chr]][i, 'scale']), log10(collapsed_cises[[chr]][i, 'scale'])),
						lty=1, col='darkgreen', lwd=3
					)
				}
				dev.off()

				# plot2:

			}

			#browser()

			#if (length(scales) > 1) {
			#	rh <- min(scales[-1] - scales[-length(scales)])  # rectangle height
			#} else {
			#	rh <- 1
			#}
			

			#plot(1, type='n', xlim=range(scale_space[[length(scale_space)]][[chr]]$x), ylim=range(scales) + c(-rh/2, rh/2), xlab='Genomic position', ylab='Scale', main='Scale space')
			
			#for(cimplr in scale_space) {
			#	if (!is.null(cimplr[[chr]]$cises)) {
			#		for (i in 1:nrow(cimplr[[chr]]$cises)) {
			#			rect(xleft=cimplr[[chr]]$cises[i, 'start'], ybottom=cimplr[[chr]]$scale-rh/2, xright=cimplr[[chr]]$cises[i, 'end'], ytop=cimplr[[chr]]$scale+rh/2, col='green', border=NA)
			#		}
			#	}
			#}

		})
		cat('done!\n')



		collapsed_cis_file <- paste(output.dir, '/collapsed_cises.csv', sep='')
		cat('writing ', collapsed_cis_file, '...', sep='')
		write.csv(collapsed.cises, file=collapsed_cis_file, row.names=TRUE)
		cat('done!\n')

		all_cis_file     <- paste(output.dir, '/all_cises.csv', sep='')
		cat('writing ', all_cis_file, '...', sep='')
		write.csv(all.cises, file=all_cis_file, row.names=TRUE)
		cat('done!\n')

	})

}


