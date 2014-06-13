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

cimplr <- function(
  insertions,

  exclude.chromosomes=c(),
  
  scales=log.seq(500, 500000, 50),
  
  biasmap = 'biasmaps_chr10/bsgenomes-mm9-TA-stepsize100-D4',
  biasmap.type = c('count', 'density'),

  alpha=.05,
  p.adjust.method=c("fdr", "BY", "bonferroni", "none"),

  kNormal = 30,
	
  n.cores = 4
  
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


	if (alpha <= 0 | alpha >= 1) {
		stop("Need to provide 0 < 'alpha.level' < 1.")
	}


	biasmap.type <- match.arg(biasmap.type)
	p.adjust.method <- match.arg(p.adjust.method)

	# determine the available chromosomes
	chromosomes <- unique(as.character(insertions$chr))
	chromosomes <- chromosomes[!chromosomes %in% exclude.chromosomes]

	# make 'insertions' a list of insertions
	insertions.org <- insertions
	insertions  <- insertions[insertions$chr %in% chromosomes, ]
	insertions  <- split(insertions$location, insertions$chr)[chromosomes]
	insertions  <- insertions[chromosomes]

	# load the biasmap info
	chr.info <- loadBiasmapInfo(biasmap)
	chr.info <- chr.info[chromosomes, ]
	chr.info$n.insertions=sapply(insertions, length)

	n.insertions.total <- sum(as.numeric(chr.info[[3]]))
	n.bgsites.total <- sum(as.numeric(chr.info[[2]]))
	n.bp.total <- sum(as.numeric(chr.info[[1]]))
  
	D <- as.numeric(substring(rev(strsplit(biasmap, '-')[[1]])[1], 2))
  
  # bundle all inputs
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


# replaces the instertions, leaving the chromosome intact.
cimplr.replace_insertions <- function(cimplr.object, new.insertions) {

	message('cimplr.replace_insertions()')
	with(cimplr.object$input, {

		insertions <- new.insertions
		insertions.org <- insertions
		insertions  <- insertions[insertions$chr %in% chromosomes, ]
		insertions  <- split(insertions$location, insertions$chr)[chromosomes]
		insertions  <- insertions[chromosomes]

		# load the biasmap info
		chr.info <- loadBiasmapInfo(biasmap)
		chr.info <- chr.info[chromosomes, ]
		chr.info$n.insertions=sapply(insertions, length)
		
		n.insertions.total <- sum(chr.info[[3]])
		n.bgsites.total <- sum(chr.info[[2]])
		n.bp.total <- sum(chr.info[[1]])
		
		D <- as.numeric(substring(rev(strsplit(biasmap, '-')[[1]])[1], 2))

		input <- as.list(environment())


		##

		## clear the fields in data

		new.data <- lapply(cimplr.object$data, function(chr.object) {
				if ('th' %in% names(chr.object)) {
          
          if (p.adjust.method == 'fdr' | p.adjust.method == 'BY') {
            message("cimplr.replace_insertions() - removing threshold because it depends on the data (in case of 'fdr' and 'BY')")
            new.chr.object <- chr.object[c("chr", "scale", "x", "bg", "n.insertions", "n.insertions.total", "n.bgsites", "n.bgsites.total", "n.bp", "n.bp.total", "n.x")]
          } else {
            new.chr.object <- chr.object[c("chr", "scale", "x", "bg", "n.insertions", "n.insertions.total", "n.bgsites", "n.bgsites.total", "n.bp", "n.bp.total", "n.x", "th", "corrected.alpha")]
          }
				} else {
				  new.chr.object <- chr.object[c("chr", "scale", "x", "bg", "n.insertions", "n.insertions.total", "n.bgsites", "n.bgsites.total", "n.bp", "n.bp.total", "n.x")]
				}
				new.chr.object$locs <- insertions[[chr.object$chr]]
				new.chr.object
    })

		list(
			input = input,
			data  = new.data,
			kse.distributions = cimplr.object$kse.distributions
		)
	})
}


cimplr.biasmap <- function(cimplr.object) {
  

	with(cimplr.object$input, {	

	  cluster <- makeForkCluster(n.cores)
	  
	  message('cimplr.biasmap() - load biasmaps')

    new.data <- clusterMap(
	    cl=cluster,
	    fun=calc.biasmap,
	    
	    rep(chromosomes, each=length(scales)),
	    rep(scales, length(chromosomes)),
	    
	    MoreArgs = list(
	      input = cimplr.object$input
	    )
	  )
	  
	  names(new.data) <- paste(rep(chromosomes, each=length(scales)), as.integer(round(rep(scales, length(chromosomes)))), sep='.')

    # cluster 2 for KSE distributions
	  
    message('cimplr.biasmap() - load KSE distributions')
	  
	  
	  kse.distributions <- clusterMap(
      cl=cluster,
      fun=function(scale, D, kNormal) {
        new('KSEDistribution', scale=scale, D=D, kNormal=kNormal, verbose=FALSE)
      },
      
      scales,
      MoreArgs = list(
        D = D,
        kNormal=kNormal
      )
    )
    
	  stopCluster(cluster)
	  
		names(kse.distributions) <- .formatScales(scales)
		
		list(
			input    = cimplr.object$input,
			data     = new.data,
			kse.distributions = kse.distributions
		)
	}) # end 'with'

}



cimplr.convolve <- function(cimplr.object) {
  message('cimplr.convolve()')
  
  with(cimplr.object$input, {	

    cluster <- makePSOCKcluster(n.cores)
    clusterExport(cl=cluster, c('whichLocalMaxima', 'pkse', 'pkseCond'))
    
    new.data <- clusterMap(
        cl=cluster,
        fun=calc.kse,
        
        cimplr.object$data,                # these are scale / chr
        cimplr.object$kse.distributions,
        
        MoreArgs = list(
          pkse.method=biasmap.type,
          n.insertions.total=n.insertions.total,
          n.bgsites.total=n.bgsites.total,
          alpha=alpha,
          p.adjust.method=p.adjust.method
        )
    )
    
    stopCluster(cluster)
    
    list(
      input = cimplr.object$input,
      data  = new.data,
      kse.distributions = cimplr.object$kse.distributions
    )
    
  }) # end 'with'
}



cimplr.ratios <- function(cimplr.object) {
  message('cimplr.ratios() - calculate corrected alpha')

  # step 1: calculate the FDR rate
  corrected.alpha <- calc.corrected.alpha(cimplr.object)
    
  message('cimplr.ratios() - calculate ratios')
  with(cimplr.object$input, {
    
    cluster <- makePSOCKcluster(n.cores)
    clusterExport(cl=cluster, c('kseThreshold', 'pkseCond'))
    
    new.data <- clusterMap(
      cl=cluster,
      fun = calc.ratios,

      cimplr.object$data,                # these are scale / chr
      cimplr.object$kse.distributions,
      
      .scheduling = 'dynamic',

      MoreArgs = list(
        corrected.alpha = corrected.alpha,
        update.threshold = TRUE
      )
    )
    
    stopCluster(cluster)
    
    list(
      input = cimplr.object$input,
      data  = new.data,
      kse.distributions = cimplr.object$kse.distributions
    )
    
  })
}




cimplr.call <- function(cimplr.object) {


	with(cimplr.object$input, {
    
    
    
	  cluster <- makePSOCKcluster(n.cores)
	  clusterExport(cl=cluster, c('segments', '.formatScales', 'GRanges', 'IRanges', 'Rle'))

    message('cimplr.call() - call CIS')
	  
	  new.data <- clusterMap(
	    cl=cluster,
	    fun = calc.cis,
	    
	    cimplr.object$data,                # these are scale / chr
	    
	    MoreArgs = list(
	    )
	  )
	  
	  stopCluster(cluster)    
    
		# merge all cises into 1 GRanges
    all_cises <- unlist(do.call('GRangesList', lapply(new.data, '[[', 'cises')), use.names=FALSE)
		
  
		message('cimplr.call() - call cross-scale CIS')
        
    cluster <- makePSOCKcluster(n.cores)
    clusterExport(cl=cluster, c('segments', 'raster', 'extent', 'extent<-', 'focal', 'xyFromCell', 'Which', 'GRangesList'))

    ### testing only
		scalechr <- lapply(chromosomes, function(chr) {
		  idx <- sapply(new.data, function(obj) obj$chr == chr)
		  new.data[idx]
		})
    ### testing only
		collapse.cis(scalechr[[1]])
		collapsed_cises <- unlist(do.call('GRangesList', clusterMap(
		  cl=cluster,
		  fun = collapse.cis,
		  
		  lapply(chromosomes, function(chr) {
		    idx <- sapply(new.data, function(obj) obj$chr == chr)
		    new.data[idx]
		  }),                # these are scalechr objects grouped by chromosome.
		  
		  MoreArgs = list(
		  )
		)), use.names=FALSE)
		
		stopCluster(cluster)    

    output <- list(
			all.cises = all_cises,
			collapsed.cises = collapsed_cises
		)

		list(
			input  = cimplr.object$input,
			data   = new.data,
			kse.distributions = cimplr.object$kse.distributions,
			output = output
		)
	}) # end 'with'
}





