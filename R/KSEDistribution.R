# KSEDistribution.R
#
# Description: An S4 object for KSE distribution.
# Author: Jelle ten Hoeve
#


#
# Initializor of the KSEDistrubition object
#   D: the number of standard deviations under the kernel devided by 2
#   nBins: the number of bins to use for calculation of the histogram of the
#     amplitude distribution
#   k.normal: After which 'k' we assume normality of the conditional KSE probability.
#

setMethod('initialize', 'KSEDistribution', function(.Object, scale, D=4, nBins=1000, kNormal=120, verbose=FALSE) {
	.Object@nBins <- nBins
	.Object@D <- D
	.Object@kNormal <- kNormal
	
	.Object@scale <- scale
	
	.Object@kernelX <- seq(-1 * scale * D, scale * D, 1)
	.Object@kernelY <- dnorm(.Object@kernelX, mean=0, sd=scale) / dnorm(0, mean=0, sd=scale)
	
	
	breaks <- seq(0, 1, length.out=.Object@nBins+1)
	h <- hist(.Object@kernelY, breaks=breaks, plot=FALSE)
	
	breaks  <- h$breaks
	density <- h$density
	mids    <- h$mids
	
	mu  <- sum(mids * density / .Object@nBins)
	var <- sum( (mids - mu)^2 * density / .Object@nBins)
	
	

	# For N=0 the histogram is the delta dirac function.
	.Object@histograms[[1]] <- list(
		N       = 1,
		breaks  = breaks,
		density = density,
		cumdens = cumsum(density / .Object@nBins),
		mids    = mids,
		mu      = mu,
		var     = var
	)
	
	if (verbose) {		
		cat('Convolving histograms...\n')
	}
	
	N <- 2
	while (N <= .Object@kNormal) {
		prev <- .Object@histograms[[N-1]] # note N is used for indexing!
		


		# uses FFT
		new.density <- c(0, stats::convolve(
			prev$density / .Object@nBins,
			rev(density / .Object@nBins),
			type='o'
		)) * .Object@nBins

		# no FFT
		new.density <- filter(
			x=c(rep(0, .Object@nBins), prev$density/.Object@nBins, rep(0, .Object@nBins)), 
			filter=density/.Object@nBins,
			side=1, method='convolution')[.Object@nBins:((N+1) * .Object@nBins - 1)] * .Object@nBins
	
		new.breaks  <- c(prev$breaks, breaks[-1]+N-1)
		new.mids    <- c(prev$mids, mids + N-1) 
		
		.Object@histograms[[N]] <- list(
			N       = N,
			breaks  = new.breaks,
			density = new.density,
			cumdens = cumsum(new.density) / .Object@nBins,
			mids    = new.mids,
			mu      = N * mu,
			var     = N * var
		)
		
		if (verbose) {
			cat('\r', N, '/', .Object@kNormal, sep='')
		}
		
		N <- N + 1
	}

	if (verbose) {
		cat('\ndone!\n')
	}
		
	.Object
})



setMethod('kseThreshold', 'KSEDistribution', function(object, n, p, max.kse, alpha=.05, max.k=120, n.unique.p=2000, n.x=1000, summed.weights.margin = 1e-15) {
	stopifnot( length(n) == 1 | length(p) == 1 )

	#
	# the discrete case, where we will loop over unique(n)
	#
	if (length(n) > 1 & length(p) == 1) {
		if (length(unique(n)) == 1) {
			u.n <- n[1]
			idx <- 1
		} else {
			u.n <- sort(unique(n))
			idx <- match(n, u.n)
		}
		
		# build the weight matrix. Rows are 'k', columns are 'u.n's
		w.mat <- sapply(u.n, function(n) dbinom(0:max.k, size=n, prob=p))
		if (length(u.n) == 1) {
			w.mat <- matrix(w.mat, ncol=1)
		}
	}


	#
	# the continous where we have multiple 'p' values
	#
	else if (length(n) == 1 & length(p) > 1) {

		if (length(unique(p)) == 1) {
			u.p <- p[1]
			idx <- 1
		} else if (length(unique(p)) <= n.unique.p) {		
			u.p <- sort(unique(p), decreasing=FALSE)
			idx <- findInterval(p, u.p, all.inside=TRUE)
		} else {
			u.p <- seq(from=min(p), to=max(p), length.out=n.unique.p)
			idx <- findInterval(p, u.p, all.inside=TRUE)
		}
		
		# build the weight matrix. Rows are 'k', columns are 'u.p's.
		n.bps.in.kernel <- length(object@kernelX)
		w.mat <- sapply(u.p, function(p) dbinom(0:max.k, size=n.bps.in.kernel, prob=p))
		if (length(u.p) == 1) {
			w.mat <- matrix(w.mat, ncol=1)
		}


	} else { stop('Incorrect lengths of p and n. Cannot determine if discrete of continious case for the binomial.') }




	#
	# irrespective of discrete or continious case, we now have the weight matrix
	# 
	# If columns of 'w.mat' don't sum to 1, max.k should larger!
	if (!all(apply(w.mat, 2, sum) > 1 - summed.weights.margin)) {
		stop('max.k to low')
	}
	


	calcThresholds <- function(x) {
		
		# Find the x, for which p <= alpha
		pcond.mat <- sapply(1:nrow(w.mat), function(i) {
			pkseCond(object, x, k = i-1)
		})


		# performance bottle neck
		vals <- apply(pcond.mat, 1, function(z) {
			colSums(z * w.mat, na.rm=TRUE)
		})


		thresholds <- apply(vals > (1-alpha), 1, function(z) {
			x.idx <- which(z)[1]
			if (is.na(x.idx)) {
				NA
			} else {
				x[ x.idx  ]
			}
		})
	
		if (any(is.na(thresholds))) {
			# increase max.kse by 1 every time no threshold boundary is found
			calcThresholds(x = seq(1, max(x) * 2, length.out=n.x))
		} else {
			thresholds[idx]
		}
	}

	thresholds <- calcThresholds(x = seq(1, max.kse, length.out=n.x))



	# second run, to get more accurate p-values.
	thr <- range(thresholds) * c(.9, 1.1) # add 1% margin
	
	calcThresholds(x = seq(thr[1], thr[2], length.out=n.x))
})




#
# This function calculates the cumulative probability density of x according to
# the KSE distibtion.
#
#   x: values of kse
#
# In the discrete case:
#   n: number of background site within the kernel
#   p: probability of a background site being inserted
#
# In the continious case:
#   n: number of base pairs within the kernel
#   p: the background probability
#
#   max.k:                 the maximum value of k to calulate the weights for.
#   n.unique.p:            in continous case, we bin the p-values into n.unique.p bins.
#   summed.weights.margin: the columns of the weight matrix should be 1 minus this margin.
#


setMethod('pkse', 'KSEDistribution', function(object, x, n, p, max.k=120, n.unique.p=2000, summed.weights.margin = 1e-15) {
	
	stopifnot( length(x) == length(n) | length(x) == length(p) )
	stopifnot( length(n) == 1 | length(p) == 1 )

	#
	# the discrete case, where we will loop over unique(n)
	#
	if (length(n) > 1 & length(p) == 1) {
		if (length(unique(n)) == 1) {
			u.n <- n[1]
			idx <- 1
		} else {
			u.n <- sort(unique(n))
			idx <- match(n, u.n)
		}
		
		# build the weight matrix. Rows are 'k', columns are 'u.n's
		w.mat <- sapply(u.n, function(n) dbinom(0:max.k, size=n, prob=p))
		if (length(u.n) == 1) {
			w.mat <- matrix(w.mat, ncol=1)
		}
	}


	#
	# the continous where we have multiple 'p' values
	#
	else if (length(n) == 1 & length(p) > 1) {

		if (length(unique(p)) == 1) {
			u.p <- p[1]
			idx <- 1
		} else if (length(unique(p)) <= n.unique.p) {		
			u.p <- sort(unique(p), decreasing=FALSE)
			idx <- findInterval(p, u.p, all.inside=TRUE)
		} else {
			u.p <- seq(from=min(p), to=max(p), length.out=n.unique.p)
			idx <- findInterval(p, u.p, all.inside=TRUE)
		}
		
		# build the weight matrix. Rows are 'k', columns are 'u.p's.
		n.bps.in.kernel <- length(object@kernelX)
		w.mat <- sapply(u.p, function(p) dbinom(0:max.k, size=n.bps.in.kernel, prob=p))
		if (length(u.p) == 1) {
			w.mat <- matrix(w.mat, ncol=1)
		}


	} else { stop('Incorrect lengths of p and n. Cannot determine if dicrete of continious case for the binomial.') }




	#
	# irrespective of discrete or continious case, we now have the weight matrix
	# 
	# If columns of 'w.mat' don't sum to 1, max.k should larger!
	if (!all(apply(w.mat, 2, sum) > 1 - summed.weights.margin)) {
		stop('max.k to low')
	}
	


	# Compute the pkse
	vals <- sapply(1:nrow(w.mat), function(i) {
		w.mat[i, idx] * pkseCond(object, x, k = i-1)
	})

	cond.mat <- sapply(1:nrow(w.mat), function(i) {
		pkseCond(object, x, k = i-1)
	})


	# Sum over the N
	if (length(x) > 1) {
		vals <- rowSums(vals, na.rm=TRUE)
	} else {
		vals <- sum(vals, na.rm=TRUE)
	}
	

	## FIXME
	if (range(vals)[1] < 0 | range(vals)[2] > 1 ) {
		warning("Not true: 0 <= range(vals) <= 1 !") 
	}
	vals[vals > 1] <- 1
	vals[vals < 0] <- 0

	return(vals)
})


setMethod('pkseCond', 'KSEDistribution', function(object, x, k) {
			
			
	if (k == 0) {
		# delta dirac
		as.numeric(x >= 0)
	} else if (k > object@kNormal) {
		# normal approximation
		h1 <- object@histograms[[1]]
		
		vals <- pnorm(x, mean = k * h1$mu, sd = sqrt(k * h1$var))

		return(vals)
	} else {
		# use convoluted cumulative density histrogram
		h <- object@histograms[[k]]	
		vals <- h$cumdens[
			findInterval(x, h$breaks, all.inside=TRUE, rightmost.closed=FALSE)
		]
		
		return(vals)
	}
})


setMethod('show', 'KSEDistribution', function(object) {
	cat('An object of class "KSEDistribution"\n')

	cat('scale:  ', object@scale, '\n')
	cat('D:      ', object@D, '\n')
	cat('nBins:  ', object@nBins, '\n')
	cat('kNormal:', object@kNormal, '\n')

}) 


setMethod('plot', 'KSEDistribution', function(x, y, type=c('kernel', 'pdf', 'cdf', 'qq'), ...) {
	
	type <- match.arg(type)
	
	nBins <- x@nBins
	
	if (type == 'kernel') {
		old.par <- par(mfrow=c(1,2))
		
		plot(x=x@kernelX, y=x@kernelY, type='l')
		legend('topright', paste('scale:', x@scale))
		hist(x@kernelY, breaks=nBins, plot=TRUE)
		
		par(old.par)
		
	} else {
	
		old.par <- par(mfrow=n2mfrow(length(x@histograms)))
		
		for (h in x@histograms) {
			
			if (type == 'qq') {
				
				qqplot(
					x    = h$density,
					y    = dnorm(h$mids, mean=h$mu, sd=sqrt(h$var)),
					type = 'l',
					main = paste('k=', h$N, sep=''),
					xlab = 'Density',
					ylab = 'Normal estimation'
				)
				
			} else if (type == 'pdf') {
				plot(
					x    = h$mids,
					y    = h$density,
					type = 'h',
					main = paste('k=', h$N, sep=''),
					xlab = paste('mu=', signif(h$mu,5), '\nvar=', signif(h$var,5), sep=''),
					ylab = 'Density'
				)

				lines(x=h$mids, y=dnorm(h$mids, mean=h$mu, sd=sqrt(h$var)), col='green')

				abline(v=h$mu, col='red')
				
			} else if (type == 'cdf') {
				plot(
					x    = h$mids,
					y    = h$cumdens,
					type = 'l',
					main = paste('X=', h$N, sep=''),
					xlab = paste('mu=', signif(h$mu,5), '\nvar=', signif(h$var,5), sep=''),
					ylab = 'Cumulative density'
				)
					
				lines(x=h$mids, y=pnorm(h$mids, mean=h$mu, sd=sqrt(h$var)), col='green')
			}			
			
		}
		
		par(old.par)
	}
})

