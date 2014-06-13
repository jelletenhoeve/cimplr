calc.biasmap <- function(chr, scale, input) {

  with(input, {
  
    genome <- strsplit(biasmap, '-')[[1]][2]
    
    bw.file <- paste(biasmap, '/', biasmap.type, '-scale', as.integer(round(scale)), '.bw', sep='')
    bgranges <- import.bw(con=bw.file, asRangedData = biasmap.type == 'density', selection=GenomicSelection(genome=genome, chrom=chr, colnames='score'))
    
    x  <- start(bgranges)
    bg <- score(bgranges)
    
    list(
      chr = chr,
      scale = scale,
      locs  = insertions[[chr]],
      x          = x,
      bg         = bg,
      n.insertions       = chr.info[chr, 3],
      n.insertions.total = n.insertions.total,
      n.bgsites          = chr.info[chr, 2],
      n.bgsites.total    = n.bgsites.total,
      n.bp               = chr.info[chr, 1],
      n.bp.total         = n.bp.total,
      n.x                = length(x)
    )
  }) # end with
}


#
# Does the kernel convolution
#
calc.kse <- function(chr.object, kse_dist, pkse.method, n.insertions.total, n.bgsites.total, alpha, p.adjust.method) {
  stopifnot(chr.object$scale == kse_dist@scale)
  with(chr.object, {
    
    # insertion density
    ftn     <- list(from=min(x), to=max(x), n=length(x))
    insdens <- density(locs, bw=scale, from=ftn$from, to=ftn$to, n=ftn$n)
    inskse  <- insdens$y * length(locs) / dnorm(0, 0, sd=scale)
    
    # n.peaks
    peak.idx <- whichLocalMaxima(inskse)
    peaks    <- x[peak.idx]
    n.peaks  <- length(peaks)
    
    # p.values
    p  <- rep(NA, length(x))
    idx <- bg > 0  # optimisation 1
    
    
    p[idx] <- 1 - pkse(kse_dist, x=inskse[idx], n=bg[idx], p=n.insertions.total / n.bgsites.total)
    
    
    ## TODO: optimize chr.object
    
    c(
      chr.object,
      list(
        inskse             = inskse,
        p                  = p,
        peaks              = peaks,
        n.peaks            = n.peaks,
        pkse.method        = pkse.method,
        p.adjust.method    = p.adjust.method,
        alpha              = alpha
      )
    )
  })
}









calc.corrected.alpha <- function(cimplr.object) {
  
  with(cimplr.object$input, {
    
    if (p.adjust.method == 'none') {
      return(alpha)
    }
    
    message('calc.corrected.alpha() - p.adjust.method = ', p.adjust.method)
    # p.adjust.method=c("fdr", "BY", "bonferroni"),
    
    # total number of p-value generated
    m <- sum(sapply(cimplr.object$data, function(obj) length(obj$p)))
    
    if (p.adjust.method == 'bonferroni') {
      corrected.alpha <- alpha / m
    } else {
      # contains all pvalues across all scales and all chromosomes
      all.pvals <- unlist(lapply(cimplr.object$data, '[[', 'p'))
      p.sorted <- sort(all.pvals, na.last=NA) # note NA's are removed. p is NA when the background is < 0
      
      #m <- length(p.sorted)
      k <- 1:m
      cm <- sum(1/(1:m))
      
      th <- switch(p.adjust.method,
                   fdr = k * alpha / m,
                   BY  = k * alpha / (m * cm)
      )
      
      corrected.alpha.idx <- which(!p.sorted < th)[1]
      corrected.alpha <- p.sorted[corrected.alpha.idx]
    }
    corrected.alpha
  })
}


# to be used for alpha's of corrected alpha's or
calc.ratios <- function(chr.object, kse_dist, corrected.alpha, update.threshold=TRUE) {
  with(chr.object, {
    
    if (update.threshold) {
      chr.object$th <- rep(NA, length(x))
      idx <- bg > 0
      chr.object$th[idx] <- kseThreshold(kse_dist, n=bg[idx], p=n.insertions.total / n.bgsites.total, max.kse=ceiling(max(inskse[idx])) * 2, alpha=corrected.alpha, max.k=120, n.unique.p=2000, n.x=1000, summed.weights.margin = 1e-15)
    }
    
    
    chr.object$ratio <- chr.object$inskse / chr.object$th
    chr.object$significant <- chr.object$inskse > chr.object$th
    chr.object$corrected.alpha <- corrected.alpha
    
    chr.object
  }) # end 'with'
}



#
# Calculates CIS for each scale
#
calc.cis <- function(chr.object) {
  chr.object$cises <- with(chr.object, {
    
    
    cat (chr, '.', scale, '\n')
    # get CIS as segments
    
    significant[is.na(significant)] <- FALSE
    segs <- segments(significant)
    
    
    
    if (segs$n == 0) {
      cises <- GRanges()
    } else {
      
      # peak is defined as the 
      peak.idx <- mapply(FUN=function(start, end) {
        which.min(ratio[start:end]) + start - 1
      }, segs$start.pos, segs$end.pos)
      
      
      peak.location   <- x[peak.idx]
      peak.kse.value  <- inskse[peak.idx]
      peak.bg         <- bg[peak.idx]
      peak.ratio      <- ratio[peak.idx]
      peak.pvalue     <- p[peak.idx]
      
      
      start <- x[segs$start.pos]
      #end   <- pmin(x[segs$end.pos], sl)
      end <- pmin(x[segs$end.pos])
      
      
      n.insertions.within.3.sigma <- mapply(FUN=function(start, end) {
        sum( locs >=  (start - 3 * scale) & locs <=  (end + 3 * scale) )
      }, start, end)
      
      
      
      
      cises <- GRanges(
        ranges    = IRanges(start, end),
        seqnames  = Rle(chr),
        peak.location               = peak.location,
        peak.kse.value              = peak.kse.value,
        peak.bg                     = peak.bg,
        peak.ratio                  = peak.ratio,
        peak.pvalue                 = peak.pvalue,
        n.insertions.within.3.sigma = n.insertions.within.3.sigma,
        scale                       = scale,
        alpha                       = alpha,
        p.adjust.method             = p.adjust.method,
        corrected.alpha             = corrected.alpha
      )
      names(cises) <- paste('CIS', substring(chr, 4), ':', peak.location, '_', .formatScales(scale), sep='')
    }
    
    cises
  }) # end 'with'
  
  chr.object
  
}





#
# Does the cross-scale CIS calling given a chromosome.
#
collapse.cis <- function(chr.objects) {
  
  library(BiocGenerics)
  library(GenomicRanges)
  
  # Calculate the min.p.value CIS
  x <- chr.objects[[1]]$x # can be optimised
  n.scales <- length(chr.objects)
  ratios <- sapply(chr.objects, function(obj) obj$ratio)
  significant.mat <- sapply(chr.objects, function(obj) obj$significant)
  significant <- apply(significant.mat, 1, any, na.rm=TRUE)
  significant[is.na(significant)] <- FALSE
  
  
  
  # get CIS
  segs <- segments(significant)
  
  
  if (segs$n == 0) {
    cises <- GRanges()
  } else {
    
    
    cises <- unlist(do.call('GRangesList', mapply(SIMPLIFY=FALSE, FUN=function(seg.start, seg.end) {
      
      
      ## Convert it to a raster object
      sub.sign <- significant.mat[seg.start:seg.end, n.scales:1, drop=FALSE]
      r_significant <- raster(t(sub.sign))
      extent(r_significant) <- extent(c(seg.start-1, seg.end, 0, n.scales) + .5)
      
      sub.ratios <- ratios[seg.start:seg.end, n.scales:1, drop=FALSE]
      sub.ratios[sub.ratios < 1] <- NA
      r <- raster(t(sub.ratios))
      
      #extent(r) <- extent(c(0, ncol(sub.p.adjusted), 0, nrow(sub.p.adjusted)) + 0.5)
      extent(r) <- extent(c(seg.start-1, seg.end, 0, n.scales) + .5)
      
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
      cises <- unlist(do.call('GRangesList', lapply(1:nrow(maxXY), function(i) {
        
        peak.location <- x[maxXY[i, 'x']]
        cises_at_selected_scale <- chr.objects[[maxXY[i, 'y']]]$cises
        
        idx <- start(cises_at_selected_scale) <= peak.location & end(cises_at_selected_scale) >= peak.location
        if (sum(idx) != 1) {
          plot(r)
          plot(r2)
          plot(r_significant)
          dev.off()
          browser()
        }
        
        stopifnot( sum(idx) == 1 )
        
        
        
        cises_at_selected_scale[idx]
      })), use.names=FALSE)
            
      # select the CIS at the global max of the cloud.
      cises[which.max(cises$peak.ratio)]
      
      
    }, segs$start.pos, segs$end.pos)), use.names=FALSE)
    
  }
  cises
}



