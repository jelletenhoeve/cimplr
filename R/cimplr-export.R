plot.chrscale <- function(cimplr.object, chr, scale=NULL, type=c('kse', 'bg', 'p', 'ratio')) {
  
  
  type <- match.arg(type)
  scale <- .formatScales(scale)
  
  with(c(cimplr.object$input, cimplr.object$data[[paste(chr, scale, sep='.')]], cimplr.object$output), {
    
    
    if (type == 'kse') {
    
      plot(x, inskse, type='l', xlab=chr, ylab='kse')
      
      if (exists('th'))
        lines(x, th, col='red')
    }

    if (type == 'bg') {
      plot(x, bg , type='l', xlab=chr, ylab=biasmap.type, main=biasmap)
    }
    
    
    if (type == 'p') {
      plot(x, -log10(p), type='l', xlab=chr, ylab='-log(p)')
      #  				abline(h=-log10(alpha.level / n.peaks.total), lwd=2, col='blue')
      #					abline(h=-log10(alpha.level / n.bgsites.total), lwd=2, col='red')
      #  abline(h=-log10(.05), col='red')
      #	abline(h=-log10(.01), col='red')
      #	abline(h=-log10(.001), col='red')
    }
      
    if (type == 'ratio') {
      plot(x, ratio, type='l', xlab=chr, ylab='ratio')
      abline(h= 1, lwd=2, col='red')
    }
    
  })
  
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


export.cis <- function(cimplr.object, what=c('collapsed', 'all'), type=c('bed', 'txt'), file, annotation.file, sep=',', score.column='peak.pvalue', maxgap=100e3) {
  what <- match.arg(what)
  type <- match.arg(type)
  
  cat('', sep='', file=file, append=FALSE) # create a new file
  cises <- switch(what,
    collapsed= cimplr.object$output$collapsed.cises,
    all = cimplr.object$output$all.cises
  )
  
  cis.peaks <- cises
  start(cis.peaks) <- cises$peak.location
  end(cis.peaks) <- cises$peak.location
  
  
  # do the annotation
  if(!missing(annotation.file) ) {
  
    message('annotating CIS using: ', annotation.file)
    annots <- import(annotation.file, asRangedData=FALSE) # contains 'genes' object.
    
    cises$nearest <- NA
    cises$others <- NA
    cises$maxgap <- maxgap
    ov <- findOverlaps(cises, annots, maxgap=maxgap)
    for (i in unique(queryHits(ov))) {
      
      genes <- annots[subjectHits(ov)[queryHits(ov) == i]]
      
      nearest.idx <- nearest(cis.peaks[i], genes)
      others.idx <- setdiff(1:length(genes), nearest.idx)
      
      cises[i]$nearest <- paste(elementMetadata(genes[nearest.idx])$name, collapse='|')
      if (length(others.idx) > 0) {
        cises[i]$others  <- paste(elementMetadata(genes[others.idx])$name, collapse='|')
      }
    }
  }
  
  
  if (type=='bed') {
    
    score.column <- match.arg(score.column, choices=colnames(elementMetadata(cises)))
    cises$score <- elementMetadata(cises)[[score.column]]
    
    export.bed(cises, con=file)
  } else {
  
    for (i in 1:length(cises)) {
      elementMetadata(cises[i])
      cis_header <- c(
        paste('id: ', names(cises[i]), sep=''),
        paste('chr: ', as.character(seqnames(cises[i])), sep=''),
        paste('start: ', start(cises[i]), sep=''),
        paste('end: ', end(cises[i]), sep=''),
        paste('width: ', width(cises[i]), sep='')
      )
      cis_header <- paste('# ', cis_header, sep='')
      cat(cis_header, sep='\n', file=file, append=TRUE)
      
      meta_data <- paste('# ', names(elementMetadata(cises[i])), ': ', unlist(elementMetadata(cises[i])), sep='')  
      cat(meta_data, sep='\n', file=file, append=TRUE)
      
      
      insertions <- cimplr.object$input$insertions.org
      insertions <- insertions[insertions$chr == as.character(seqnames(cises[i])), ]
      
      start <- start(cises[i])
      end <- end(cises[i])
      scale <- elementMetadata(cises[i])$scale
      insertions <- insertions[insertions$location >=  (start - 3 * scale) & insertions$location <=  (end + 3 * scale), ]
      
      insertion_header <- paste('# insertions: ', paste(colnames(insertions), collapse=sep), sep='')  
      cat(insertion_header, sep='\n', file=file, append=TRUE)
      write.table(insertions, file=file, sep=sep, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    }
  }
}

