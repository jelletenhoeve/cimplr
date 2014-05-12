getArtificialIsetIds <- function(conn, artificial_model_id=1) {
	sql <- paste("	SELECT *
					FROM `artificial_iset`
					WHERE
						artificial_model_id= ", artificial_model_id, sep='')
	
	artificial_isets <- dbGetQuery(conn, sql)
	artificial_isets$artificial_iset_id
}

getArtificialCIS <- function(conn, artificial_model_id=1) {
	sql <- paste("SELECT * FROM `artificial_cis` WHERE artificial_model_id='", artificial_model_id,"'", sep='')
	
	dbGetQuery(conn, sql)

	artificial_cis <- dbGetQuery(conn, sql)
	
	GRanges(
		seqnames = Rle(artificial_cis$seqname),
		ranges = IRanges(
			start = artificial_cis$start,
			end = artificial_cis$end
		),
		artificial_cis_id = artificial_cis$artificial_cis_id,
		gene_id = artificial_cis$gene_id,
		strand     = Rle('*')
	)
}

getArtificialIset <- function(conn, artificial_iset_id) {
	
	sql <- paste("SELECT * FROM `artificial_iset` WHERE artificial_iset_id=", artificial_iset_id, sep="")
	
	db_iset <- dbGetQuery(conn, sql)
	if (nrow(db_iset) != 1) {
		stop("'artificial_iset_id'=", artificial_iset_id, " not found!")
	}
		
	data.frame(
		chr  = fromJSON(db_iset$seqname),
		location = fromJSON(db_iset$location),
		stringsAsFactors = FALSE
	)
}



writeEvaluation <- function(conn, actual, predicted, artificial_iset_id, caller) {
	

	I <- overlapsAny(predicted, actual)
	I2 <- overlapsAny(actual, predicted)
	
	tp <- sum(I2)
	fp <- sum(!I)
	fn <- sum(!I2)

	evaluation <- data.frame(
		artificial_iset_id = artificial_iset_id,
		caller = toJSON(caller),
		tp=tp,
		fp=fp,
		fn=fn
	)

	# Table 1: evaluation
	dbWriteTable(conn, name="evaluation", row.names=FALSE, append=TRUE, overwrite=FALSE, evaluation)
	evaluation_id <- dbGetQuery(conn, 'SELECT LAST_INSERT_ID()')[[1]]
	
	
	# Table 2: evaluation_predicted_cis
	if (length(predicted) != 0) {
		evaluation_predicted_cis <- data.frame(
			evaluation_id = evaluation_id,
			seqname = as.character(seqnames(predicted)),
			start = start(predicted),
			end = end(predicted),
			actual = I
		)	
		dbWriteTable(conn, name="evaluation_predicted_cis", row.names=FALSE, append=TRUE, overwrite=FALSE, evaluation_predicted_cis)
	}				

	
	# Table 3: evaluation_artificial_cis
	evaluation_artificial_cis <- data.frame(
		evaluation_id = evaluation_id,
		artificial_cis_id = elementMetadata(actual)$artificial_cis_id,
		predicted = I2
	)
	dbWriteTable(conn, name="evaluation_artificial_cis", row.names=FALSE, append=TRUE, overwrite=FALSE, evaluation_artificial_cis)

}


getEvaluations <- function(conn, callers, artificial_model_id=1) {

	sql <- paste(
		"SELECT
			caller,
			lambda,
			tp,
			fp,
			fn,
			tp / (tp + fn) as tpr,
			tp / (tp + fp) as ppv,
			fp / (tp + fp) as fdr,
			IF(fp > 0, 1, 0) as fwe
		FROM evaluation
		LEFT JOIN artificial_iset USING (artificial_iset_id)
		WHERE
			artificial_model_id = ", artificial_model_id, " AND
			caller IN (", paste("'", sapply(callers, toJSON), "'", collapse=', ', sep=''), ")", sep='')
#browser()
	dbGetQuery(conn, sql)
}



plotPerformance <- function(conn, callers, artificial_model_id=1, pi='tp', x='lambda', legend='bottomright', ...) {
	



	all.evals <- getEvaluations(conn, callers, artificial_model_id=artificial_model_id)


	xlim <- range(all.evals[[x]])
	ylim <- c(0, max(all.evals[[pi]], na.rm=TRUE))
	plot(1, type='n', xlab=x, ylab=pi, xlim=xlim, ylim=ylim, bty='n', ...)
	

	
	cols <- rainbow(length(callers))
	lwds <- rep(1, length(callers))

	for (i in 1:length(callers)) {
		evals <- all.evals[all.evals$caller == toJSON(callers[[i]]), ]
		means <- tapply(evals[[pi]], evals[[x]], mean, na.rm=TRUE)
		lines(x=as.numeric(names(means)), y=means, col=cols[i], lwd=lwds[i])
	}
	
	if (pi == 'fwe' | pi == 'fdr') {
		abline(h=0.05, col='black', lty=2)
	}
	
	if (!is.na(legend)) {
		legend(legend, names(callers), col=cols, lwd=lwds, title='Callers')
	}
}
