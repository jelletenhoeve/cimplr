# TODO: Add comment
# 
# Author: jthoeve
###############################################################################


#setClass(
#	Class = 'InsertionSet',
#	representation = representation(
#		seqname    = 'character',
#		location   = 'integer',
#		strand     = 'character',
#		metadata   = 'data.frame',
#		phenoData  = 'data.frame',
#		phenoIndex = 'integer'
#	)
#)

setClass(
	Class = 'KSEDistribution',
	representation = representation(
		scale   = 'numeric',
		D       = 'numeric',
		nBins   = 'numeric',
		kNormal = 'numeric',
		
		# derived
		kernelX     = 'numeric',
		kernelY     = 'numeric',
		histograms  = 'list'
	)
)
