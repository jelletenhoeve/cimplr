# TODO: Add comment
# 
# Author: jthoeve
###############################################################################

# InsertionSet
setGeneric('nInsertions', def=function(object) {standardGeneric('nInsertions')})
setGeneric('nSamples', def=function(object) {standardGeneric('nSamples')})
setGeneric('metadata', def=function(object) {standardGeneric('metadata')})
setGeneric('write.bed', def=function(object, ...) {standardGeneric('write.bed')})



# KSEDistribution
setGeneric('pkse', function(object, x, n, p, ...) {standardGeneric('pkse')})
setGeneric('pkseCond', function(object, x, k) {standardGeneric('pkseCond')})
setGeneric('kseThreshold', function(object, n, p, alpha, ...) {standardGeneric('kseThreshold')})


