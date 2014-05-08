Using CIMPLR to call Common Insertion Sites using bias maps without permutation testing
========================================================


In preparation for running the actual CIMPLR commands, we start by preparing an R environment containing CIMPLR functions and insertion data.

First set your working directory to the root folder of the package. Here we do that by setting the 'Knit' options to create this Vignette; use `setwd('../')` at the console.


```r
opts_knit$set(root.dir = "../")
```


Load the necessary libraries and sources of the package.


```r

library(Biobase)
```

```
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     Filter, Find, Map, Position, Reduce, anyDuplicated,
##     as.data.frame, cbind, colnames, duplicated, eval, get,
##     intersect, lapply, mapply, match, mget, order, paste, pmax,
##     pmax.int, pmin, pmin.int, rank, rbind, rep.int, rownames,
##     sapply, setdiff, sort, table, tapply, union, unique, unlist
## 
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```r
library(GenomicFeatures)
```

```
## Loading required package: IRanges
## Loading required package: GenomicRanges
## Loading required package: AnnotationDbi
```

```r
library(rjson)
library(multicore)
```

```
## 
## Attaching package: 'multicore'
## 
## The following objects are masked from 'package:parallel':
## 
##     mclapply, mcparallel, pvec
```

```r
library(rtracklayer)
library(raster)
```

```
## Loading required package: sp
## 
## Attaching package: 'sp'
## 
## The following object is masked from 'package:IRanges':
## 
##     %over%
## 
## 
## Attaching package: 'raster'
## 
## The following object is masked from 'package:rtracklayer':
## 
##     values
## 
## The following objects are masked from 'package:AnnotationDbi':
## 
##     direction, select
## 
## The following objects are masked from 'package:GenomicRanges':
## 
##     distance, shift, values, values<-
## 
## The following objects are masked from 'package:IRanges':
## 
##     distance, shift, trim, values, values<-
```

```r

source("R/AllClasses.R")
source("R/AllGenerics.R")
```

```
## Creating a new generic function for 'metadata' in the global environment
```

```r
source("R/KSEDistribution.R")
source("R/core.R")
source("R/cimplr.R")
```


Load your data as a `data.frame` containing at least a `chr` and `location` column.


```r
load("data/test_iset.rda")

head(iset)
```

```
##     chr location
## 1 chr10 14550904
## 2 chr10 34266702
## 3 chr10 45748052
## 4 chr10 35169341
## 5 chr10 19266925
## 6 chr10 52070550
```




Now we can start with CIMPLR.


Step 1: Initialise a CIMPLR object
-------------------------

We construct a CIMPLR object based on the insertion data and the parameters for the CIMPLR run.


```r
co <- cimplr(insertions = iset, reference = "mm9", pattern = "TA", scales = log.seq(5000, 
    5e+05, 100)[seq(1, 50, 5)], pkse.method = "count", alpha.level = 0.05, p.adjust.method = "bonferroni", 
    p.adjust.n.method = "n.x", mtf = 1, test.chromosomes.seperately = TRUE, 
    n.cores = c(scales = 10, chromosomes = 1), genes.file = NULL)
```

```
## *** CIMPLR *** 
## 
## List of 22
##  $ n.bp.total                 : int 129993255
##  $ n.bgsites.total            : int 8228738
##  $ n.insertions.total         : int 1000
##  $ chr.info                   :'data.frame':	1 obs. of  3 variables:
##   ..$ length      : int 129993255
##   ..$ TA.count    : int 8228738
##   ..$ n.insertions: int 1000
##  $ insertions.org             :'data.frame':	1000 obs. of  2 variables:
##   ..$ chr     : chr [1:1000] "chr10" "chr10" "chr10" "chr10" ...
##   ..$ location: num [1:1000] 14550904 34266702 45748052 35169341 19266925 ...
##  $ chromosomes                : chr "chr10"
##  $ insertions                 :List of 1
##   ..$ chr10: num [1:1000] 14550904 34266702 45748052 35169341 19266925 ...
##  $ reference                  : chr "mm9"
##  $ pattern                    : chr "TA"
##  $ exclude.chromosomes        : NULL
##  $ scales                     : num [1:10] 5000 6309 7961 10046 12677 ...
##  $ pkse.method                : chr "count"
##  $ alpha.level                : num 0.05
##  $ p.adjust.method            : chr "bonferroni"
##  $ p.adjust.n.method          : chr "n.x"
##  $ mtf                        : num 1
##  $ test.chromosomes.seperately: logi TRUE
##  $ stepsize                   : num 1000
##  $ D                          : num 4
##  $ genes.file                 : NULL
##  $ reference.dir              : chr "./references"
##  $ n.cores                    : Named num [1:2] 10 1
##   ..- attr(*, "names")= chr [1:2] "scales" "chromosomes"
```



Step 2: Load a bias map
-------------------------

As a bias map, we load the TA count files under the './references' folder. This bias map is simply the TA site background, but more sophisticated bias maps can be supplied here (link to CreateBiasMap Vignette).


```r
co <- cimplr.reference(co)
```

```
## cimplr.reference()
```




Step 3: Gaussian Kernel Convolution
-------------------------

With the insertion data and the bias map loaded and pre-processed within the CIMPLR object, we can now start the kernel convolution, resulting in Kernel Smoothed Estimates (KSEs) for all scales.

```r
co <- cimplr.convolve(co)
```

```
## cimplr.convolve()
```




Step 4: Calculate the threshold
-------------------------

A seperate step calculates the threshold for the KSE at each scale. This is done using the analytical KSE distribution at a controlled FDR and taking into account multiple scales (link to KSEdistribution Vignette).

```r
co <- cimplr.threshold(co)
```

```
## cimplr.threshold()
## calc.fdr.pval()
```




Step 5: Call CISes
-------------------------

With the thresholds in place, the last step is to do the CIS calling. This returns two lists of CISes: (1) all CIS for each scale, and (2) the cross scale CIS using the 'cis-cloud-global-maximum' procedure (link to CrossScaleCIScalling Vignette)


```r
co <- cimplr.call(co)
```

```
## cimplr.call()
## cimplr.call() - merge CIS list
## cimplr.call() - call CIS across scales (chr = chr10)
```

```r


head(co$output$all.cises)
```

```
##                     chromosome    start      end width peak.location
## CIS10:5941001_5000       chr10  5919001  5941001 22000       5941001
## CIS10:11623001_5000      chr10 11620001 11623001  3000      11623001
## CIS10:57203001_5000      chr10 57203001 57234001 31000      57203001
## CIS10:71380001_5000      chr10 71380001 71388001  8000      71380001
## CIS10:72966001_5000      chr10 72959001 72966001  7000      72966001
## CIS10:76959001_5000      chr10 76959001 76978001 19000      76959001
##                     peak.kse.value peak.bg peak.ratio
## CIS10:5941001_5000           2.317    2262      1.149
## CIS10:11623001_5000          2.319    2739      1.041
## CIS10:57203001_5000          2.312    2905      1.008
## CIS10:71380001_5000          2.238    2722      1.009
## CIS10:72966001_5000          2.350    3032      1.004
## CIS10:76959001_5000          1.958    1625      1.003
##                     n.insertions.within.3.sigma scale alpha.level
## CIS10:5941001_5000                           13  5000        0.05
## CIS10:11623001_5000                           3  5000        0.05
## CIS10:57203001_5000                          12  5000        0.05
## CIS10:71380001_5000                           4  5000        0.05
## CIS10:72966001_5000                           4  5000        0.05
## CIS10:76959001_5000                           5  5000        0.05
##                     p.adjust.n p.adjust.method p.adjust.n.method mtf
## CIS10:5941001_5000      129994      bonferroni               n.x   1
## CIS10:11623001_5000     129994      bonferroni               n.x   1
## CIS10:57203001_5000     129994      bonferroni               n.x   1
## CIS10:71380001_5000     129994      bonferroni               n.x   1
## CIS10:72966001_5000     129994      bonferroni               n.x   1
## CIS10:76959001_5000     129994      bonferroni               n.x   1
##                     adjusted.alpha.level
## CIS10:5941001_5000             3.846e-07
## CIS10:11623001_5000            3.846e-07
## CIS10:57203001_5000            3.846e-07
## CIS10:71380001_5000            3.846e-07
## CIS10:72966001_5000            3.846e-07
## CIS10:76959001_5000            3.846e-07
```

```r
head(co$output$collapsed.cises)
```

```
##                      chromosome    start      end width peak.location
## CIS10:5941001_5000        chr10  5919001  5941001 22000       5941001
## CIS10:11623001_5000       chr10 11620001 11623001  3000      11623001
## CIS10:57243001_12676      chr10 57196001 57243001 47000      57243001
## CIS10:71379001_10046      chr10 71379001 71392001 13000      71379001
## CIS10:72957001_10046      chr10 72957001 72969001 12000      72957001
## CIS10:76498001_6309       chr10 76498001 76500001  2000      76498001
##                      peak.kse.value peak.bg peak.ratio
## CIS10:5941001_5000            2.317    2262      1.149
## CIS10:11623001_5000           2.319    2739      1.041
## CIS10:57243001_12676          3.339    7757      1.031
## CIS10:71379001_10046          2.967    5605      1.031
## CIS10:72957001_10046          3.058    6097      1.041
## CIS10:76498001_6309           2.030    2198      1.016
##                      n.insertions.within.3.sigma scale alpha.level
## CIS10:5941001_5000                            13  5000        0.05
## CIS10:11623001_5000                            3  5000        0.05
## CIS10:57243001_12676                          14 12677        0.05
## CIS10:71379001_10046                           5 10046        0.05
## CIS10:72957001_10046                           4 10046        0.05
## CIS10:76498001_6309                            3  6309        0.05
##                      p.adjust.n p.adjust.method p.adjust.n.method mtf
## CIS10:5941001_5000       129994      bonferroni               n.x   1
## CIS10:11623001_5000      129994      bonferroni               n.x   1
## CIS10:57243001_12676     129994      bonferroni               n.x   1
## CIS10:71379001_10046     129994      bonferroni               n.x   1
## CIS10:72957001_10046     129994      bonferroni               n.x   1
## CIS10:76498001_6309      129994      bonferroni               n.x   1
##                      adjusted.alpha.level
## CIS10:5941001_5000              3.846e-07
## CIS10:11623001_5000             3.846e-07
## CIS10:57243001_12676            3.846e-07
## CIS10:71379001_10046            3.846e-07
## CIS10:72957001_10046            3.846e-07
## CIS10:76498001_6309             3.846e-07
```


