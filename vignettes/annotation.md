CIMPLR - Generate a CIS annotation file (genes)
========================================================




Create a annation BED file from Ensemble
-------------------------



```r
library(biomaRt)
library(rtracklayer)
```

```
## Loading required package: GenomicRanges
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
## Loading required package: IRanges
```



```r
chromosomes <- c(paste("chr", 1:19, sep = ""), "chrX", "chrY")

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_id", "chromosome_name", 
    "start_position", "end_position", "strand"), filters = "chromosome_name", 
    values = substring(chromosomes, 4), mart = mart)


genes$strand <- c("-", "+")[(genes$strand + 3)/2]
genes$chromosome_name <- paste("chr", genes$chromosome_name, sep = "")

str(genes)
```

```
## 'data.frame':	30431 obs. of  6 variables:
##  $ ensembl_gene_id : chr  "ENSMUSG00000061024" "ENSMUSG00000025909" "ENSMUSG00000098034" "ENSMUSG00000051285" ...
##  $ external_gene_id: chr  "Rrs1" "Sntg1" "Serpinb10" "Pcmtd1" ...
##  $ chromosome_name : chr  "chr1" "chr1" "chr1" "chr1" ...
##  $ start_position  : int  9545408 8361475 107535508 7088920 6487231 162647203 162639150 162570515 162532127 162477680 ...
##  $ end_position    : int  9547454 9299878 107547072 7173628 6860940 162658173 162649693 162599084 162548551 162479943 ...
##  $ strand          : chr  "+" "-" "+" "+" ...
```




```r

granges <- with(genes, GRanges(ranges = IRanges(start_position, end_position), 
    strand = Rle(factor(strand)), seqnames = Rle(factor(chromosome_name)), name = external_gene_id))

granges
```

```
## GRanges with 30431 ranges and 1 metadata column:
##           seqnames                 ranges strand   |        name
##              <Rle>              <IRanges>  <Rle>   | <character>
##       [1]     chr1 [  9545408,   9547454]      +   |        Rrs1
##       [2]     chr1 [  8361475,   9299878]      -   |       Sntg1
##       [3]     chr1 [107535508, 107547072]      +   |   Serpinb10
##       [4]     chr1 [  7088920,   7173628]      +   |      Pcmtd1
##       [5]     chr1 [  6487231,   6860940]      +   |        St18
##       ...      ...                    ...    ... ...         ...
##   [30427]     chr5 [ 29481525,  29481631]      -   |     Gm24538
##   [30428]     chr5 [ 23710991,  23711061]      -   |     Mir3096
##   [30429]     chr7 [ 59362628,  59362706]      -   |     Gm24417
##   [30430]     chr9 [113776182, 113776298]      -   |     Gm25059
##   [30431]     chr9 [ 84649496,  84649659]      +   |     Gm26126
##   ---
##   seqlengths:
##     chr1 chr10 chr11 chr12 chr13 chr14 ...  chr7  chr8  chr9  chrX  chrY
##       NA    NA    NA    NA    NA    NA ...    NA    NA    NA    NA    NA
```




```r
export(granges, con = "ensembl_genes.bed")
```




