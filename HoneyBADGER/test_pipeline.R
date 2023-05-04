if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
BiocManager::install("ggplot2")
BiocManager::install("GenomicFeatures")
BiocManager::install("AnnotationDbi")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

setwd("C:/Users/Ahmed Lone/Desktop/COMS/final_project/git")
devtools::install("HoneyBADGER")

library("HoneyBADGER")





data(gexp) ## tumor cells
data(ref) ## reference

require(biomaRt) ## for gene coordinates
#mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org") ## version used in manuscript
mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl') ## current version 20210603- default version hg38

mart.obj <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", dataset="hsapiens_gene_ensembl") # specify the version hg19

print(gexp[1:5,1:5])

## Object of class 'Mart':
##  Using the ENSEMBL_MART_ENSEMBL BioMart database
##  Using the hsapiens_gene_ensembl dataset

hb <- new('HoneyBADGER', name='MGH31')

hb$setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE, verbose=TRUE)

hb$plotGexpProfile() ## initial visualization

hb$setMvFit(verbose=TRUE) ## model variance

hb$setGexpDev(verbose=TRUE) ## model necessary expression deviation to identify CNVs

hb$calcGexpCnvBoundaries(init=TRUE, verbose=FALSE) ## HMM

## double check what CNVs were identified
bgf <- hb$bound.genes.final
genes <- hb$genes
regions.genes <- range(genes[unlist(bgf)])
print(regions.genes)

hb$retestIdentifiedCnvs(retestBoundGenes = TRUE, retestBoundSnps = FALSE, verbose=FALSE)

## look at final results
results <- hb$summarizeResults(geneBased=TRUE, alleleBased=FALSE)
print(head(results[,1:5]))

## visualize as heatmap 
trees <- hb$visualizeResults(geneBased=TRUE, alleleBased=FALSE, details=TRUE, margins=c(5,3))

## order cells
hc <- trees$hc
order <- hc$labels[hc$order]
## plot all chromosomes
hb$plotGexpProfile(cellOrder=order)

## plot just identified cnvs
hb$plotGexpProfile(cellOrder=order, region=hb$cnvs[['gene-based']][['amp']])

hb$plotGexpProfile(cellOrder=order, region=hb$cnvs[['gene-based']][['del']])












#Helper functions










#Local honeybadger running w/ copy 



rm(list = ls(all = TRUE))
remove.packages("HoneyBADGER")


setwd("C:/Users/Ahmed Lone/Desktop/COMS/final_project/rstudio")
devtools::install("HoneyBADGER")

library(HoneyBADGER)
#library(HoneyBADGER_local)




#simulation pre-processing
r_simulated <- simulate_deletions(r, deletion_size, 1, clonality)




hb <- new('HoneyBADGER', name='MGH31')

data(r) ## alternate allele
data(cov.sc) ## total coverage

library(TxDb.Hsapiens.UCSC.hg19.knownGene) ## in order to map SNPs to genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#hb <- new('HoneyBADGER', name='MGH31')
## Add to existing hb object
hb$setAlleleMats(r.init=r, n.sc.init=cov.sc, het.deviance.threshold=0.1, n.cores=1)

## Initializing allele matrices ... 
## Creating in-silico bulk ... 
## using 75 cells ... 
## Filtering for putative heterozygous snps ... 
## allowing for a 0.1 deviation from the expected 0.5 heterozygous allele fraction ... 
## must have coverage in at least 3 cells ... 
## 5359 heterozygous SNPs identified 
## Setting composite lesser allele count ... 
## Done setting initial allele matrices!

hb$setGeneFactors(txdb) ## map SNPs to genes

## Mapping snps to genes ... 
## >> preparing features information...		 2018-01-06 12:50:05 
## >> identifying nearest features...		 2018-01-06 12:50:05 
## >> calculating distance from peak to TSS...	 2018-01-06 12:50:05 
## >> assigning genomic annotation...		 2018-01-06 12:50:05 
## >> assigning chromosome lengths			 2018-01-06 12:50:07 
## >> done...					 2018-01-06 12:50:07 
## Done mapping snps to genes!

hb$plotAlleleProfile() ## visualize individual SNPs

hb$calcAlleleCnvBoundariesDBSCAN(init=TRUE, verbose=FALSE) ## HMM

## double check what CNVs were identified
bsf <- get('bound.snps.final', slot(hb, '.xData'))
snps <- get('snps', slot(hb, '.xData'))
regions.snp <- range(snps[unlist(bsf)])
print(regions.snp)

hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, retestBoundSnps=TRUE, verbose=FALSE)

## look at final results
results <- hb$summarizeResults(geneBased=FALSE, alleleBased=TRUE)
print(head(results[,1:5]))

## visualize as heatmap 
trees2 <- hb$visualizeResults(geneBased=FALSE, alleleBased=TRUE, details=TRUE, margins=c(5,3))

## order cells
hc2 <- trees2$hc
#order2 <- hc2$labels[hc2$order]
order2 = hc2$order
## plot all chromosomes
hb$plotAlleleProfile() ## order cells by same order as previously

## plot just identified cnvs
hb$plotAlleleProfile(region=hb$cnvs[['allele-based']][['del.loh']])

print(hb$cnvs[['allele-based']][['del.loh']])
print(results)
print(r)


## compare to new order
hb$plotAlleleProfile(cellOrder=order2) 

hb$plotGexpProfile(cellOrder=order2) 

print(r)
print(cov.sc)
print(results)


# Combined analysis 
hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=TRUE, verbose=FALSE)
results <- hb$summarizeResults(geneBased=TRUE, alleleBased=TRUE)

print(results)
