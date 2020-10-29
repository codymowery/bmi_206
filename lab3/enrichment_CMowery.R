#load packages
library(regioneR)
set.seed(10)
setwd('~/github/bmi_206/lab3/')
#read in BED formatted region files: all tested regions and two sets of positives
# these are in hg19 human genome assembly coordinates
all=toGRanges(read.table("all.bed",sep="\t"))
hits1=toGRanges(read.table("hits1.bed",sep="\t"))
hits2=toGRanges(read.table("hits2.bed",sep="\t"))

#Q1. How many regions are in hits1? How many in hits2? 
length(hits1) ## 140
length(hits2) ## 322

#Q2. Are hits1 and hits2 proper subsets of all the tested regions?
#    Check how many of each set overlaps a region in all.
hits1xall = overlapRegions(hits1, all)
overlapGraphicalSummary(hits1, all, regions.labels=c('hits1', 'all'))
dim(hits1xall)
table(hits1xall$type)
## Looks like every hits1 region is in all (equal=140), but there's an additional region in all that's inside a hits1 region (BinA)
hits2xall = overlapRegions(hits2, all)
overlapGraphicalSummary(hits2, all, regions.labels=c('hits2', 'all'))
dim(hits2xall)
table(hits2xall$type)
## Similarly, every hits2 region is in all (equal=322) 
## but there are additional all regions in hits2 regions (i.e. AinB, + vice versa) and partially overlapping regions (i.e. ArightB)

#The next few questions explore the overlap of genomic regions in hits1 and hits2.
hitsxhits = suppressWarnings(overlapRegions(hits1, hits2))
overlapGraphicalSummary(hits1, hits2, regions.labels=c('hits1', 'hits2'))
#Q3. How many regions overlap? 
## 6 regions overlap
#    How many regions are exactly identical? 
## All 6 overlapping regions are exactly identical


#Q4. Generate a set of random genomic regions of the same size as hits1.
#     Match these to the mean and sd of the genomic length of the hits1 regions.
random = createRandomRegions(nregions = length(hits1), 
                             length.mean = mean(width(hits1)),
                             length.sd = sd(width(hits1)),
                             genome = 'hg19')
#    Do the random genomic regions overlap hits2 more or less than hits1 does?
overlapRegions(random, hits2)
## None of the randomly generated regions overlap with hits2
#    How much do the random genomic regions overlap all tested regions?
overlapRegions(random, all) 
## None of the randomly generated regions overlap with all
#    Repeatedly generate genomic regions to compute a z-score for the overlap of hits1 with hits2
#     Use the set of overlaps with random regions to test the null hypothesis
#     that hits2 overlaps hits1 more than expected compared to totally random regions
randomize = function(x){createRandomRegions(genome = all)}
pt1 = suppressWarnings(permTest(A=hits1, B=hits2, ntimes=100, randomize.function = randomize, evaluate.function = numOverlaps, genome = all,
                                force.parallel = FALSE))
#    What is the smallest p-value you could have gotten? 
#    How do the results change with number of resamples? Random seed?


#Q5. Repeat Q4 switching the roles of hits1 and hits2. Are conclusions similar? 


#Q6. Create a random bootstrap sample of regions from all tested regions. 
#    Do these random regions overlap hits2 more or less than hits1 does?
#    How does this test differe from the one in Q4? Look at the z-score and p-value.


#Q7. Repeat Q6 switching the role of hits1 and hits2. Are conclusions similar?


#Q8. Which null distribution would you use in your own research and why?  


#The next few questions involve downloading genomics data
# You can choose sets of regions, e.g, gene annotation, ChIPseq, RNAseq, ATACseq, GWAS SNPs

#Q9. Using data you download, can you infer what function was tested in the assay that discovered
#     hits1 and hits2? Choose data sets that will be informative about candidate functions.
#     Compute overlaps or mean values of the downloaded data for the union of hits1 and hits2
#     Guess what type of genomic element these hits are (i.e., what assay was performed)


#BONUS Q10. Do you think hits1 and hits2 function in the same cell type? 
#     Build on your analysis in Q9 by separately testing overlaps with hits1 and hits2
#     Choose datasets that are from several different cell types


#BONUS Q11: Try matching the random regions more closely to regions in hits1
#          On what variables will you match them?
#            e.g., spacing, chromosome, GC-content, distance to nearest gene, masking
#          How does matching affect the z-score and p-value?