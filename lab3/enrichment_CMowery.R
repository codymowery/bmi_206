#load packages
library(regioneR)
library(BSgenome.Hsapiens.UCSC.hg19)
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
## Here I am using the overlapRegions() function to look for shared and/or overlapping regions in the two region sets, hits1 and all
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
hitsxhits = overlapRegions(hits1, hits2)
overlapGraphicalSummary(hits1, hits2, regions.labels=c('hits1', 'hits2'))
#Q3. How many regions overlap? 
## 6 regions overlap
#    How many regions are exactly identical? 
## All 6 overlapping regions are exactly identical


#Q4. Generate a set of random genomic regions of the same size as hits1.
#     Match these to the mean and sd of the genomic length of the hits1 regions.
## Here I am using the createRandomRegions() function to pull random regions from the hg19 genome
## As coded, these random region set will be the same size as hits1 with the same (or similar) mean and SD
## Pulling random regions from hg19 with characteristics similar to hits1 will test whether hits1 x hits2 overlap is simply random and expected by chance
random = createRandomRegions(nregions = length(hits1), 
                             length.mean = mean(width(hits1)),
                             length.sd = sd(width(hits1)),
                             genome = 'hg19')

#    Do the random genomic regions overlap hits2 more or less than hits1 does?
overlapRegions(random, hits2)
## None of the randomly generated regions overlap with hits2, so less than with hits1

#    How much do the random genomic regions overlap all tested regions?
overlapRegions(random, all) 
## None of the randomly generated regions overlap with all, so less than with hits1

## My thought is that randomly selecting regions throughout the genome is a sufficiently large search space that we don't get any overlap
## But we know that hits1 and hits2 are subsets of all, so it makes more sense to randomly sample from all and test that random subset against hits2 overlap


#    Repeatedly generate genomic regions to compute a z-score for the overlap of hits1 with hits2
#     Use the set of overlaps with random regions to test the null hypothesis
#     that hits2 overlaps hits1 more than expected compared to totally random regions
overlapPermTest(A=hits1, B=hits2, ntimes=100)
#### OUTPUTS ####
# P-value: 0.0099009900990099
# Z-score: 24.8865
# Alternative: greater
# Number of iterations: 100
# Alternative: greater
# Evaluation of the original region set: 6
# Evaluation function: numOverlaps
# Randomization function: randomizeRegions
overlapPermTest(A=random, B=hits2, ntimes=100)
# P-value: 0.0693069306930693
# Z-score: 3.9383
# Number of iterations: 100
# Alternative: greater
# Evaluation of the original region set: 1
# Evaluation function: numOverlaps
# Randomization function: randomizeRegions

## Based on the larger z-score for the first perm test, it appears that the hits1 x hits2 overlap is more frequent than is expected from a random hg19 sample
## Further, there is insignificant (p=0.0693069306930693) overlap between random and hits2

#    What is the smallest p-value you could have gotten? 
#    How do the results change with number of resamples? Random seed?
overlapPermTest(A=hits1, B=hits2, ntimes=500)
# P-value: 0.00199600798403194
# Z-score: 21.335
overlapPermTest(A=hits1, B=hits2, ntimes=1000)
# P-value: 0.000999000999000999
# Z-score: 21.0876
set.seed(30)
overlapPermTest(A=hits1, B=hits2, ntimes=1000)
# P-value: 0.000999000999000999
# Z-score: 20.224

## The z-score appears relatively stable, but the p-value expectedly gets smaller with a larger number of permutations

set.seed(10) ## setting seed back to 10 as beginning


#Q5. Repeat Q4 switching the roles of hits1 and hits2. Are conclusions similar? 
#     Match these to the mean and sd of the genomic length of the hits2 regions.
random = createRandomRegions(nregions = length(hits2), 
                             length.mean = mean(width(hits2)),
                             length.sd = sd(width(hits2)),
                             genome = 'hg19')

#    Do the random genomic regions overlap hits2 more or less than hits1 does?
overlapRegions(random, hits1)
## None of the randomly generated regions overlap with hits1

#    How much do the random genomic regions overlap all tested regions?
overlapRegions(random, all) 
## None of the randomly generated regions overlap with all

overlapPermTest(A=random, B=hits1, ntimes=100)
# P-value: 0.910891089108911
# Z-score: -0.3129
# Number of iterations: 100
# Alternative: less
# Evaluation of the original region set: 0
# Evaluation function: numOverlaps
# Randomization function: randomizeRegions

## This random sample with hits2 characteristics is even less significantly overlapping with hits1 than was the previous random sample
## It's z-score is closer to zero and the P-value is very close to 1



#Q6. Create a random bootstrap sample of regions from all tested regions. 
random = resampleRegions(A = hits1, universe = all)
length(random) ## 140, matches hits1
#    Do these random regions overlap hits2 more or less than hits1 does?
randomxhits2 = overlapRegions(random, hits2)
overlapGraphicalSummary(random, hits2, regions.labels=c('random','hits2'))
dim(randomxhits2)
table(randomxhits2$type)
## This random sample of all overlaps with hits2 more frequently than does hits1
## There are 28 regions that are exactly the same, and 1 region with partial overlap (AleftB)
permTest(A=hits1, ntimes=100, randomize.function=resampleRegions, universe=all,
         evaluate.function=numOverlaps, B=hits2, verbose=FALSE)
# P-value: 0.0099009900990099
# Z-score: -4.4091
# Number of iterations: 100
# Alternative: less
# Evaluation of the original region set: 6
# Evaluation function: numOverlaps
# Randomization function: resampleRegions

#    How does this test differe from the one in Q4? Look at the z-score and p-value.
## It's a relatively close-to-zero and negative Z-score, though the P-value is still significant 

## When I run a permutation test with a random sample of all matching hits1 characteristics and tested against hits2,
##  we get a negative z-score and significant p-value


#Q7. Repeat Q6 switching the role of hits1 and hits2. Are conclusions similar?
random = resampleRegions(A = hits2, universe = all)
length(random) ## 322, matches hits2
#    Do these random regions overlap hits2 more or less than hits1 does?
randomxhits1 = overlapRegions(random, hits1)
overlapGraphicalSummary(random, hits1, regions.labels=c('random','hits1'))
dim(randomxhits2)
table(randomxhits2$type)
## This random sample of all overlaps with hits1 more frequently than does hits2
## There are 20 regions that are exactly the same ('equal')
permTest(A=hits2, ntimes=100, randomize.function=resampleRegions, universe=all,
         evaluate.function=numOverlaps, B=hits1, verbose=FALSE)
# P-value: 0.0099009900990099
# Z-score: -4.1402
# Number of iterations: 100
# Alternative: less
# Evaluation of the original region set: 6
# Evaluation function: numOverlaps
# Randomization function: resampleRegions

## Again, the p-value is significant but the z-score is negative and closer to zero

#Q8. Which null distribution would you use in your own research and why?  
## I think the random sample of "all" is a better null distribution than is the random sample of "hg19"
## because the first questions above reveal that hits1 and hits2 are subsets of 'all', 
##  it makes sense to compare their overlap to that of a random sample of 'all' with the same characteristics (size, mean, sd) as one of the groups
## This allows us to more confidently conclude whether or not to reject the null hypothesis that hits1 and hits2 overlap the same amount as expected by chance


#The next few questions involve downloading genomics data
# You can choose sets of regions, e.g, gene annotation, ChIPseq, RNAseq, ATACseq, GWAS SNPs
cpg <- toGRanges("http://www.haowulab.org/software/makeCGI/model-based-cpg-islands-hg19.txt")

#Q9. Using data you download, can you infer what function was tested in the assay that discovered
#     hits1 and hits2? Choose data sets that will be informative about candidate functions.
#     Compute overlaps or mean values of the downloaded data for the union of hits1 and hits2
#     Guess what type of genomic element these hits are (i.e., what assay was performed)
cpgxhits1 = overlapRegions(cpg, hits1)
table(cpgxhits1$type)
cpgxhits2 = overlapRegions(cpg, hits2)
table(cpgxhits2$type)
cpgxall = overlapRegions(cpg, all)
table(cpgxall$type)
## There is quite a bit of overlap between this cpg islands dataset (pulled from regioneR documentation) and the all, hits1, hits2 datasets
## Because CpG islands are enriched in gene promoters, and there is a bit of overlap as mentioned above, 
## and this overlap is predominantly "AinB" meaning the CpG region(s) are within the larger all/hits1/hits2 regions (like CpG islands are enriched within promoters)
##    I believe that "all" is a subset of predicted or experimentally measured gene promoters

#BONUS Q10. Do you think hits1 and hits2 function in the same cell type? 
#     Build on your analysis in Q9 by separately testing overlaps with hits1 and hits2
#     Choose datasets that are from several different cell types


#BONUS Q11: Try matching the random regions more closely to regions in hits1
#          On what variables will you match them?
#            e.g., spacing, chromosome, GC-content, distance to nearest gene, masking
#          How does matching affect the z-score and p-value?
