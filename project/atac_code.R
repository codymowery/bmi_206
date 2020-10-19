library(tidyverse)

infile = read_tsv('~/github/bmi_206/project/GSE99287_Retina_ATACSeq_peak_counts.txt')
infile
nrow(infile)
length(grep('Retina', colnames(infile), value = TRUE))
length(grep('NOR', grep('Retina', colnames(infile), value = TRUE)))
length(grep('AMD', grep('Retina', colnames(infile), value = TRUE)))

long = pivot_longer(infile, cols = grep('Retina', colnames(infile), value = TRUE), names_to = 'Sample', values_to = 'counts')
long
table(is.na(long$counts))

table(is.na(infile[,9:ncol]))
table(infile$Type)
table(is.na(infile))
