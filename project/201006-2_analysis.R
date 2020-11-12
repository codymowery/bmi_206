## 201029
## Analyzing .txt files downloaded from /wynton/group/ye/cody/syntrons/data/201006-2_Jake/fastq/200818-1/ after running "zcat input1*R1* | grep CCTACCTCTGTGCTGTGAGGGTATG | uniq > input1_r1.txt"
setwd('~/github/syntrons/experiments/201006-2/')

library(tidyverse)
library(Biostrings)
library(ggseqlogo)

flist = list.files('~/Downloads/', pattern = 'r1.txt.gz', full.names = T)
overhangs = read_csv('~/github/syntrons/gene_sequences/overhangSets_HamediRadEtAl2019.csv')

data = tibble(sample = character(), read = character())
data
for(file in flist){
  infile = read_table(file, col_names = F) %>% dplyr::rename(read=X1)
  infile$sample = str_replace(file, pattern = '/Users/codymowery/Downloads//', replacement = '')
  data = data %>% bind_rows(infile)
}
write_rds(data %>% filter(grepl('pos',sample)), '~/Box/Projects/Syntrons/201006-2/201006-2_posRawReads.rds', compress = 'gz')
write_rds(data %>% filter(grepl('input',sample)), '~/Box/Projects/Syntrons/201006-2/201006-2_inputRawReads.rds', compress = 'gz')

split = data %>% 
  mutate(start = substr(read, start=1, stop=20),
         intron = substr(read, start=21, stop=105),
         end = substr(read, start=106, stop=nchar(read))) %>% 
  dplyr::select(-read)
write_rds(split %>% filter(grepl('pos',sample)), '~/Box/Projects/Syntrons/201006-2/201006-2_posSplitReads.rds', compress = 'gz')
write_rds(split %>% filter(grepl('input',sample)), '~/Box/Projects/Syntrons/201006-2/201006-2_inputSplitReads.rds', compress = 'gz')

## ggSeqLogos for intron conservation in NYESO+ vs input
split_pos = split %>% filter(grepl('pos',sample))
plot = DNAStringSet(split_pos$intron)
cons = consensusMatrix(plot, as.prob = T)
ggseqlogo(cons, seq_type='dna', method = 'bits', col_scheme='auto') + theme(axis.text.x = element_blank())

split_input = split %>% filter(grepl('input',sample))
plot = DNAStringSet(split_input$intron)
cons = consensusMatrix(plot, as.prob = T)
ggseqlogo(cons, seq_type='dna', method = 'bits', col_scheme='auto') + theme(axis.text.x = element_blank())


##########
## Look at improper splicing
##########
test = read_rds('~/Box/Projects/syntrons/201006-2/201006-2_posSplitReads_fixed.rds')
test
fivep_wrong = test %>% filter(!(substr(intron, start = 1, stop = 5)=='GTATG'))
plot = DNAStringSet(fivep_wrong$intron)
cons = consensusMatrix(plot, as.prob = T)
ggseqlogo(cons, seq_type='dna', method = 'bits', col_scheme='auto') + theme(axis.text.x = element_blank())

threep_wrong = test %>% filter(!(substr(intron, start = nchar(intron)-3, stop = nchar(intron))=='ATAG'))
