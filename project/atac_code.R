library(tidyverse)

infile = read_tsv('~/github/bmi_206/project/GSE99287_Retina_ATACSeq_peak_counts.txt')
infile
nrow(infile)
length(grep('Retina', colnames(infile), value = TRUE))
length(grep('NOR', grep('Retina', colnames(infile), value = TRUE)))
length(grep('AMD', grep('Retina', colnames(infile), value = TRUE)))

long = pivot_longer(infile, cols = grep('Retina', colnames(infile), value = TRUE), names_to = 'Sample', values_to = 'counts')
long
save_rds(long, 'GSE99287_Retina_ATACSeq_peak_counts.txt')
table(is.na(long$counts))

table(is.na(infile[,9:ncol]))
table(infile$Type)
table(is.na(infile))

long
table(long$Sample)

rho = long %>% filter(Gene == 'RHO')

## Fig 2D
library(ggpubr)
library(ggExtra)
fig2d = infile %>% 
  mutate(fc = AMD2_Retina_MacL/AMD2_Retina_MacR,
         avg = rowMeans(select(., grep('AMD2_Retina_Mac', colnames(infile))), na.rm = TRUE))
         # avg = rowMeans(select(., grep('Retina', colnames(infile))), na.rm = TRUE))
p = ggplot(fig2d, aes(log2(avg), log2(fc)))+
  geom_point(color='grey')+
  geom_smooth(se = FALSE)+
  theme_pubr()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  xlab('Log2 Average Peak Signal for Each ATAC Fragment')+
  ylab('Log2(FC) of Left vs Right Eye in  AMD2')+
  scale_x_continuous(limits = c(0,17), breaks = c(5,10,15))+
  scale_y_continuous(limits = c(-4,4), labels = c(-4,-2,0,2,4))
p3 <- ggMarginal(p, margins = 'y', color="red", size=4)
p3
ggsave('~/github/bmi_206/project/fig2d.pdf')
## Fig 2E
fig2e = infile %>% 
  mutate(fc = AMD1_Retina_MacL/AMD1_Retina_MacR,
         avg = rowMeans(select(., grep('AMD1_Retina_Mac', colnames(infile))), na.rm = TRUE))
# avg = rowMeans(select(., grep('Retina', colnames(infile))), na.rm = TRUE))
p2 = ggplot(fig2e, aes(log2(avg), log2(fc)))+
  geom_point(color='grey')+
  geom_smooth(se = FALSE)+
  theme_pubr()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  xlab('Log2 Average Peak Signal for Each ATAC Fragment')+
  ylab('Log2(FC) of Left vs Right Eye in  AMD1')+
  scale_x_continuous(limits = c(0,17), breaks = c(5,10,15))+
  scale_y_continuous(limits = c(-4,4), labels = c(-4,-2,0,2,4))
p4 <- ggMarginal(p2, margins = 'y', color="red", size=4)
p4
ggsave('~/github/bmi_206/project/fig2e.pdf')
