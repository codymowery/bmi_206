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

long
table(long$Sample)

rho = long %>% filter(Gene == 'RHO')


## Load in normalized count table
infile = read.csv('~/github/bmi_206/project/GSE99287_normed_counts.csv')
## Fig 2D
library(ggpubr)
library(ggExtra)
library(tidyverse)
library(gridExtra)
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
## calculate percent of dec fc
pct = paste0(round((nrow(fig2d %>% filter(log2(fc)<0))/nrow(fig2d))*100, 1), "%")
## Plot density with pct text above
p2 = ggplot(fig2d, aes(log2(fc)))+
  # geom_histogram(color='red', bins=1000)+
  geom_density(color='red')+
  theme_pubr()+
  coord_flip()+
  geom_vline(xintercept=0, linetype="dashed", color = "black")+
  scale_x_continuous(limits = c(-4,4), labels = c(-4,-2,0,2,4))+
  theme(axis.title.y = element_blank())+
  xlab('Density')+
  annotate("text", x = -2, y = 0.3, label = pct, size = 5)
g = grid.arrange(arrangeGrob(p, p2, ncol=2, widths=c(2,1)))
g
ggsave('~/github/bmi_206/project/fig2d_normed_gridArrange.png', g)
# p3 <- ggMarginal(p, margins = 'y', color="red", size=4)
# p3
# ggsave('~/github/bmi_206/project/fig2d_normed.pdf')

## Fig 2E
fig2e = infile %>%
  mutate(fc = AMD1_Retina_MacL/AMD1_Retina_MacR,
         avg = rowMeans(select(., grep('AMD1_Retina_Mac', colnames(infile))), na.rm = TRUE))
# avg = rowMeans(select(., grep('Retina', colnames(infile))), na.rm = TRUE))
p3 = ggplot(fig2e, aes(log2(avg), log2(fc)))+
  geom_point(color='grey')+
  geom_smooth(se = FALSE)+
  theme_pubr()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  xlab('Log2 Average Peak Signal for Each ATAC Fragment')+
  ylab('Log2(FC) of Left vs Right Eye in  AMD1')+
  scale_x_continuous(limits = c(0,17), breaks = c(5,10,15))+
  scale_y_continuous(limits = c(-4,4), labels = c(-4,-2,0,2,4))
## calculate percent of dec fc
pct = paste0(round((nrow(fig2e %>% filter(log2(fc)<0))/nrow(fig2e))*100, 1), "%")
## Plot density with pct text above
p4 = ggplot(fig2e, aes(log2(fc)))+
  # geom_histogram(color='red', bins=1000)+
  geom_density(color='red')+
  theme_pubr()+
  coord_flip()+
  geom_vline(xintercept=0, linetype="dashed", color = "black")+
  scale_x_continuous(limits = c(-4,4), labels = c(-4,-2,0,2,4))+
  theme(axis.title.y = element_blank())+
  xlab('Density')+
  annotate("text", x = -2, y = 0.3, label = pct, size = 5)
g2 = grid.arrange(arrangeGrob(p3, p4, ncol=2, widths=c(2,1)))
g2
ggsave('~/github/bmi_206/project/fig2e_normed_gridArrange.png', g2)
# p5 <- ggMarginal(p3, margins = 'y', color="red", size=4)
# p4
# ggsave('~/github/bmi_206/project/fig2e_normed.pdf')

## https://www.geeksforgeeks.org/kolmogorov-smirnov-test-in-r-programming/
library(dgof)
ks.test(fig2d$fc, fig2e$fc)
ref2d = rnorm(nrow(fig2d), mean=0, sd = sd(fig2d$fc))
t.test(ref2d, log2(fig2d$fc), paired = TRUE, alternative = "two.sided")
wilcox.test(ref2d, log2(fig2d$fc), paired = TRUE, alternative = "two.sided")
ks.test(ref2d, log2(fig2d$fc), alternative = "l")
plot(ecdf(ref2d),  
     xlim = range(c(ref2d, log2(fig2d$fc))),  
     col = "blue") 
plot(ecdf(log2(fig2d$fc)),  
     add = TRUE,  
     lty = "dashed", 
     col = "red")
ref2e = rnorm(nrow(fig2e), mean=0, sd = sd(fig2e$fc))
t.test(ref2e, log2(fig2e$fc), paired = TRUE, alternative = "two.sided")
wilcox.test(ref2e, log2(fig2e$fc), paired = TRUE, alternative = "two.sided")
ks.test(ref2e, log2(fig2e$fc), alternative = "l")
plot(ecdf(ref2e),  
     xlim = range(c(ref2e, log2(fig2e$fc))),  
     col = "blue") 
plot(ecdf(log2(fig2e$fc)),  
     add = TRUE,  
     lty = "dashed", 
     col = "red")

wilcox.test(log2(fig2d$fc), log2(fig2e$fc), paired = TRUE, alternative = "two.sided")
wilcox.test(ref2e, log2(fig2e$fc), paired = TRUE, alternative = "two.sided")
wilcox.test(ref2d, log2(fig2d$fc), paired = TRUE, alternative = "two.sided")