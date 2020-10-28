##201024
## Cody Mowery

setwd('~/github/bmi_206/lab2/')
#Read in the genotype and phenotype matrices
genos = as.matrix(read.table("./genos.txt"))
phenos = as.matrix(read.table("./phenos.txt"))

#Make a histogram of the phenotypes. Do they look normally distributed? 
hist(phenos)

#How are the genotypes encoded? 
table(genos)

#How many individuals are there in the data set and how many SNPs?
dim(genos)
dim(phenos)
N = 1500 ## There are 1500 individuals
M = 10000 ## There are 10k SNPs

#Compute the *minor* allele frequency of every SNP
MAFs = array(0,M)
## Here, I am calculating the number of non-0 alleles for each variant
## Then I divide by the total number of alleles counted per variant (N individuals * 2 alleles/individual) 
## Then I ensure I am reporting the MAF by storing the minimum value of the non-0 allele frequency, or 1 minus that value
for(i in 1:M) {
      af = sum(genos[,i])/(N*2)
      MAFs[i] = min(af, 1-af)
}
MAFs

#Run a GWAS on the data under an additive model and save the p-values, z-scores, and effect sizes
pvalues = array(0,M)
zscores = array(0,M)
betas = array(0,M)
for(i in 1:M) {
	g = genos[,i]
	res = lm(phenos ~ g)
	zscores[i] = cor(g, phenos) * sqrt(N) ## Calculate Z-scores as we did in live coding session, 
	pvalues[i] = summary(res)$coefficients[2,4] ## Extract p-value from lm() output
	betas[i] = summary(res)$coefficients[2,1] ## Extract beta coefficient from lm() output
}

#Draw a QQ plot for log10(p) values
obsLogPvs = sort(-log10(pvalues))
expLogPvs = sort(-log10(seq(1/M,1,1/M)))

plot(expLogPvs,obsLogPvs,main='QQ plot')
abline( a=0, b=1 )

#### PLOT: shows that most variants lie along expected abline but that a handful of significant variants deviate from expectations

#Compute the chi2 test statistics
chis = zscores^2

#Is there inflation?
lambdaGC = median(chis)/0.454 # why .454?
## I'm a little confused by the 0.454
## Googling (http://genometoolbox.blogspot.com/2014/08/how-to-calculate-genomic-inflation.html) gives me qchisq(0.5,1) = 0.454
## So it's just the expected median of chi2 values? and we look for deviation from that distribution? 
abline( a=0, b=lambdaGC, col=2 )
#### PLOT: new abline almost superimposable on previous abline, so no significant inflation from handful of variants

#Are there any signficantly associated SNPs? If so, which SNPs are they?
assoc = which(pvalues<(0.5/M)) 

#Build a linear predictor of the phenotype using these associated SNPs. 
ypred = array(0,N)
for(i in 1:N) {
      ypred[i] = genos[i,assoc] %*% betas[assoc] 
      ## matrix multiplication of the sig genos and their betas returns a predicted value
      ## that can be compared to the measured value in phenos
}
plot(ypred,phenos)

#### PLOT: actual and predicted phenotypes lie pretty much on the diagonal, so predictor appears to perform well

#What is the correlation between the predicted phenotype and the true phenotype 
cor(ypred,phenos) ## 0.958 is highly correlated

#Test each of the associated SNPs for non-linearity.
hp = array(0,length(assoc))
for (i in 1:length(assoc)) {
  g = genos[,assoc[i]]
  h = g
  h[h==2]=0
  #Hint: can use anova(lm(?),lm(?)) or summary(lm(?))
  hp[i] <- summary(lm(phenos~g+h))$coefficients[3,4]
}

is_lin = hp < 0.05 ## Determine which variants are non-linear by looking for p<0.05
is_lin
lin = min(which(is_lin==TRUE)) ## Get first example of non-linear variant
not_lin = min(which(is_lin==FALSE)) ## Get first example of linear variant

#Visualize a linear SNP and a non-linear SNP
par( mfrow=c(1,2) )
plot(genos[,assoc[lin]], phenos)
points( c(0,1,2), tapply(phenos, genos[,assoc[lin]], mean ), col=2, pch=16, cex=3 )
lines( c(0,1,2), tapply(phenos, genos[,assoc[lin]], mean ), col=2, lwd=2  )
plot(genos[,assoc[not_lin]], phenos)
points( c(0,1,2), tapply(phenos, genos[,assoc[not_lin]], mean ), col=2, pch=16, cex=3 )
lines( c(0,1,2), tapply(phenos, genos[,assoc[not_lin]], mean ), col=2, lwd=2  )

#### PLOT: For non-linear plot, there is a steeper dropoff from 1-2 than 0-1, whereas expression for linear plot changes the same 0-1-2

##########################################################
######### Simulating genos with LD #######################
##########################################################

### NOTE ###
## I didn't really understand what's going on here, so I went into Katie's Solutions code and went through step by step, below

N = 1000 #number of individuals
M = 30   #number of non-causal SNPs
gs = matrix(0,nrow=N,ncol=M)

MAF = .5
gC = rbinom(N,1,MAF) #causal variant

MAF = 0.5 #minor allele frequency of all SNPs
set.seed = (42) #set random seed so we all get the same numbers

#Generate 10 tight LD partners
rho = 0.9
for(i in 1:10) {
  idx = rbinom(N,1,rho)
  gs[,i]=gC*idx+rbinom(N,1,MAF)*(1-idx)
  # test they have the right LD empirically
  cat( 'Observed LD = ', cor( gs[,i], gC ), '\n' )
  # Bonus: prove they have the right LD theoretically
}

#Do the same for 10 moderate LD partners (rho=0.6) 
rho = 0.6
for(i in 11:20) {
  idx = rbinom(N,1,rho)
  gs[,i]=gC*idx+rbinom(N,1,MAF)*(1-idx)
  # test they have the right LD empirically
  cat( 'Observed LD = ', cor( gs[,i], gC ), '\n' )
  # Bonus: prove they have the right LD theoretically
}

#Do the same for 10 independent SNPs (rho=0)
rho = 0
for(i in 21:30) {
  idx = rbinom(N,1,rho)
  gs[,i]=gC*idx+rbinom(N,1,MAF)*(1-idx)
  # test they have the right LD empirically
  cat( 'Observed LD = ', cor( gs[,i], gC ), '\n' )
  # Bonus: prove they have the right LD theoretically
}


beta = 0.3
pheno = gC*beta + rnorm(N)
zsC = summary(lm(pheno~gC))$coef[2,3]
zs = sapply( 1:M, function(i) summary(lm(pheno~gs[,i]))$coef[2,3] )


#Visualize the relationship between the mean z-scores at the tag SNPs and the z-score at the causal SNP
par( mfrow=c(2,2) )
breaks = hist(c(0,zsC,zs),plot=F)$breaks
hist(zs[1:10],breaks=breaks, col=1, main='LD partners')
abline(v=zsC)
hist(zs[11:20],breaks=breaks, col=2, main='Low-LD partner SNPs')
abline(v=zsC)
hist(zs[21:30],breaks=breaks, col=3, main='Independent SNPs')
abline(v=zsC)

#### PLOT: Partners in LD are closer to causal variant t-value than unrelated SNPs

#What is the empirical LD between all pairs of SNPs including the causal SNP? 
lds = cor(cbind(gs,gC)) ## look at the correlation of calculated SNPs + causal SNP gC set at beginning

#BONUS: LD score regression
#Calculate the LD scores. There should be M+1 of them
ldscores = sapply(1:(M+1), function(m) sum(lds[m,]^2))
ldscores
#Visualize LD score regression
plot( ldscores, c( zsC, zs )^2, ylab=expression(chi^2))
#### PLOT: Demonstrates 3 distinct groups in differing LD from causal SNP?

#Use LD score regression to test for inflation and to estimate heritability
lambdaGC = median(chis)/0.454 # why .454?
summary( lm( c(zsC,zs)^2 - 1 ~ ldscores ) )$coef[2,1] * M/N
#What is the true heritability?
var(gC*beta ) / var(pheno)