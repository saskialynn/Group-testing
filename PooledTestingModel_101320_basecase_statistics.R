setwd("/Users/User/Desktop")
library(tidyverse)
library(ggplot2)
library(fitdistrplus)
library(DescTools)

# Incorporate decrement to positive percent agreement expected based on 95% LoD of assay and probit regression curve shape

# Get all tests
tests <- data.table::fread("alltests_1mar24jun_v1.csv")

# Keep only valid positive ct values for first tests
tests %>% filter(result == "positive", firsttest==TRUE,!is.na(cttarget)) %>% pull(cttarget) -> cts

# Fit distribution
# Plot distribution and skew/kurtosis of ct values of real first tests
plotdist(cts, histo = TRUE, demp = TRUE, breaks=18)
descdist(cts, boot = 1000)

# Fit distributions
fw <- fitdist(cts, "weibull")
fno <- fitdist(cts, "norm")
fg <- fitdist(cts, "gamma")
fln <- fitdist(cts, "lnorm")

# Plot fit distributions and generate goodness-of-fit statistics
par(mfrow = c(2, 2))
plot.legend <- c("weibull","normal","gamma", "lnorm")
denscomp(list(fw,fno,fg,fln), legendtext = plot.legend,xlim=c(0,50),breaks=25)
qqcomp(list(fw,fno,fg,fln), legendtext = plot.legend)
cdfcomp(list(fw,fno,fg,fln), legendtext = plot.legend)
ppcomp(list(fw,fno,fg,fln), legendtext = plot.legend)
gofstat(list(fw,fno,fg,fln))

# Weibull has best fit. Obtain fit parameters
summary(fw)

# Add variable of proportion Ct value >LoD to model, as surrogate for differences in viral load distribution in different populations

# Number of replicates to generate ecdf offset values for above.
set.seed(42)
x<-50000

# Determine ct value quantiles for real first tests vs. fitted Weibull distribution
quantile(cts,c(0.5,0.25,0.75,0.95,0.99))
quantile(fw,c(0.5,0.25,0.75,0.95,0.99))

# Create simulated vector of ct values with same shape as that from Weibull distribution fit to ct values from real first tests
ct_fake_input <- rweibull(x,fw[[1]][1],fw[[1]][2])

# Create matrix of desired input parameters
# Change above.lod to % samples with ct value >LoD to reflect actual population of interest. Changing "lod" itself has no effect on model output.
lod <- 35
above.lod <- seq(0.05,0.3,0.05)
sequence <- seq(-10,15,0.01)
mat<-matrix(ncol=3,nrow=length(sequence)*length(lod))

for(v in 1:length(sequence)){
  fn<-ecdf(subset(pmax(5,ct_fake_input+sequence[v]),(pmax(5,ct_fake_input+sequence[v]))<45))
		for(i in 1:length(lod)){
			a<-1-fn(lod[i])
			mat[(v-1)*length(lod)+i,1]<-lod[i]
			mat[(v-1)*length(lod)+i,2]<-a
			mat[(v-1)*length(lod)+i,3]<-sequence[v]
			}
}

# Shift ct values for each lod and %above lod
mat2<-matrix(ncol=3,nrow=length(lod)*length(above.lod))
for(i in 1:length(lod)){
for(j in 1:length(above.lod)){
	tmp<-subset(mat,mat[,1]==lod[i])
	u<-tmp[which.min(abs(above.lod[j]-tmp[,2])),3]
	mat2[(i-1)*length(above.lod)+j,1]<-lod[i]
	mat2[(i-1)*length(above.lod)+j,2]<-above.lod[j]
	mat2[(i-1)*length(above.lod)+j,3]<-u
	}
}

# probit data input
probit_input <- read.csv("probit_zscores_cts_tissue_agnostic.csv")
probit_t<-subset(probit_input,probit_input$z_value<=1.96 & probit_input$z_value>=-1.96)
probit<-probit_t[,c(2,4,3)]
probit<-probit[order(probit$z_value,probit$ct_value),]
z_scores<-as.numeric(unlist(distinct(probit,z_value)))

# Number of replicates for model
set.seed(42)
n <- 10000
pool.max<-20
probit.mode<-c("base","dsa.lower","dsa.upper","psa") 
probit.mode.index<-1 # 1 = no variation, 2, = LLN, 3 = ULN, 4 = probabilistic
probit.z.indices<-c(488,1,length(z_scores)) # 488 is a z score of 0 (base case) in the z index vector
dilution.vary.index<-1 # 1 = no variation, 2 = probabilistic

# Model loop, vary % above each LOD

for (llod in 1:length(lod)){

allfirst.poolct <- data.table::rbindlist(lapply(1:length(above.lod),function(above) {

		ct_set<-pmax(5,subset(ct_fake_input+mat2[(llod-1)*length(above.lod)+above,3],ct_fake_input+mat2[(llod-1)*length(above.lod)+above,3]<45))
		threshold.ct <- c(sample(ct_set,n,replace=T))
	
		data.table::rbindlist(lapply(1:pool.max, function (p) {
		  data.table::rbindlist(lapply(c(0.001,0.01,0.03,0.05,0.1,0.15), function(prevalence) {
		    data.table::rbindlist(lapply(0:p, function(positives) {

		      if (positives == 0) data.frame(limit=lod[llod], pool=p, pos=positives, prevalence=prevalence, above.llod = above.lod[above], 
		              concentration=0, probability = dbinom(positives, p, prevalence), random=0, z.index=0, call.each.conc=FALSE, tests=1, tn=1,tp=0,fn=0,fp=0)
		      
		      else {
				dat <- matrix(sample(threshold.ct, positives * n, replace=T), nrow=positives)
		      
				each.conc = -log2(colSums(2^-dat)/p)+ifelse(dilution.vary.index==1,0,rnorm(mean=0,sd=1.1,n=ncol(dat))) # sd of 1.1 reflects confidence interval for deviation from perfect log2 dilution in assays
			
				data.frame(
					limit=lod[llod],
					pool=p,
					pos=positives,
					prevalence=prevalence,
					above.llod=above.lod[above],
					concentration=each.conc,
					probability = dbinom(positives, p, prevalence)/n,
					random=sample(n,n)/n,
					z.index=ifelse(probit.mode.index<4,probit.z.indices[probit.mode.index],sample(1:length(z_scores),n,replace=T))) %>%
					mutate(
						call.each.conc=probit[1+(z.index-1)*571+each.conc*10-(lod[llod]-35.9)*10,2]>random,
						tests=1 + p * (call.each.conc),
						tn=0,
						tp=1 * (call.each.conc),
						fn=1 * (!call.each.conc),
						fp=0)
        }
      }))
    }))
  }))
}))

group_by(allfirst.poolct, pool, prevalence, above.llod, limit) %>% 
  summarize(pos=weighted.mean(pos, w=probability),
            total.tests=weighted.mean(tests, w=probability), 
            tests.per.sample=weighted.mean(tests, w=probability)/mean(pool),
            tn=weighted.mean(tn, w=probability), 
            tp=weighted.mean(tp, w=probability), 
            fn=weighted.mean(fn, w=probability), 
            fp=weighted.mean(fp, w=probability)) -> apw

rm(allfirst.poolct)

tmp <- as.data.frame(apw)
ifelse(llod==1,all<-apw,all<-rbind(all,tmp))
rm(tmp)
rm(apw)
}

all2 <- all %>%
  mutate(ppa=tp/(tp+fn))
cpi <- as_tibble(BinomCI(all2$ppa*n,n,conf.level=0.95,method="clopper-pearson"))
allci <- bind_cols(all2,cpi)

write.csv(allci,"model_output.csv")




# Generate plots

# Print plot of pool size vs. PPA, panel grid of % samples with Ct > LoD, color by proportion test pos
pdf("plot.ppa.above.colprev.pdf",width=8, height=5)
ggplot(allci, aes(x=pool, y=tp/(tp + fn), color=as.factor(prevalence))) + 
  geom_line() +
  facet_wrap(vars(factor(above.llod, labels=c("5% samples Ct>LoD","10% samples Ct>LoD","15% samples Ct>LoD","20% samples Ct>LoD","25% samples Ct>LoD","30% samples Ct>LoD")))) +
  theme_bw() + 
  scale_color_brewer(palette="RdYlBu", name="Prevalence") + 
  xlab("Pool size") + ylab("Expected positive percent agreement") +
  scale_x_continuous(limits=c(1,20), breaks=seq(0,20,2)) + scale_y_continuous(limits=c(0.5,1),breaks=seq(0.5,1,.1)) +
  guides(color=guide_legend(title="Proportion of\ntests positive"))
dev.off()

# Print plot of pool size vs. PPA, panel grid of proportion test pos, color by % samples with Ct > LoD
pdf("plot.ppa.prev.colabove.pdf",width=8, height=5)
ggplot(allci, aes(x=pool, y=tp/(tp + fn), color=as.factor(above.llod))) + 
  geom_line() +
  facet_wrap(vars(factor(prevalence, labels=c("0.1% Tests positive", "1% Tests positive", "3% Tests positive", "5% Tests positive", "10% Tests positive","15% Tests positive")))) +
  theme_bw() + 
  scale_color_brewer(palette="RdYlBu", name="Prevalence") + 
  xlab("Pool size") + ylab("Expected positive percent agreement") +
  scale_x_continuous(limits=c(1,20), breaks=seq(0,20,2)) + scale_y_continuous(limits=c(0.5,1),breaks=seq(0.5,1,.1)) +
  guides(color=guide_legend(title="Proportion\nsamples\nCt > LoD"))
dev.off()

# Print plot of pool size vs. average tests/sample, panel grid of %Ct > LoD, color by proportion test pos
pdf("plot.eff.above.colprev.pdf",width=8, height=5)
ggplot(allci, aes(x=pool, y=tests.per.sample,color=as.factor(prevalence))) + 
  geom_line() +
  facet_wrap(vars(factor(above.llod, labels=c("5% samples Ct>LoD","10% samples Ct>LoD","15% samples Ct>LoD","20% samples Ct>LoD","25% samples Ct>LoD","30% samples Ct>LoD")))) +
  theme_bw() + 
  scale_color_brewer(palette="RdYlBu", name="Prevalence") + 
  xlab("Pool size") + ylab("Average tests per sample") + 
  scale_x_continuous(limits=c(1,20), breaks=seq(0,20,2)) + scale_y_continuous(trans="reverse", limits=c(1.2,0),breaks=seq(0,1.2,0.2)) +
  guides(color=guide_legend(title="Proportion of\ntests positive"))
dev.off()

# Print plot of pool size vs. PPA, gridded showing PPA changes with %>LoD but is constant at different LoDs
grid <- read.csv("model_output_varylod.varydist.csv")
grid <- as_tibble(grid)
pdf("AppendixFigure2.pdf",width=9, height=8)
ggplot(grid %>% filter((above.llod=="0.05"|above.llod=="0.15"|above.llod=="0.25"),(limit==34|limit==35|limit==36)), aes(x=pool, y=tp/(tp + fn), color=as.factor(prevalence))) + 
  geom_line() + 
  facet_grid(rows=vars(factor(limit,levels=c(36,35,34),labels=c("LoD Ct=36","LoD Ct=35","LoD Ct=34"))),cols=vars(factor(above.llod,labels=c("5% Ct>LoD","15% Ct>LoD","25% Ct>LoD")))) + 
  theme_bw() + 
  scale_color_brewer(palette="RdYlBu", name="Prevalence") + 
  xlab("Pool size") + ylab("Expected positive percent agreement") + 
  scale_x_continuous(limits=c(1,20),breaks=seq(0,20,2)) + scale_y_continuous(limits=c(0.5,1),breaks=seq(0.5,1,.1)) + 
  guides(color=guide_legend(title="Proportion of\ntests positive"))
dev.off()

