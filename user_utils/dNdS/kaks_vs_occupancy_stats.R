# Script to analyze ka/ks and occupancy of pan-genome sequences
# Tested with codon-aligned barley transcripts and A. thaliana CDS sequences
# P Vinuesa, Carlos P Cantalapiedra, B Contreras Moreira 2016

#setwd("/path/to/file")
library(dplyr)
library(multcomp)
library(gplots)


dfr <- read.table(file="kaks4R.tab", header=T, sep="\t", stringsAsFactors = FALSE)
str(dfr)
dfr$n_of_seqs <- as.factor(dfr$n_of_seqs)
str(dfr)

#------------#
# >>> EDA <<<#
#------------#
# 1, lets get some summary stats for the data; look for outliers
dfr %>% group_by(n_of_seqs) %>% summarise(min(omega), max(omega), mean(omega), 
                                          median(omega), sd(omega))

# very obvious in the table, and in the graph
# Ooops, some dN/dS > 50 or even > 100, i.e. 1 and 2 orders of magnitude larger than
# most maxima !!!
plot(omega ~ n_of_seqs, data=dfr)

# So, remove outliers; start gently, based on the values on the plot
# Note that conserved aligned sequences have omega = 0 and will be also excluded

dfr <- dfr %>% filter(omega > 0 & omega < 20)
plot(omega ~ n_of_seqs, data=dfr)

dfr <- dfr %>% filter(omega > 0 & omega < 2)
plot(omega ~ n_of_seqs, data=dfr)

dfr <- dfr %>% filter(omega > 0 & omega < 1.5)
plot(omega ~ n_of_seqs, data=dfr)

# Extra check: inspection of some alignments with omega values > 1.5
#4.23	10	5	my_137701_TR9142-c0_g1_i3_aln.fna
#26	45	10	my_19169_TR18123-c0_g1_i1_aln.fna
#59.82	45	10	my_67440_TR23435-c0_g1_i1_aln.fna
#4.72	6	4	my_25328_TR25835-c0_g1_i1_aln.fna
#112.33	15	6	my_172475_TR22069-c0_g1_i1_aln.fna
#1.82	6	4	my_95720_TR13780-c0_g1_i1_aln.fna
# They all have poorly aligned regions, usualy at the termini, or very few aligned positions

# ok, with omega < 1.5 we have now similar numbers of outliers per occupancy class
# notches are fine now
boxplot(omega ~ n_of_seqs, data=dfr, notch=T)



#-----------------#
# >>> ANALYSIS <<<#
#-----------------#

# 2. now we can use standard 1-way analysis of variance
fit <- aov(omega ~ n_of_seqs, data=dfr)
summary(fit)

opar <- par(no.readonly = FALSE)
par(mar=c(5,5,5,2))

plotmeans(dfr$omega ~ dfr$n_of_seqs, cex = 0.8, xlab="cluster occupancy", ylab="dN/dS", 
          main="Mean plot with 95% CI")

pdf(file="Mean_plot_with_95percent_CI.pdf")
par(mar=c(5,5,5,2))
plotmeans(dfr$omega ~ dfr$n_of_seqs, cex = 0.8, xlab="cluster occupancy", ylab="dN/dS", 
          main="Mean plot with 95% CI", las=2)
dev.off()
par(opar)

# I prefer the meanplots or the boxplots that will be plotted below, but here are the
# original TukeyHSD plots; bars for treatments not crossing 0 are significant
plot(TukeyHSD(fit))

# prepare a nicer display of Tukey's honestly significant differences (HSD) as parallel
# boxplots, with significance letters plotted ontop
pancolors <- c(rep('orange',11),'yellow','white')

opar <- par(no.readonly = FALSE)
par(mar=c(5,3,10,2))
tuk <- glht(fit, linfct=mcp(n_of_seqs="Tukey"))
plot(cld(tuk, level=0.5), col=pancolors)

legend("topright", inset=.05, c("shell","soft-core","core"), 
       fill=c('orange','yellow','white'), horiz=TRUE, cex=1)

par(opar)

# now to file
pdf(file="tukeyHSD_tests_as_boxplots.pdf")
par(mar=c(5,6,10,2))
tuk <- glht(fit, linfct=mcp(n_of_seqs="Tukey"))
plot(cld(tuk, level=0.5), col=pancolors, ylab="dN/dS", xlab="cluster occupancy")

legend("topright", inset=.05, c("shell","soft-core","core"), 
       fill=c('orange','yellow','white'), horiz=TRUE, cex=0.7)
dev.off()
par(opar)

# 3. we expect a significant negative correlation
cor(dfr$omega,as.numeric(dfr$n_of_seqs))
cor.test(dfr$omega,as.numeric(dfr$n_of_seqs))

sink(file="correlation_analysis_omega_vs_occupancy.out")
cor(dfr$omega,as.numeric(dfr$n_of_seqs))
cor.test(dfr$omega,as.numeric(dfr$n_of_seqs))
sink()

# 4. lets fit a simple linear model
fit.lm <- lm(dfr$omega ~ as.numeric(dfr$n_of_seqs))
summary(fit.lm)

sink(file="linear_model_summary.out")
summary(fit.lm)
sink()

# 4.1 diagnose it; lm is certainly not the best fit ...
par(mfrow = c(2,2))
plot(fit.lm)
dev.off()
par(opar)

# 4.2 plot the lm
plot(as.numeric(dfr$n_of_seqs), dfr$omega, xlab="cluster occcupancy", ylab="dN/dS",
     main="y = 0.476 -0.023x")
abline(fit.lm)
dev.off()


# 5. plot Tukey's HSD tests for the 1-way anova results as boxplots with the linear model 
# plotted ontop
par(mar=c(5,6,10,2))
tuk <- glht(fit, linfct=mcp(n_of_seqs="Tukey"))
plot(cld(tuk, level=0.5), col=pancolors, ylab="dN/dS", xlab="cluster occupancy")

legend("topright", inset=.05, c("shell","soft-core","core"), 
       fill=c('orange','yellow','white'), horiz=TRUE, cex=0.6)
abline(fit.lm)

par(opar)
dev.off()

# ok, now to file
pdf(file="tukeyHSD_tests_as_boxplots_with_fitted_lm.pdf")
par(mar=c(5,6,10,2))
tuk <- glht(fit, linfct=mcp(n_of_seqs="Tukey"))
plot(cld(tuk, level=0.5), col=pancolors, ylab="dN/dS", xlab="cluster occupancy")

legend("topright", inset=.05, c("shell","soft-core","core"), 
       fill=c('orange','yellow','white'), horiz=TRUE, cex=0.7)
abline(fit.lm)
dev.off()

par(opar)


### And now the boxplots
par(mar=c(5,5,2,2))
boxplot(omega ~ n_of_seqs, data=dfr, notch=T, col=pancolors, ylab="dN/dS", 
        xlab="cluster occupancy" )

legend("topright", inset=.05, c("shell","soft-core","core"), 
       fill=c('orange','yellow','white'), horiz=TRUE, cex=0.6)

abline(fit.lm)

par(opar)
dev.off()

# plot to file  with and without the lm
pdf(file="fig4_parallel_boxplots_with_fitted_lm.pdf")
par(mar=c(5,5,2,2))
boxplot(omega ~ n_of_seqs, data=dfr, notch=T, col=pancolors, ylab="dN/dS", 
        xlab="cluster occupancy" )

legend("topright", inset=.05, c("shell","soft-core","core"), 
       fill=c('orange','yellow','white'), horiz=TRUE, cex=0.8)

abline(fit.lm)

par(opar)
dev.off()

pdf(file="parallel_boxplots_without_fitted_lm.pdf")
par(mar=c(5,5,2,2))
boxplot(omega ~ n_of_seqs, data=dfr, notch=T, col=pancolors, ylab="dN/dS", 
        xlab="cluster occupancy" )

legend("topright", inset=.05, c("shell","soft-core","core"), 
       fill=c('orange','yellow','white'), horiz=TRUE, cex=0.8)

par(opar)
dev.off()
