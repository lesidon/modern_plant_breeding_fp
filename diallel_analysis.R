graphics.off()
rm(list = ls())
options(contrasts=c('contr.sum','contr.poly'))

library(Hmisc)
library(gplots)
library(lattice)
library(car)
library(agricolae)
library(MASS)
library(lmDiallel)
library(tidyverse)

hmwg <-read.table('./data/Table3_Diallele_HMWG.csv', sep=";", header=T, stringsAsFactors = TRUE)
regeneration_rate <-read.table('./data/Table2_Diallel_RegenerationRate.csv', sep=";", header=T, stringsAsFactors = TRUE)
describe(hmwg)
describe(regeneration_rate)

#boxlpots
x11()
boxplot(HMWG ~ Genotype, las=2, data = df)

tri<-sort(tapply(HMWG, Genotype, mean))
x11()
boxplot( HMWG ~ factor(Genotype, levels=names(tri)), las=2, xlab="Progeny")

### Other informative graph
library(ggplot2)
library(reshape2)
library(gridExtra)


x11()
P1XP1<- df[c(1:4),]
P1XP1.m<- melt(P1XP1)
pP1XP1<- ggplot(P1XP1.m,aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P1 vs P1") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))

P1vsP2<- df[c(21:24),]
P1vsP2.m<- melt(P1vsP2)
pP1vsP2<- ggplot(P1vsP2.m,aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P1 vs P2") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))

P1vsP3 <- df[c(25:28),]
P1vsP3.m <- melt(P1vsP3)
pP1vsP3<- ggplot(P1vsP3.m,aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P1 vs P3") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))

P1vsP4<- df[c(29:32),]
P1vsP4.m<- melt(P1vsP4)
pP1vsP4<- ggplot(P1vsP4.m,aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P1 vs P4") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))

P1vsP5<- df[c(33:36),]
P1vsP5.m<- melt(P1vsP5)
pP1vsP5<- ggplot(P1vsP5.m,aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P1 vs P5") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))

#------

pP2XP2<- ggplot(melt(hmwg[hmwg$Genotype == 'P2XP2',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P2 vs P2") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))
pP2XP1<- ggplot(melt(hmwg[hmwg$Genotype == 'P2XP1',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P2 vs P1") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))
pP2XP3<- ggplot(melt(hmwg[hmwg$Genotype == 'P2XP3',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P2 vs P3") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))
pP2XP4<- ggplot(melt(hmwg[hmwg$Genotype == 'P2XP4',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P2 vs P4") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))
pP2XP5<- ggplot(melt(hmwg[hmwg$Genotype == 'P2XP5',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P2 vs P5") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))


#------

pP3XP3<- ggplot(melt(hmwg[hmwg$Genotype == 'P3XP3',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P3 vs P3") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))
pP3XP1<- ggplot(melt(hmwg[hmwg$Genotype == 'P3XP1',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P3 vs P1") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))
pP3XP2<- ggplot(melt(hmwg[hmwg$Genotype == 'P3XP2',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P3 vs P2") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))
pP3XP4<- ggplot(melt(hmwg[hmwg$Genotype == 'P3XP4',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P3 vs P4") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))
pP3XP5<- ggplot(melt(hmwg[hmwg$Genotype == 'P3XP5',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P3 vs P5") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))


#-----

pP4XP2<- ggplot(melt(hmwg[hmwg$Genotype == 'P4XP2',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P4 vs P2") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))
pP4XP1<- ggplot(melt(hmwg[hmwg$Genotype == 'P4XP1',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P4 vs P1") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))
pP4XP3<- ggplot(melt(hmwg[hmwg$Genotype == 'P4XP3',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P4 vs P3") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))
pP4XP4<- ggplot(melt(hmwg[hmwg$Genotype == 'P4XP4',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P4 vs P4") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))
pP4XP5<- ggplot(melt(hmwg[hmwg$Genotype == 'P4XP5',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P4 vs P5") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))

#-----

pP5XP2<- ggplot(melt(hmwg[hmwg$Genotype == 'P5XP2',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P5 vs P2") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))
pP5XP1<- ggplot(melt(hmwg[hmwg$Genotype == 'P5XP1',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P5 vs P1") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))
pP5XP3<- ggplot(melt(hmwg[hmwg$Genotype == 'P5XP3',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P5 vs P3") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))
pP5XP4<- ggplot(melt(hmwg[hmwg$Genotype == 'P5XP4',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P5 vs P4") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))
pP5XP5<- ggplot(melt(hmwg[hmwg$Genotype == 'P5XP5',]),aes(x=Genotype,y=value)) + geom_boxplot(fill="black",colour="red") + ggtitle("P5 vs P5") + theme(plot.title = element_text(lineheight=.8, face="bold"),axis.text.x =element_text(size=7))


x11()
grid.arrange(pP1XP1, pP1vsP2, pP1vsP3, pP1vsP4, pP1vsP5,ncol=5)
grid.arrange(pP2XP1, pP2XP2, pP2XP3, pP2XP4, pP2XP5,ncol=5)
grid.arrange(pP3XP1, pP3XP2, pP3XP3, pP3XP4, pP4XP5,ncol=5)
grid.arrange(pP4XP1, pP4XP2, pP4XP3, pP4XP4, pP4XP5,ncol=5)
grid.arrange(pP5XP1, pP5XP2, pP5XP3, pP5XP4, pP5XP5,ncol=5)


# ANOVA: naive approach

anova_naive <- aov(HMWG ~ Genotype + Rep, data = hmwg)
summary(anova_naive)

## ANOVA did not provide us with the evidence that there is a significant difference between the genotypes

x11()
par(mfrow=c(2,2))
plot(anova_naive)

# Checking required conditions
shapiro.test(residuals(anova_naive))

# Genuine Levene's test
leveneTest(aov(HMWG ~ Genotype, data = hmwg), 'median')

x11()
boxcox(aov(HMWG ~ Genotype, data = hmwg))

## the transformation would be useless, because 1 in the CI

comparisons <- TukeyHSD(anova_naive) ## from base R
comparisons

print(SNK.test(anova_naive, "Genotype"))  ## from agricolae
print(HSD.test(anova_naive, "Genotype"))

## Post-hoc analysis
MultCompTest  <- SNK.test( anova_naive, trt = c("Genotype", "Rep"), console = TRUE ) 
x11()
plot(MultCompTest, las=2)  

# The Griffing's models (1956)-full diallel experiments (including selfs and reciprocals; mating scheme 1)

dMod2_G1 <- lm.diallel(HMWG ~ Femelle + Male, Block = Rep,
                       data = hmwg, fct = "GRIFFING1")
anova(dMod2_G1)
summary(dMod2_G1)

#In order to obtain the full list of genetical parameters
gh_G1 <- glht(linfct = diallel.eff(dMod2_G1))
summary(gh_G1, test = adjusted(type = "none")) 

#The "Hayman" model (type 1):

dMod2_H1 <- lm.diallel(HMWG ~ Femelle + Male, Block = Rep,
                       data = hmwg, fct = "HAYMAN1")
summary(dMod2_H1)
anova(dMod2_H1)

gh_H1 <- glht(linfct = diallel.eff(dMod2_H1))
summary(gh_H1, test = adjusted(type = "none")) 


