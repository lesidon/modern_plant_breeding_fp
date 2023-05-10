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


grid.arrange(pP1XP1,pP1vsP2,pP1vsP3,pP1vsP4,pP1vsP5,ncol=5) #rajouter l'argument pour avoir la m?me ?chelle.

dMod_H1 <- lm(H ~ Block + GCA(Par1, Par2) + tSCA(Par1, Par2) +
                RGCA(Par1, Par2) + RSCA(Par1, Par2), data = df)
summary(dMod_H1)
anova(dMod_H1)