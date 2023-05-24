library(tidyverse)
library(FactoMineR)
library(vcd)
library(factoextra)

pheno_remerged_data <- read.csv('./data/Re_Merged_data_fill.csv', sep = ',')

pheno_remerged_data_2 <- pheno_remerged_data %>% 
  separate(col = 8, into = c('month', 'day', 'year'), sep = '\\/') %>%
  unite(col = 'year_location', c(7, 10)) %>%
  dplyr::select(-c('month', 'day')) %>%
  column_to_rownames(var = 'plot_id') %>%
  mutate(entry = str_replace_all(entry, "~", "-"))


pheno_remerged_data_2$trt <- as.factor(pheno_remerged_data_2$trt)
pheno_remerged_data_2$entry <- toupper(pheno_remerged_data_2$entry)
pheno_remerged_data_2$trt <- toupper(pheno_remerged_data_2$trt)
pheno_remerged_data_2$trt <- as.factor(pheno_remerged_data_2$trt)
pheno_remerged_data_2$year_location <- toupper(pheno_remerged_data_2$year_location)

str(pheno_remerged_data_2)
pheno_remerged_data_2$trt <- as.factor(pheno_remerged_data_2$trt)
pheno_remerged_data_2$entry <- as.factor(pheno_remerged_data_2$entry)

rows_untr <- as.data.frame(which(pheno_remerged_data_2=="UNTREATED",arr.ind=TRUE))$row



ggplot(pheno_remerged_data_2, aes(x = trt, y = yield)) + 
  geom_boxplot() + 
  facet_wrap(~ year_location) + 
  ggtitle('Yields of treated and untreated genotypes with instecides') + 
  xlab('Treated/Untreated') + 
  ylab('Yields') + 
  theme_classic()

res.pca <- PCA(pheno_remerged_data_2[,-c(1:6)], quali.sup = 1, 
               ind.sup = rows_untr, scale.unit=TRUE, ncp=5, graph=T)

res.pca$ind.sup

eig.val <- get_eigenvalue(res.pca)
eig.val

fviz_pca_biplot(res.pca, repel = TRUE, axes = 2:3)


x11()
fviz_pca_biplot(res.pca, 
                # Individuals
                geom.ind = "point",
                # fill.ind = pheno_remerged_data_2$trt,
                col.ind = "black",
                pointshape = 21, pointsize = 2,
                palette = "jco",
                addEllipses = FALSE,
                # Variables
                alpha.var ="contrib", col.var = "contrib",
                # gradient.cols = "RdYlBu",
                
                legend.title = list(fill = "Species", color = "Contrib",
                                    alpha = "Contrib")
)

x11()
p <- fviz_pca_biplot(res.pca, col.ind.sup = "blue", repel = TRUE, geom.ind = 'point', axes = 2:3, 
                     gradient.cols = "RdYlBu",
                     
                     legend.title = list(fill = "Species", color = "Contrib",
                                         alpha = "Contrib"))
p <- fviz_add(p, res.pca$quali.sup$coord, color = "red", geom.ind = 'point', axes = 2:3)
p

x11()
p <- fviz_pca_biplot(res.pca, col.ind.sup = "blue", repel = TRUE, geom.ind = 'point', axes = 1:2, 
                     gradient.cols = "RdYlBu",
                     
                     legend.title = list(fill = "Species", color = "Contrib",
                                         alpha = "Contrib"))
p <- fviz_add(p, res.pca$quali.sup$coord, color = "red", geom.ind = 'point', axes = 1:2)
p

x11()
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

var <- get_pca_var(res.pca)
var

var$coord
var$cos2
var$contrib

x11()
fviz_pca_var(res.pca, col.var = "black")

x11()
fviz_pca_var(res.pca, axes=c(2,3), col.var = "black")

library("corrplot")
x11()
corrplot(var$cos2, is.corr=FALSE)

x11()
fviz_cos2(res.pca, choice = "var", axes = 1:2)

x11()
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

x11()
fviz_pca_var(res.pca, axes=c(2, 3), col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

x11()
fviz_pca_var(res.pca, alpha.var = "cos2")

x11()
corrplot(var$contrib, is.corr=FALSE)    

x11()
par(mfrow=c(1,2))
fviz_contrib(res.pca, choice = "var", axes = 1)
fviz_contrib(res.pca, choice = "var", axes = 2)

x11()
fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 10)

x11()
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

#The most important (or, contributing) variables on axes 2-3 can be highlighted on the correlation plot as follow:
x11()
fviz_pca_var(res.pca, axes=c(2,3), col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

res.km <- kmeans(var$coord, centers = 2)
grp2<- as.factor(res.km$cluster)

x11()
fviz_pca_var(res.pca, col.var = grp2, 
             palette = c("#0073C2FF","#AA98CB"),
             legend.title = "Cluster", axes = c(2, 3))

x11()
fviz_pca_biplot(res.pca, label="var", habillage = pheno_remerged_data_2$trt,
                addEllipses=TRUE, ellipse.level=0.95, select.ind = list(contrib = 500), axes = c(1,2))

x11()
fviz_pca_biplot(res.pca, label="var", habillage = pheno_remerged_data_2$trt,
                addEllipses=TRUE, ellipse.level=0.95, axes = c(2,3))

dd1 <- dist(ind$coord[,c(1:3)],method="euclidean")#square root of the square sum of values (vector norm)
## distance matrix from all original centered-scaled variables.
dd2 <- dist(scale(pheno_merged_data2[, c(9:12)]),method="euclidean")

# calculation of the Hierarchical Clustering Assignment (HCA) on the factorial coordinates
#classif1 <- hclust(dd1,method="average")               # average of the distance
classif1 <- hclust(dd1,method="complete")               # Maximization of the differences
#classif1 <- hclust(dd1,method="single")                # Minimization of differences
classif2 <-hclust(dd1,"ward.D")                           # Method of ward (varianceintra/variance inter) ; proche kmeans en philosophie


# calculation of the Hierarchical Clustering Assignment (HCA) on the original variables
#classif1 <- hclust(dd1,method="average")               # average of the distance
classif1b <- hclust(dd2,method="complete")              # Maximization of the differences
#classif1 <- hclust(dd1,method="single")                # Minimization of differences
classif2b <-hclust(dd2,"ward.D")                        # Method of ward (varianceintra/variance inter) ; proche kmeans en philosophie

# Simple graph for 'complete linkage'. For more control, use dend (see below)
x11()
plot(classif1,hang=-1,
     main='Factorial coordinates',
     sub='HCA, complete linkage')  
x11()
plot(classif1b,hang=-1,
     main='Original variables',
     sub='HCA, complete linkage') 
dend1 <- as.dendrogram(classif1)
dend2 <- as.dendrogram(classif2)
# Dendrogram plot
x11();par(mfrow=c(1,2), mar=c(2,2,2,6)) 
plot(dend1, nodePar=list(pch = 19,cex=.1, col = 4), edgePar=list(col="red",lwd=1,lty=1),horiz = TRUE, main=paste("Factorial Coordinates"),sub="Average link, distance: Euclidean",cex=0.01, cex.lab = 0.01)
plot(dend2, nodePar=list(pch = 19,cex=.1, col = 4), edgePar=list(col="blue",lwd=1,lty=1),horiz = TRUE, main=paste("Original variables"),sub="Average link, distance: Euclidean",cex=0.01, cex.lab = 0.01)


# Simple graph for 'ward'
x11()
plot(classif2,hang=-1,
     main='Factorial coord. ',
     sub='HCA, Ward')  # Factorial coordinates
x11()
plot(classif2b,hang=-1,
     main='Original variables',
     sub='HCA, Ward') # Original variables

# we can wonder how many groups we could define, based on the fact. coord.
# 4 groups for 'complete linkage'
groupeslmax <- cutree(classif1, k=4)

# 4 groups for Ward
groupesward <- cutree(classif2, k=4)


# We compare the groups given by the two methods
table(groupeslmax,groupesward)

#Let's plot the clusters on PCA biplots
x11()
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             #pointsize = "cos2",
             col.ind = as.factor(groupeslmax), # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "#868686FF"),
             addEllipses = TRUE, # Concentration ellipses. If you want confidence ellipses instead of concentration ellipses, use ellipse.type = "confidence".
             legend.title = "Groups", 
             select.ind = list(contrib = 1000)
)

x11()
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = as.factor(groupeslmax), # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "#868686FF"),
             addEllipses = TRUE, ellipse.type="confidence", # Concentration ellipses. If you want confidence ellipses instead of concentration ellipses, use ellipse.type = "confidence".
             legend.title = "Groups"
)

moys.cah <- aggregate(pheno_merged_data2[, c(9:12)],list(groupeslmax),FUN=mean)
moys.cah


moys.ward <- aggregate(pheno_merged_data2[, c(9:12)],list(groupesward),FUN=mean)
moys.ward

anadelai1 <- aov(pheno_merged_data2[, c(9:12)]$yield ~ as.factor(groupeslmax))
plot(anadelai1)
shapiro.test(residuals(anadelai1))
car::leveneTest(anadelai1)

x11()
b <- MASS::boxcox(anadelai1)

lambda <- b$x[which.max(b$y)]
lambda
anadelai1transfo <- aov((pheno_merged_data2[, c(9:12)]$yield^lambda -1) ~ as.factor(groupeslmax))
summary(anadelai1transfo)


shapiro.test(residuals(anadelai1transfo))
# Variance homogeneity of the ANOVA residuals
leveneTest(anadelai1transfo)
print(agricolae::SNK.test(anadelai1transfo,"as.factor(groupeslmax)"))
x11()
plot(agricolae::SNK.test(anadelai1transfo,"as.factor(groupeslmax)"))


x11()
par(mfrow=c(2,2),pty="s")
hist(byd,breaks="Sturges")
hist(yield,breaks="Sturges")
hist(Ndvi,breaks="Sturges")
hist(Ptht,breaks="Sturges")

#------------------------------------------------------------------------------#

########
###Graph of individuals
########
ind <- get_pca_ind(res.pca)
ind

head(ind$coord)
# Quality of individuals
head(ind$cos2)
head(ind$contrib)

#color individuals by their cos2 values:
x11() 
fviz_pca_ind(res.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) # Avoid text overlapping (slow if many points)
)


#To change both point size and color by cos2:
x11()
fviz_pca_ind(res.pca, col.ind = "cos2", pointsize = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07") , select.ind = list(contrib = 500), repel = TRUE)# Avoid text overlapping (slow if many points)

#To create a bar plot of the quality of representation (cos2) of individuals on the factor map
x11()
fviz_cos2(res.pca, choice = "ind")

# Total contribution on PC1 and PC2
x11()
fviz_contrib(res.pca, choice = "ind", axes = 1:2)

#-------------------------------------------------------------------------------------------------------------#

res.pca <- PCA(pheno_remerged_data_2[,-c(1:4)], quali.sup = c(1: 3),  scale.unit=TRUE, ncp=5, graph=T)

