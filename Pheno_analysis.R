library(tidyverse)
library(lme4)
library(lmerTest)
library(caret)

install.packages('igr')

rm(list=ls())

#------------------------------------------------------------------------------#
pev_rel_estimation <- function(model) {
  # Using Henderson's equations
  
  X <- getME(model,'X')
  Z <- getME(model,'Z')
  Y <- getME(model,'y')
  
  # model w/o genetic relationships
  varianceestimates <- as.data.frame(VarCorr(model))[,4]
  se <- varianceestimates[2]
  su <- varianceestimates[1]
  
  lambda <- se/su
  # Iu <- diag(length(levels(as.factor(pheno_merged_data2$trt))) * 
  #              length(levels(as.factor(pheno_merged_data2$entry))) +
  #              length(levels(as.factor(pheno_merged_data2$trt))) +
  #              length(levels(as.factor(pheno_merged_data2$entry))) +
  #              length(levels(as.factor(pheno_merged_data2$rep))) * 
  #              length(levels(as.factor(pheno_merged_data2$column))))
  
  XpX <- crossprod(X)
  XpZ <- crossprod(X, Z)
  ZpX <- crossprod(Z, X)
  ZpZ <- crossprod(Z)
  XpY <- crossprod(X, Y)
  ZpY <- crossprod(Z, Y)
  
  Iu <- diag(dim(ZpZ)[1])
  
  ## LHS
  LHS <- rbind(cbind(XpX, XpZ),
               cbind(ZpX, ZpZ + Iu * lambda))
  
  
  ## RHS
  RHS <- rbind(XpY, 
               ZpY)
  
  # Inverse of LHS
  InvLHS <- round(solve(LHS), digit=3) ; size <- dim(InvLHS)[[1]]
  PEVpre <- diag(InvLHS)[ 6 : size]
  
  #' and then corresponding values:
  PEV <- PEVpre * se
  SEP <- sqrt(PEV)
  REL <- 1 - PEV / su
  
  solutions <- as.data.frame(ranef(model)) 
  
  return(data.frame(entry = solutions[, 3],
                    BLUP = solutions[, 4],
                    PEV = PEV,
                    SEP = SEP,
                    REL = REL))
}
#------------------------------------------------------------------------------#



sum(is.na(pheno_merged_data2$entry))

rm(list = ls())

pheno_merged_data <- read.csv('./data/Re_Merged_data_fill.csv', sep = ',')
str(pheno_merged_data)
pheno_merged_data$location <- as.factor(pheno_merged_data$location)
pheno_merged_data$trt <- as.factor(pheno_merged_data$trt)
pheno_merged_data$entry <- toupper(pheno_merged_data$entry)
pheno_merged_data$trt <- toupper(pheno_merged_data$trt)
pheno_merged_data$location <- toupper(pheno_merged_data$location)

pheno_merged_data2 <- pheno_merged_data %>% 
  separate(col = 8, into = c('month', 'day', 'year'), sep = '\\/') %>%
  unite(col = 'year_location', c(7, 10)) %>%
  mutate(entry = str_replace_all(entry, "~", "-"))

str(pheno_merged_data2)
pheno_merged_data2$trt <- as.factor(pheno_merged_data2$trt)
pheno_merged_data2$entry <- as.factor(pheno_merged_data2$entry)


out <- irg(pheno_merged_data2[,c('entry', 'Ndvi')])

pheno_merged_data2 


x11()
chart.Correlation(pheno_merged_data2[,c(11:14)], histogram = TRUE, method = "pearson")

#------------------------------------------------------------------------------#


#########

MLM_yields_ayn <- lmerTest::lmer(yield ~ entry
                                 + (1|trt) + (1|trt:entry)  ## interactions with a random factor are random
                                 + (1|rep/block),  
                                 data = ayn_byd0_16)
MLM_byd_ayn <- lmerTest::lmer(byd ~ entry
                                 + (1|trt) + (1|trt:entry)  ## interactions with a random factor are random
                                 + (1|rep/block),  
                                 data = ayn_byd0_16)


MLM_byd_ayn <- lme4::lmer(byd ~ entry + (1|trt) + (1|trt:entry) + (1|rep/block), data = ayn_byd0_16)
MLM__ayn <- lme4::lmer(yield ~ entry + (1|trt) + (1|trt:entry) + (1|rep/block), data = ayn_byd0_16)
MLM_yields_ayn <- lme4::lmer(yield ~ entry + (1|trt) + (1|trt:entry) + (1|rep/block), data = ayn_byd0_16)

ayn_byd0_16$mean_ndvi <- mean(c(ayn_byd0_16$))

rel_df_2 <- pev_rel_estimation(MLM_yields_ayn)

summary(MLM_yields_ayn)

mean(pheno_merged_data2)

MLM_yields_2 <- lme4::lmer(yield ~ year_location + (1|trt) + (1|entry) + (1|trt:entry)  ## interactions with a random factor are random
                           + (1|rep/block) + (1|range/rep/block) + (1|column/rep/block), data = pheno_merged_data2)

MLM_byd_2 <- lme4::lmer(byd ~ year_location + (1|trt) + (1|entry) + (1|trt:entry)  ## interactions with a random factor are random
                        + (1|rep/block) + (1|range/rep/block) + (1|column/rep/block), data = pheno_merged_data2)

MLM_Ndvi_2 <- lme4::lmer(Ndvi ~ year_location + (1|trt) + (1|entry) + (1|trt:entry)  ## interactions with a random factor are random
                         + (1|rep/block) + (1|range/rep/block) + (1|column/rep/block), data = pheno_merged_data2)

MLM_Ptht_2 <- lme4::lmer(100 * Ptht ~ year_location + (1|trt) + (1|entry) + (1|trt:entry)  ## interactions with a random factor are random
                           + (1|rep/block) + (1|range/rep/block) + (1|column/rep/block), data = pheno_merged_data2)

pev_rel_estimation(MLM_yields_2)

rel_df_2 <- pev_rel_estimation(MLM_yields_2)
rel_df_2 <- cbind(rel_df_2, dplyr::select(pev_rel_estimation(MLM_byd_2), -entry))
# rel_df_2 <- cbind(rel_df_2, dplyr::select(pev_rel_estimation(MLM_Ndvi_2), -entry))
rel_df_2 <- cbind(rel_df_2, dplyr::select(pev_rel_estimation(MLM_Ptht_2), -entry))

rel_df_2 %>% filter(BLUP_yields == '')


colnames(rel_df_2) <- c("variety_name", "BLUP_yields","PEV_yields","SEP_prod",
                      "REL_yields","BLUP_byd","PEV_byd","SEP_byd","REL_byd",
                      "BLUP_Ptht","PEV_Ptht","SEP_Ptht","REL_Ptht")

write.csv(rel_df_2, file = 'new_blups_final.csv')

rel_def_2_sep <- rel_df_2 %>% separate(col = 1, into = c('trt', 'variety'), sep = ':')

write.csv(rel_df, file = 'new_reldf.csv')

x11()
ggplot(rel_def_2_sep, aes(x = c(), y = )) + 
  geom_boxplot() + 
  facet_wrap(~ trt)

cor_blues <- cor(rel_df[, c(2, 6, 10 ,14)])
corrplot::corrplot(cor_blues)

install.packages('PerformanceAnalytics')
library(PerformanceAnalytics)

x11()
chart.Correlation(rel_df_2[c(1:760), c(2, 6, 10)], histogram = TRUE, method = "pearson")

lsmeans::lsmeans(MLM_byd)

rel_def_2_sep <- rel_df_2[c(1:760), c(1, 2, 6, 10)] %>% separate(col = 1, into = c('trt', 'variety'), sep = ':')

rel_def_2_sep %>% filter(trt == 'TREATED') %>% arrange(desc(BLUP_yields)) %>% head(n = 5) %>% knitr::kable()

pheno_merged_data2.active <- pheno_merged_data2 %>% 
                          filter(trt == 'UNTREATED') %>% 
                          dplyr::select(c(1,11:13)) %>% 
                          column_to_rownames(var = 'plot_id')

res.pca <- PCA(rel_def_2_sep %>% dplyr::select(-c(variety)), ind.sup = 1:380, quali.sup = 1, graph=TRUE)
res.pca$quali.sup

fviz_pca_var(res.pca)

p <- fviz_pca_ind(res.pca, col.ind.sup = "blue", repel = TRUE, axes = c(2:3))
p <- fviz_add(p, res.pca$quali.sup$coord, color = "red", axes = c(2:3))
p

fviz_pca_ind(res.pca, habillage = 13,
             addEllipses =TRUE, ellipse.type = "confidence",
             palette = "jco", repel = TRUE) 

eig.val <- get_eigenvalue(res.pca)
eig.val

fviz_eig(res.pca, addlabels = TRUE)

var <- get_pca_var(res.pca)
var

library("corrplot")
corrplot(var$cos2, is.corr=FALSE)

fviz_pca_biplot(res.pca, 
                # Individuals
                geom.ind = "point",
                # fill.ind = iris$Species, col.ind = "black",
                pointshape = 21, pointsize = 2,
                palette = "jco",
                addEllipses = TRUE,
                # Variables
                alpha.var ="contrib", col.var = "contrib",
                gradient.cols = "RdYlBu",
                
                legend.title = list(fill = "Species", color = "Contrib",
                                    alpha = "Contrib")
)

#------------------------------------------------------------------------------#
