# Correlation matrix: pairwise variables
## We will be looking at PCA from steps A-E using the R base,
## and then using prcomp()
##   the splitting of eigenvalues?
## anyways

# 0. Setup
library(tidyverse)
library(knitr)
library(datasets)
library(ggplot2)
install.packages("devtools")
pak::pak("arleyc/PCAtest", force = TRUE)
devtools::install_github("arleyc/PCAtest", force = TRUE)
library(pcatest)

# 2. Read in the frogs’ data as a dataframe."
fd_raw <- read.csv("C:/DLSU/25-26/[4] LBBBI27/activities/a09/Meneses_et.al_2026_Platymantis.csv")
class(fd_raw)

# 3. Subset your data to retain only the “Species” column and morphological
#    features (e.g., SVL, HL, SNL, etc.) from columns 10 through the last one
#    (i.e., 10:X)."
fd_indivs <- fd_raw[,c(1,10:34)]

# 4. Check the class and structure of the subset data. How many morphological
#    features are present in the dataset?"
class(fd_indivs)
ncol(fd_indivs)
##   given that one column is Spp., it's 25 morphological features

# 5. Ensure the dataframe has no NAs."
complete.cases(fd_indivs)
fd_indivs_g <- fd_indivs[complete.cases(fd_indivs),]

# 6. Scale the data using only the numeric columns,
#    then compute the pairwise correlations. Do not further subset the data.
fd_scaled <- scale(fd_indivs_g[,2:26])

fd_cor <- cor(fd_scaled)

# 7. Use the correlation matrix to calculate the eigenvectors and the proportion
#    of variance (expressed in percentage). How many principal components were
#    produced? What is the proportion of variance explained by components 1–4?
fd_eig <- eigen(fd_cor)

## Proportion of variance:
fd_propvar <- fd_eig$values/sum(fd_eig$values)
fd_propvar

fd_propvar_percent <- (fd_eig$values/sum(fd_eig$values))*100
fd_propvar_percent

## Proportion of variance explained by components 1-4
fd_propvar_percent_1to4 <- sum(fd_propvar_percent[1:4])

# 8. Now, create a 1x2 plot showing the screeplot and the PCA biplot. You can
#    use either base R or ggplot2. You can be creative here by coloring the
#    point shapes per species and adding a legend. Consider the high data-ink
#    ratio you learned from the lecture.

fd_eigvecs <- fd_eig$vectors
fd_scores <- fd_scaled %*% fd_eigvecs

## Plotting a elbow plot for Scree I guess
plot(
 1:25,
 fd_propvar_percent,
 xlab = "Eigenvalues",
 ylab = "% Variance explained",
 main = "Elbow plot",
 lwd = 2,
 cex = 2,
 pch = 20,
 type = "b"
)

## PCA Biplot
plot(
  fd_scores[,1],
  fd_scores[,2],
  xlab = "PC1",
  ylab = "PC2",
  main = "PCA Biplot",
  pch = 20,
  col = as.factor(fd_indivs$Species)
)
legend("topright",
       legend = levels(as.factor(fd_indivs_g$Species)),
       col = 1:2,
       pch = 19)

###########
# B
###########

# 1. Examine the variance explained by principal components using the
#    Broken-stick method and the PCAtest package. What is the difference between
#   these two methods for assessing the importance of principal components?


## Broken stick method
broken_stick <- function(p){
  sapply(1:p, function(k){
    sum(1 / (k:p)) / p
  })
}

fd_bs <- broken_stick(ncol(fd_indivs_g[,2:26]))
fd_bs

## Comparison

fd_varcompare <- data.frame(
  PC = paste0("PC", 1:25),
  Observed = fd_propvar,
  BrokenStick = fd_bs
)
fd_varcompare


## PCAtest Package
fd_pcatest <- pcatest(fd_scaled, 100, 100, 0.5, varcorr=FALSE, counter=FALSE, plot=TRUE)

