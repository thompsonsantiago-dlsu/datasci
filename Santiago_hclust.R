# Setup
library(tidyverse)
library(knitr)
library(datasets)
library(factoextra)
library(cluster)
library(ggplot2)

######################
# Data initialization
######################

# 2. Read in the morphology data as a dataframe.
morph_data <- read.csv("C:/DLSU/25-26/[4] LBBBI27/activities/a07/morphology_data.csv")
class(morph_data)

# 3. Remove the sex and population columns.
morph_ambi <- morph_data[c(3:ncol(morph_data))]

# 4. Ensure that the dataframe has no NAs. Standardize with scale().
morph_good <- morph_ambi[complete.cases(morph_ambi),]
morph_scaled <- scale(morph_good)
rownames(morph_scaled) <- morph_data$Pop

#############################
# A. Hierarchical Clustering
#############################

# 5. Before performing hierarchical clustering, calculate the “Euclidean distance” 
#    and use it as your input.
euc_d <- dist(morph_scaled, method = "euclidean")


# 6. Perform hierarchical clustering using Ward’s and Single Linkage methods.

## 6a. Ward's
morph_ward <- hclust(euc_d, method = "ward.D2")

plot(morph_ward,
     main = "Ward's Method",
     labels = rownames(morph_scaled),
     cex = 0.5)

## 6b. Single Linkage
morph_sl <- hclust(euc_d, method = "single")

plot(morph_sl,
     main = "Single Linkage",
     labels = rownames(morph_scaled),
     cex = 0.6)


# 7. Create a for loop function to explore Ks from 2 to 10 for each method using
#    the Silhouette score and display the results with a plot.
morph_avg_sil_widths <- matrix(nrow = 9, ncol = 2)

for (k_val in 2:10) {
  ward_clusters <- cutree(morph_ward, k = k_val)
  ward_sil <- silhouette(ward_clusters, euc_d)
  morph_avg_sil_widths[(k_val-1),1] <- mean(ward_sil[,3])
  morph_avg_sil_widths[(k_val-1),1]
  
  sl_clusters <- cutree(morph_sl, k = k_val)
  sl_sil <- silhouette(sl_clusters, euc_d)
  morph_avg_sil_widths[(k_val-1),2] <- mean(sl_sil[,3])
  morph_avg_sil_widths[(k_val-1),2]
}

plot(
  2:10,
  morph_avg_sil_widths[,1],
  type = "b",
  lwd = 2,
  xlab = "Number of Clusters",
  ylab = "Average Silhouette Widths of morph data (Ward)"
)

plot(
  2:10,
  morph_avg_sil_widths[,2],
  type = "b",
  lwd = 2,
  xlab = "Number of Clusters",
  ylab = "Average Silhouette Widths of morph data (Single Linkage)"
)

# 8. Now, generate a dendrogram for each method,
#    highlighting the optimal K with colored boxes.
## Ward's
plot(morph_ward,
     main = "Ward's Method",
     labels = rownames(morph_scaled),
     cex = 0.5,
     )
rect.hclust(morph_ward, k = 2,
            border = c("red", "blue"))

## Sinkle Ligages
plot(morph_sl,
     main = "SInkle Linkages",
     labels = rownames(morph_scaled),
     cex = 0.5,
)
rect.hclust(morph_ward, k = 2,
            border = c("red", "blue"))