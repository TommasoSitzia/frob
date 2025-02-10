# Script used in the manuscript ..

# libraries

library(rgdal)
library(ClustGeo)
library(tmap)

####

setwd("..")  # Set your project root as the working directory

####

# input data
data = readOGR(dsn="data", layer="frob")

# geographic distances between municipalities centroids
D1 = as.dist(dist(coordinates(data),diag=T,upper=T))
# assign spatial weights
# neighbours
data.nb <- poly2nb(data)
A = nb2mat(data.nb,style="B")
diag(A) = 1
D1n = 1-A
D1n = as.dist(D1n)
# partitioning municipalities in homogeneous clusters based on forest expansion by black locust (frob) (see https://arxiv.org/pdf/1707.03897.pdf)
# black locust
D0.frob = dist(data@data$frob) # frob distance
# Ward dendrogram based on frob
tree.frob=hclustgeo(D0.frob)
plot(tree.frob,hang=-1,label=FALSE,xlab="",sub="",main="")
# visual inspection suggests to retain k = 4 clusters
rect.hclust(tree.frob, k=4, border=c(4, 3, 2, 1))
legend("topright", legend=paste("cluster", 1:4), fill=1:4, bty="n", border="white")
# cut the dendrogram to get the partition in 4 clusters
p4frob <- cutree(tree.frob, 4) 
plot(data, border=NA, col=p4frob) # plot italy
legend("topleft", legend=paste("cluster", 1:4), fill=1:4, bty="n", border="white")
# take into account neighborhood constraints
crfrobD1n <- choicealpha(D0.frob, D1n, range.alpha=seq(0, 1, 0.1), K=4, graph=TRUE) # takes time!
crfrobD1n$Q # proportion of explained pseudo-inertia
crfrobD1n$Qnorm # normalized proportion of explained pseudo-inertia
# D1 always very small
# for alfa = 0.1 loss of 8% and 50% increase in D1
# Ward dendrogram based on frob plus neighborhood constraints
tree.frobn <- hclustgeo(D0.frob, D1n, alpha=0.1)
# cut the dendrogram to get the partition in 4 clusters
p4frobn<- cutree(tree.frobn, 4)
plot(data, border=NA, col=p4frobn) # plot italy
legend("topleft", legend=paste("cluster", 1:4), fill=1:4, bty="n", border="white")
# Create frob map
tmfrob <- tm_shape(data) + 
  tm_polygons(col = "frob", 
              style = "fisher", 
              n = 5, 
              border.alpha = 0, 
              title = expression(f[rob])) +
  tm_layout(legend.outside = FALSE)
# Create fnat map
tmfnat <- tm_shape(data) + 
  tm_polygons(col = "fnat", 
              style = "fisher", 
              n = 5, 
              border.alpha = 0, 
              title = expression(f[nat])) +
  tm_layout(legend.outside = FALSE)
# Combine maps into one frame with subframe labels
combined_map <- tmap_arrange(
  tmfrob + tm_layout(title = "a)", title.position = c("left", "top")),
  tmfnat + tm_layout(title = "b)", title.position = c("left", "top")),
  ncol = 2
)
# Display the combined map
combined_map # Fig. 1 
# Plot the cluster map
tmfrobncl <- tm_shape(data) + 
  tm_fill("clp4frobn", 
          style = "cat", 
          border.alpha = 0, 
          palette = "Paired", 
          title = expression("Clusters of " * f[rob])) +  # Subscript in legend title
  tm_layout(frame = FALSE, title = NULL)  # No map title, clean layout
# Display the map
tmfrobncl # Fig. 2
# Creat cluster boxplot # Fig. 3
data$clp4frobn = as.factor(p4frobn)
boxplot (frob ~ clp4frobn, 
         data = data, 
         xlab = "Clusters", 
         ylab = "Forest expansion by black locust")
# analysis was limited to cluster 2-4 from clp4frobn
