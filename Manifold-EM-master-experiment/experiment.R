require('dbscan')
require('MASS')
require(ggplot2)
require('Rcpp')
require(mclust)
require(scatterplot3d)
require(fcd)
sourceCpp('.cpp')

##Read Data
load_image_file <- function(filename) {
  return_list = list()
  f = file(filename, 'rb')
  readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  return_list$n = readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  nrow = readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  ncol = readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  x = readBin(f, 'integer',n = return_list$n*nrow*ncol, size = 1, signed = F)
  print(x)
  return_list$x = matrix(x, ncol = nrow*ncol, byrow=T)
  close(f)
  return(return_list)
}

load_label_file <- function(filename) {
  f = file(filename, 'rb')
  readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  n = readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  the_label = readBin(f, 'integer', n = n, size = 1, signed = F)
  close(f)
  return(the_label)
}

imagename = 'dataset/mnist/t10k-images-idx3-ubyte'
labelname = 'dataset/mnist/t10k-labels-idx1-ubyte'
manifold_data = load_image_file(imagename)$x
label_data = load_label_file(labelname)

##Nearest Neighbour Graph Construction
knng = kNN(manifold_data, 10)

##Geodesic Distance Approximation using Dijkstra's Method
require(igraph)
g <- make_empty_graph() %>%
  add_vertices(nrow(manifold_data)) 
for(i in 1:nrow(manifold_data)){
  for(j in 1:10){
    g = g+edges(c(i,knng$id[i,j]),weight = knng$dist[i,j])
  }
}
pathdist = distances(g,v=V(g),to = V(g),mode ='all',algorithm = 'dijkstra')

##Implementation of our algorithm and other algorithms
##Manifold-EM Clustering 
cats = cats_EM(pathdist,seq(1,1000,100),3,10,5)
ggplot()+
  geom_point(aes(x = as.vector(xall[which(cats==0),1]),y = as.vector(xall[which(cats ==0),2]), color = 'red'))+
  geom_point(aes(x = as.vector(xall[which(cats == 1),1]),y = as.vector(xall[which(cats == 1),2]), color = 'green'))+
  geom_point(aes(x = as.vector(xall[which(cats == 2),1]),y = as.vector(xall[which(cats == 2),2]), color = 'blue'))+
  geom_point(aes(x = as.vector(xall[which(cats == 3),1]),y = as.vector(xall[which(cats == 3),2]), color = 'yellow'))+
  xlim(-5,5)+
  ylim(-5,5)

##EM Clustering 
fit_em = Mclust(manifold_data,G=10)
cats_std = fit_em$classification
ggplot()+
  geom_point(aes(x = as.vector(xall[which(cats_std==4),1]),y = as.vector(xall[which(cats_std ==4),2]), color = 'red'))+
  geom_point(aes(x = as.vector(xall[which(cats_std == 1),1]),y = as.vector(xall[which(cats_std == 1),2]), color = 'green'))+
  geom_point(aes(x = as.vector(xall[which(cats_std == 2),1]),y = as.vector(xall[which(cats_std == 2),2]), color = 'blue'))+
  geom_point(aes(x = as.vector(xall[which(cats_std == 3),1]),y = as.vector(xall[which(cats_std == 3),2]), color = 'yellow'))+
  xlim(-5,5)+
  ylim(-5,5)

# Spectral Clustering
spec = spectral.clustering(pathdist, K=10)
ggplot()+
  geom_point(aes(x = as.vector(xall[which(spec == 4),1]),y = as.vector(xall[which(spec == 4),2]), color = 'red'))+
  geom_point(aes(x = as.vector(xall[which(spec == 1),1]),y = as.vector(xall[which(spec == 1),2]), color = 'green'))+
  geom_point(aes(x = as.vector(xall[which(spec == 2),1]),y = as.vector(xall[which(spec == 2),2]), color = 'blue'))+
  geom_point(aes(x = as.vector(xall[which(spec == 3),1]),y = as.vector(xall[which(spec == 3),2]), color = 'yellow'))+
  xlim(-5,5)+
  ylim(-5,5)

# k-means Clustering
kmeans_all_result = kmeans(manifold_data, centers = 10)
kmeans_result = kmeans_all_result$cluster
ggplot()+
  geom_point(aes(x = as.vector(xall[which(kmeans_result == 4),1]),y = as.vector(xall[which(kmeans_result == 4),2]), color = 'red'))+
  geom_point(aes(x = as.vector(xall[which(kmeans_result == 1),1]),y = as.vector(xall[which(kmeans_result == 1),2]), color = 'green'))+
  geom_point(aes(x = as.vector(xall[which(kmeans_result == 2),1]),y = as.vector(xall[which(kmeans_result == 2),2]), color = 'blue'))+
  geom_point(aes(x = as.vector(xall[which(kmeans_result == 3),1]),y = as.vector(xall[which(kmeans_result == 3),2]), color = 'yellow'))+
  xlim(-5,5)+
  ylim(-5,5)


## check the true percentage of each cluster
summary(as.factor(label_data == cats))
# visualization
dim2_matrix = layout_with_fr(g)
colnames(dim2_matrix) <- c("x1", "x2")
newlabel <- as.factor(spec)
originlabel <- as.factor(label_data)
vis_matrix = cbind(dim2_matrix, newlabel, originlabel)
ggplot(as.data.frame(vis_matrix), aes(x = x1, y = x2, colour = as.factor(newlabel))) + 
  geom_point(shape = originlabel)


# input:origin_data, label_result
plot_algo_cluster <- function(origin_data, label_result){
  label1 <- cbind(rep(1, length(origin_data$x1)), rep(2, length(origin_data$x2)), 
                  rep(3, length(origin_data$x3)), rep(4, length(origin_data$x4)))
  label2 <- new_result$cate[,1000]
  plot_matrix = cbind(xall, label1, label2)
  colnames(plot_matrix) <- c("x", "y", "x1", "x2")
}


