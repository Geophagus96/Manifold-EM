require('dbscan')
require('MASS')
require(ggplot2)
require('Rcpp')
require(mclust)
require(scatterplot3d)
sourceCpp('manifoldEM.cpp')

mean1 = c(-2,2)
mean2 = c(2,2)
mean3 = c(-2,-2)
mean4 = c(2,-2)
rad = 2

x1 = runif(600,-2.2,2.2)
x1 = matrix(x1,300,2)
radius2 = apply(x1^2,1,sum)
x1 = x1[which(radius2 <= rad^2),]
x2 = runif(600,-2.2,2.2)
x2 = matrix(x2,300,2)
radius2 = apply(x2^2,1,sum)
x2 = x2[which(radius2 <= rad^2),]
x3 = runif(600,-2.2,2.2)
x3 = matrix(x3,300,2)
radius2 = apply(x3^2,1,sum)
x3 = x3[which(radius2 <= rad^2),]
x4 = runif(600,-2.2,2.2)
x4 = matrix(x4,300,2)
radius2 = apply(x4^2,1,sum)
x4 = x4[which(radius2 <= rad^2),]
x1[,1] = x1[,1]-2
x1[,2] = x1[,2]+2
x2[,1] = x2[,1]+2
x2[,2] = x2[,2]+2
x3[,1] = x3[,1]-2
x3[,2] = x3[,2]-2
x4[,1] = x4[,1]+2
x4[,2] = x4[,2]-2
xall = rbind(x1,x2,x3,x4)

ggplot()+
  geom_point(aes(x = x1[,1],y = x1[,2], color = 'red'))+
  geom_point(aes(x = x2[,1],y = x2[,2], color = 'green'))+
  geom_point(aes(x = x3[,1],y = x3[,2], color = 'blue'))+
  geom_point(aes(x = x4[,1],y = x4[,2], color = 'yellow'))+
  xlim(-5,5)+
  ylim(-5,5)

r = (max(xall[,1])-min(xall[,1])+2)/(2*pi)
theta = (xall[,1]-min(xall[,1]))/r
x = r*cos(theta)
y = r*sin(theta)
z = xall[,2]
manifold_data = cbind(x,y,z)
scatterplot3d(x,y,z,color = c(rep(1,nrow(x1)),rep(2,nrow(x2)),rep(3,nrow(x3)),rep(4,nrow(x4))))

knng = kNN(manifold_data,5)

require(igraph)
g <- make_empty_graph() %>%
  add_vertices(nrow(xall)) 
for(i in 1:nrow(xall)){
  for(j in 1:5){
    g = g+edges(c(i,knng$id[i,j]),weight = knng$dist[i,j])
  }
}
pathdist = distances(g,v=V(g),to = V(g),mode ='all',algorithm = 'dijkstra')


cats = cats_EM(pathdist,c(10,120,230,340),2,4,10)
ggplot()+
  geom_point(aes(x = as.vector(xall[which(cats==0),1]),y = as.vector(xall[which(cats ==0),2]), color = 'red'))+
  geom_point(aes(x = as.vector(xall[which(cats == 1),1]),y = as.vector(xall[which(cats == 1),2]), color = 'green'))+
  geom_point(aes(x = as.vector(xall[which(cats == 2),1]),y = as.vector(xall[which(cats == 2),2]), color = 'blue'))+
  geom_point(aes(x = as.vector(xall[which(cats == 3),1]),y = as.vector(xall[which(cats == 3),2]), color = 'yellow'))+
  xlim(-5,5)+
  ylim(-5,5)

fit_em = Mclust(manifold_data,G=4)
cats_std = fit_em$classification
ggplot()+
  geom_point(aes(x = as.vector(xall[which(cats_std==4),1]),y = as.vector(xall[which(cats_std ==4),2]), color = 'red'))+
  geom_point(aes(x = as.vector(xall[which(cats_std == 1),1]),y = as.vector(xall[which(cats_std == 1),2]), color = 'green'))+
  geom_point(aes(x = as.vector(xall[which(cats_std == 2),1]),y = as.vector(xall[which(cats_std == 2),2]), color = 'blue'))+
  geom_point(aes(x = as.vector(xall[which(cats_std == 3),1]),y = as.vector(xall[which(cats_std == 3),2]), color = 'yellow'))+
  xlim(-5,5)+
  ylim(-5,5)
