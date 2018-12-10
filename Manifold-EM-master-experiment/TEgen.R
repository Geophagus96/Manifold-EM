require('dbscan')
require('MASS')
require(ggplot2)
require('Rcpp')
require(mclust)
require(scatterplot3d)
require(fcd)
sourceCpp('manifoldKmod.cpp')

TEgenerator = function(mean = c(-2, 2, 2, 2, -2, -2, 2, -2), rad = 2, range = 2.2, nnode = 300){
  # function to generate testing data for our toy example
  # in the toy example we are generating four clusters in 2D and then bend it into 3D
  # mean: center of 4 clusters
  # rad: radius of each cluster
  # range: parameter for random generator
  # nnode: default number of points to generate
  
  x1 = runif(2*nnode, -range, range)
  x1 = matrix(x1, nnode, 2)
  radius2 = apply(x1^2, 1, sum) #removing outliers
  x1 = x1[which(radius2 <= rad^2),]
  x2 = runif(2*nnode, -range, range)
  x2 = matrix(x2, nnode, 2)
  radius2 = apply(x2^2, 1, sum)
  x2 = x2[which(radius2 <= rad^2),]
  x3 = runif(2*nnode, -range, range)
  x3 = matrix(x3, nnode, 2)
  radius2 = apply(x3^2, 1, sum)
  x3 = x3[which(radius2 <= rad^2),]
  x4 = runif(2*nnode, -range, range)
  x4 = matrix(x4, nnode, 2)
  radius2 = apply(x4^2, 1, sum)
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
  
  TEx = list()
  TEx$x1 = x1
  TEx$x2 = x2
  TEx$x3 = x3
  TEx$x4 = x4
  TEx$xall = xall
  return(TEx)
}

TEMP = function(TElist){
  # constructing manifold and ploting the toy example
  r = (max(TElist$xall[,1]) - min(TElist$xall[,1])+2)/(2*pi)
  theta = (TElist$xall[,1] - min(TElist$xall[,1]))/r
  x = r*cos(theta)
  y = r*sin(theta)
  z = TElist$xall[,2]
  manifold_data = cbind(x,y,z)
  
  scatterplot3d(x, y, z, color = c(rep(1,nrow(TElist$x1)), rep(2,nrow(TElist$x2)), rep(3,nrow(TElist$x3)), rep(4,nrow(TElist$x4))))
  return(manifold_data)
}

TEplot = function(TElist, TEmanifold, TEclust){
  O2 = ggplot()+
    geom_point(aes(x = TElist$x1[,1],y = TElist$x1[,2], color = 'red'))+
    geom_point(aes(x = TElist$x2[,1],y = TElist$x2[,2], color = 'green'))+
    geom_point(aes(x = TElist$x3[,1],y = TElist$x3[,2], color = 'blue'))+
    geom_point(aes(x = TElist$x4[,1],y = TElist$x4[,2], color = 'yellow'))+
    xlim(-5,5)+
    ylim(-5,5)
  
  cats = TEclust$cate
  clust = ggplot()+
    geom_point(aes(x = as.vector(TElist$xall[which(cats == 0), 1]), y = as.vector(TElist$xall[which(cats == 0), 2]), color = 'red'))+
    geom_point(aes(x = as.vector(TElist$xall[which(cats == 1), 1]), y = as.vector(TElist$xall[which(cats == 1), 2]), color = 'green'))+
    geom_point(aes(x = as.vector(TElist$xall[which(cats == 2), 1]), y = as.vector(TElist$xall[which(cats == 2), 2]), color = 'blue'))+
    geom_point(aes(x = as.vector(TElist$xall[which(cats == 3), 1]), y = as.vector(TElist$xall[which(cats == 3), 2]), color = 'yellow'))+
    xlim(-5,5)+
    ylim(-5,5)
  
  IG = ggplot()+
    geom_point(aes(x = as.vector(TElist$xall[c(TEclust$initials)+1,1]), y = as.vector(TElist$xall[c(TEclust$initials)+1,2]), color = 'initials'),shape = 25, fill = 'red',size= 1.55)+
    geom_point(aes(x = as.vector(TElist$xall[-(c(TEclust$initials)+1),1]), y = as.vector(TElist$xall[-(c(TEclust$initials)+1),2])),size = 0.45)+
    xlim(-5,5)+
    ylim(-5,5)
  
  plotlist = list()
  plotlist$O2 = O2
  plotlist$clusters =clust
  plotlist$initg = IG
  
  return(plotlist)
}
