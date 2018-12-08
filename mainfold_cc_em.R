require('dbscan')
require(igraph)
require(Rcpp)
sourceCpp('ccdist.cpp')

knng = kNN(manifold_data,3)

g <- make_empty_graph() %>%
  add_vertices(nrow(xall)) 
for(i in 1:nrow(xall)){
  for(j in 1:3){
    g = g+edges(c(i,knng$id[i,j]),weight = knng$dist[i,j])
  }
}

ccs = clusters(g)

pathdist = distances(g,v=V(g),to = V(g),mode ='all',algorithm = 'dijkstra')
l2dist = as.matrix(dist(manifold_data,method = 'euclidean'))

ccdis = matrix(vector(length = (ccs$no)^2),ccs$no, ccs$no)
ccdis = ccdis + t(ccdis)
ccrep = matrix(vector(length = (ccs$no)^2),ccs$no, ccs$no)
for (i in 1:(ccs$no-1)){
  for (j in (i+1):ccs$no){
    temp_dist = as.matrix(l2dist[which(ccs$membership == i), which(ccs$membership == j)])
    min_pair = minccdist(temp_dist)
    ccdis[i,j] = temp_dist[(min_pair[1]+1),(min_pair[2]+1)]
    ccrep[j,i] = which(ccs$membership == j)[(min_pair[2]+1)]
    ccrep[i,j] = which(ccs$membership == i)[(min_pair[1]+1)]
  }
}


