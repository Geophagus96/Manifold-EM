require('dbscan')
require(igraph)
require(Rcpp)
sourceCpp('manifoldEM.cpp')
sourceCpp('nnrank.cpp')
sourceCpp('ccdist.cpp')
sourceCpp('geodist.cpp')
sourceCpp('nnstomerge.cpp')
sourceCpp('nnmindistmerge.cpp')
sourceCpp('nnmaxdistmerge.cpp')

Manifold_EM = function(manifold_data, n_manifolds, knns, categories, max_iter, method, thresh){
  g = graph_construction(manifold_data, knns)
  
  nnrank = sigpointsampling(g$knns$id)
  gdm = gdmgenerator(g$g)#geodesic distance matrix generator
  
  if (method == 1){
    initp = initialpointsamp(g$knns$id,nnrank,categories)
    # greedy nn-rank
  }
  else if (method == 2){
    initp = IPstomerge(g$knns$id, gdm, nnrank, categories, thresh)
    # with node merge method (stochastically)
  }
  else if (method == 3){
    initp = IPmindistmerge(g$knns$id, gdm, nnrank, categories, thresh)
    # with node merge method (min distance)
  }
  else{
    initp = IPmaxdistmerge(g$knns$id, gdm, nnrank, categories, thresh)
    # with node merge method (max distance)
  }
  
  cats = list();
  cats$cate = cats_EM(gdm,initp,n_manifolds,categories,max_iter)
  cats$initials = initp
  return(cats)
}

graph_construction = function(manifold_data, knns){
  knng = kNN(manifold_data,knns)
  g <- make_empty_graph() %>%
    add_vertices(nrow(manifold_data)) 
  for(i in 1:nrow(manifold_data)){
    for(j in 1:knns){
      g = g+edges(c(i,knng$id[i,j]), weight = knng$dist[i,j])
    }
  }
  graph_info = list()
  graph_info$g = g
  graph_info$knns = knng
  return(graph_info)
}

gdmgenerator = function(g){
  ccs = clusters(g)
  if (ccs$no == 1){
    gdm = distances(g,v=V(g),to = V(g),mode ='all', algorithm = 'dijkstra')
  }
  else{
    pathdist = distances(g,v=V(g),to = V(g),mode ='all', algorithm = 'dijkstra')
    l2dist = as.matrix(dist(manifold_data,method = 'euclidean'))
    ccdis = matrix(vector(length = (ccs$no)^2),ccs$no, ccs$no)
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
    ccdis = ccdis + t(ccdis)
    gdm = geodist(pathdist,l2dist,ccdis,ccrep,ccs$membership)
  }
  return(gdm)
}
