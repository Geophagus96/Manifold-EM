source('TEgen.R')
source('manifoldEMwrapper.R')

## generate toy example
# For merge method (stocahstic)
Toy = TEgenerator() #generate data
Toymd = TEMP(Toy)# generate manifold
Toymerge = Manifold_EM(Toymd, 2, 4, 4, 3, 2, 5)# run the alg (stochastic)
Toyplots = TEplot(Toy, Toymd, Toymerge)

Toyplots$O2
Toyplots$clusters
Toyplots$initg


Toymerge2 = Manifold_EM(Toymd, 2, 4, 4, 2, 3, 5)# run the alg (min dist with their neighbors)
Toyplots2 = TEplot(Toy, Toymd, Toymerge2)
Toyplots2$O2
Toyplots2$clusters
Toyplots2$initg

Toymerge3 = Manifold_EM(Toymd, 2, 4, 4, 2, 4, 5)# run the alg (max dist with their neighbors)
Toyplots3 = TEplot(Toy, Toymd, Toymerge3)
Toyplots3$O2
Toyplots3$clusters
Toyplots3$initg
