library(ape)
library(ggtree)
tree <- read.tree("ClustalOmegaTreeDataCAPN1toWormCalpains.ph")
plot(tree)

ggtree(tree)

ggtree(tree) +
  geom_treescale() +
  # theme_tree2() + # horizontal scale shows the amount of amino acid change
  geom_tiplab() # add tip label layer for a tree

?geom_treescale
