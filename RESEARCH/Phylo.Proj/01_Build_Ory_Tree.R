# Load needed Libraries
library(msa)
library(phangorn)
library(ape)
library(tidyverse)
library(ggtree)
library(DECIPHER)
library(seqinr)
library(castor)

# Read in DNA alignment (MEGA: MUscle: gap open= -400, gap extend= 0).
seqs<- readDNAMultipleAlignment("./DATA_PROJECT/cleaned_aligned_oryza(final).fas")

print(seqs)


#convert to phangorn format
phydat <- msaConvert(seqs,type="phangorn::phyDat")
print(phydat)
# Quick and dirty look at tree.
  dm <- dist.ml(phydat)
  treeUPGMA  <- upgma(dm)
  treeNJ  <- NJ(dm)

  plot(treeUPGMA, main="UPGMA")
  plot(treeNJ, main="NJ")
#UPGMA tree looks good.

# Maximum likelihood:

# preper with modelTest 
mt <- modelTest(phydat, tree=treeUPGMA) 

# choose best model from the table according to AICc. 
bestmodel <- mt$Model[which.min(mt$AICc)]
# best model:
bestmodel
env = attr(mt, "env")
fitStart = eval(get(bestmodel, env), env) 


# optimize for bes model picked (GTR+G)
fit = optim.pml(fitStart, rearrangement = "stochastic",
                optGamma=TRUE, optInv=TRUE, model="GTR") 
# apply (BS)bootsrap test w/ 500 samples:
set.seed(123)
bs = bootstrap.pml(fit, bs=500, optNni=TRUE)
# plot BS tree
treeBS <- plotBS(midpoint(fit$tree), bs, p = 30, type="p")

# add BS values to fit tree (only reads in fig tree).
fit$node.label <- bs

# ggtree a pretty tree
ORY<- ggtree(fit$tree,root.position =48, branch.length='none')+
  geom_tiplab(align = TRUE, linesize = 0.5)

ORY

ORY + coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 240, 6, 6))

#export tree
write.tree(fit$tree, file = "ORZY_TREE.nwk", append = FALSE,
           digits = 10, tree.names = FALSE)

                


