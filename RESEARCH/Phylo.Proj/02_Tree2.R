library(Biostrings)
library(ggplot2)
library(ape)
library(tidyverse)
library(ggtree)

nwk <- system.file("./TREE/ORY_R_tree.nwk")

read.tree("./TREE/ORY_R_tree.nwk")


ggtree(tree)
