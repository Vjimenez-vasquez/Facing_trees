# Facing_trees
Facing_trees

A code compiled by Victor Jimenez-Vasquez - vr.jimenez.vs@gmail.com

![facing_trees](https://user-images.githubusercontent.com/89874227/228942739-a3fe1217-654a-4a74-ab98-5a26c2b622f3.png)

# The code 
```r
setwd("C:/Users/USUARIO/Documents/facing")
dir()

#install.packages("ape")#
#install.packages("phytools")#
#install.packages("dplyr")#

library(ape)
library(phytools)
library(dplyr)

snp <- read.tree("snp_RAxML_bestTree_rooted.nwk")
plot(snp,main="SNPs")

mls <- read.tree("mlst_RAxML_bipartitions_rooted.nwk")
plot(mls,main="MLST")

par(mfrow=c(1,2))
plot(snp,main="SNPs")
plot(mls,main="MLST")

association <- cbind(sort(snp$tip.label), sort(mls$tip.label))

## examples ##
cophyloplot(snp, mls,use.edge.length = T, col=c("red","blue"),length.line=2,assoc=association)
cophyloplot(snp, mls,use.edge.length = T, col=c("red","blue"),length.line=2,assoc=association, space=200, gap=10)

m <- read.csv("metadata.tsv", header=T, sep="\t")
dim(m)
head(m)
m1 <- m[,c(1,6)]

ct <- read.csv("metadata_colors.tsv", header=T, sep="\t")
dim(ct)
head(ct)

a <- 0
a1 <- 0
b <- 0
b1 <- 0
for (i in unique(ct$ST)){
  a <- ct[ct$ST == i , 5]
  b <- unique(a)
  a1 <- append(a1,b)
  b1 <- append(b1,i)
}

cl <- data.frame(st = b1, col = a1)
colors <- cl[2:14,]


## the loop ##
colp <- colors$col
a <- 0 
a1 <- 0 
b <- 0 
b1 <- 0
c <- 0
c1 <- 0

gen <- unique(colors$st)
for (i in 1:length(gen)){
  a <- rep(colp[i],length(m1[m1$ST == gen[i],1]))
  b <- m1[m1$ST == gen[i],1]
  c <- m1[m1$ST == gen[i],2]
  a1 <- append(a1,a)
  b1 <- append(b1,b)
  c1 <- append(c1,c)
}
colq <- data.frame(accession=b1,genotype=c1,col=a1)
colq

ass1 <- as.data.frame(association)
ass2 <- data.frame(order=1:87,accession=gsub("'","",ass1$V1))
ass2

cr <- merge(colq,ass2,by="accession",all.x = F)
cr
colors <- cr[order(cr$order),]

h1<-max(nodeHeights(snp))
h2<-max(nodeHeights(mls))

## plot 1 ##
cophyloplot(snp, mls, use.edge.length = F, col=colors$col, length.line=10, 
            assoc=association, space=100, gap=0.5, show.tip.label=F, lwd = 2.5, lty = 1)

st <- cl[2:14,1]

legend(10,85,legend=st,fill=colp,cex=0.65,title="Sequence Type")

## plot 2 ##
wasp.cophylo <- cophylo(force.ultrametric(snp), force.ultrametric(mls),assoc=association)
plot(wasp.cophylo,link.type="curved",link.lwd=2,
     link.lty="solid",link.col=colors$col,main="fasdf",use.edge.length=T,ftype="off", scale.bar=round(0.5*c(h1,h2),2))
legend(0,0,legend=gen,fill=colp,cex=0.5,title="plots")

## plot 3 ##
cotangleplot(snp,mls,lwd=2,tangle="tree1",use.edge.length=F,cex=0.7,type="phylogram",link.col=colors$col, assoc=association)
