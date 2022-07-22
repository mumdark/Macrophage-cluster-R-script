#Activate R package
library(GSEABase)
library(GSVA)
library(limma)
library(pheatmap)
library(ggplot2)
library(NMF)

#load data
load("Original data.RData")
#read exp matrix "xx.exp"
rt=as.matrix(xx.exp)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=normalizeBetweenArrays(avereps(rt))
#read macrophage gene sets
keggSet = macrophage.gene.sets
#do the gsva analysis
kegg <- gsva(expr=as.matrix(rt), keggSet, kcdf="Gaussian",method = "gsva",parallel.sz=1)
#do the nmf analysis,Data is the GSVA score matrix of macrophage gene set filtered by single factor COX regression and difference analysis
res=nmf(exp(data), rank=2:10, nrun=10, seed=123)