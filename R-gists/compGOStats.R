# these three need to be installed
# do this in R first before you run this script
# source("http://bioconductor.org/biocLite.R")
# biocLite("GSEABase")
# biocLite("AnnotationDbi")
# biocLite("GOstats")

library("AnnotationDbi")
library("GSEABase")
library("GOstats")


# to run this script simply make you file Target_gene_list.txt and run it like
# R --no-save < compGOStats.R > compGOStats.out

# read in the table
# assumes a table that looks like this listing all the genes which are in the Target set, one per line
# GENE1
# GENE2
genetable <-read.table("Target_gene_list.txt",header=F,sep="\t",stringsAsFactors=F, quote="")
genes <- genetable$V1

# read in the name of all the genes to compute  the background frequencies
allgenes <- read.csv("MYORG.maker.proteins.list.txt",
	 header=F,stringsAsFactors=F,sep="\t",quote="")
universe <- allgenes$V2

# problem matching mode of this before
mode(universe)
mode(genes)

# this is the Gene -> GO associations
# here I am using the GOslim associated terms, you could also use the more detailed .go file which
# is the full GO terms
godat <- read.table("MYORG.maker.proteins.INTERPRO.go_slim",header=F);
goframeData <- data.frame(godat$V1, godat$V2, godat$V3)
goFrame <- GOFrame(goframeData,organism="MYORG")
goAllFrame=GOAllFrame(goFrame)

gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

params <- GSEAGOHyperGParams(name="Molecular Function, overrepresentation test",
          geneSetCollection=gsc,
	  geneIds = genes,
	  universeGeneIds = universe,
	  ontology = "MF",
	  pvalueCutoff = 0.05,
	  conditional = FALSE,
	  testDirection = "over")


OverMF <- hyperGTest(params)
summary(OverMF)
OverMF

paramsCC <- GSEAGOHyperGParams(name="Cellular Compartment, overrepresentation test",
          geneSetCollection=gsc,
	  geneIds = genes,
	  universeGeneIds = universe,
	  ontology = "CC",
	  pvalueCutoff = 0.05,
	  conditional = FALSE,
	  testDirection = "over")

OverCC <- hyperGTest(paramsCC)
summary(OverCC)
OverCC

paramsBP <- GSEAGOHyperGParams(name="Biological Process, overrepresentation test",
          geneSetCollection=gsc,
	  geneIds = genes,
	  universeGeneIds = universe,
	  ontology = "BP",
	  pvalueCutoff = 0.05,
	  conditional = FALSE,
	  testDirection = "over")

OverBP <- hyperGTest(paramsBP)
summary(OverBP)
OverBP



params <- GSEAGOHyperGParams(name="Molecular Function, underrepresentation test",
          geneSetCollection=gsc,
	  geneIds = genes,
	  universeGeneIds = universe,
	  ontology = "MF",
	  pvalueCutoff = 0.05,
	  conditional = FALSE,
	  testDirection = "under")


UnderMF <- hyperGTest(params)
summary(UnderMF)
UnderMF

paramsCC <- GSEAGOHyperGParams(name="Cellular Compartment, underrepresentation test",
          geneSetCollection=gsc,
	  geneIds = genes,
	  universeGeneIds = universe,
	  ontology = "CC",
	  pvalueCutoff = 0.05,
	  conditional = FALSE,
	  testDirection = "under")

UnderCC <- hyperGTest(paramsCC)
summary(UnderCC)
UnderCC

paramsBP <- GSEAGOHyperGParams(name="Biological Process, underrepresentation test",
          geneSetCollection=gsc,
	  geneIds = genes,
	  universeGeneIds = universe,
	  ontology = "BP",
	  pvalueCutoff = 0.05,
	  conditional = FALSE,
	  testDirection = "under")

UnderBP <- hyperGTest(paramsBP)
summary(UnderBP)
UnderBP

