snps<-read.table("SNPs.vcf",sep="\t",header=F,blank.lines.skip=TRUE,
                 comment.char = "#")
colnames(snps)<-c("chr","start","id","refallele","altallele","qual",
                  "filter","info","format")
summary(snps)
# put the chromosomes in the good order: chr1, chr2, chr22, chrX
goodChrOrder <- c(paste("Supercontig_1.",c(1:8),sep=""),"MT_CBS_6936")
snps$chr <- factor(snps$chr,levels=goodChrOrder)

# Plot the densities of snps in the bed file for each chr seperately
library(ggplot2)
snpDensity<-ggplot(snps) + 
geom_histogram(aes(x=start),binwidth=1e4) + # pick a binwidth that is not too small 
facet_wrap(~ chr,ncol=2) + # seperate plots for each chr, x-scales can differ from chr to chr
ggtitle("Density of SNPs") +
xlab("Position in the genome") + 
ylab("SNP density") + 
theme_bw() # I prefer the black and white theme

# save the plot to .pdf file
pdf("snp_density.pdf",800,1000)
print(snpDensity)
#dev.off()
