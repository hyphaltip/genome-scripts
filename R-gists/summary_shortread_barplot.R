size <- read.table("%INPUT%",header=T,sep="\t",row.names=1)

pdf("%INPUT%.pdf")
barplot(t(size),xlab="Read Length", ylab="Total # Reads", main="Reads mapped by size",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=size$V1,legend=T,col=rainbow(5, start=.1, end=.91),beside=F))


