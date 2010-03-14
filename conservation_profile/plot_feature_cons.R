kb1<-read.table("promoter_alignments/gb_sordaria.gene-EMBL.1.stats.dat",header=T,sep="\t");
kb2<-read.table("promoter_alignments/gb_sordaria.gene-EMBL.2.stats.dat",header=T,sep="\t");
kb3<-read.table("promoter_alignments/gb_sordaria.gene-EMBL.3.stats.dat",header=T,sep="\t");
kb4<-read.table("promoter_alignments/gb_sordaria.gene-EMBL.4.stats.dat",header=T,sep="\t");


cds<-read.table("feature_alignments/gb_sordaria.CDS.stats.dat",header=T,sep="\t");
intron<-read.table("feature_alignments/gb_sordaria.intron.stats.dat",header=T,sep="\t");

                                        # summarize the data
summary(kb1);
summary(kb2);
summary(kb3);
summary(kb4);
summary(cds);
summary(intron);

# only take where %ID > 60
clean_1kb<-subset(kb1,kb1$PW_ID > 60);
clean_2kb<-subset(kb2,kb2$PW_ID > 60);
clean_3kb<-subset(kb3,kb3$PW_ID > 60);
clean_4kb<-subset(kb4,kb4$PW_ID > 60);
clean_cds<-subset(cds,cds$PW_ID > 60);
clean_intron<-subset(intron,intron$PW_ID > 60);

# summarize the cleaned data
summary(clean_1kb);
summary(clean_2kb);
summary(clean_3kb);
summary(clean_4kb);
summary(clean_cds);
summary(clean_intron);

# feature data

# Percent ID plots
pdf("promoter_ID.pdf");
boxplot(clean_1kb$PW_ID,clean_2kb$PW_ID,clean_3kb$PW_ID,clean_4kb$PW_ID,main="Boxplot % Identity");
hist(clean_1kb$PW_ID,100,main="% Identity for 1kb 5'-Upstream");
hist(clean_2kb$PW_ID,100,main="% Identity for 2kb 5'-Upstream");
hist(clean_3kb$PW_ID,100,main="% Identity for 3kb 5'-Upstream");
hist(clean_4kb$PW_ID,100,main="% Identity for 4kb 5'-Upstream");

# do some KS-TESTs
# only the 1-vs-[2-4] are significant, 2,3,4 are equivalent
ks.test(clean_1kb$PW_ID,clean_2kb$PW_ID);
ks.test(clean_1kb$PW_ID,clean_3kb$PW_ID);
ks.test(clean_1kb$PW_ID,clean_4kb$PW_ID);
ks.test(clean_2kb$PW_ID,clean_3kb$PW_ID);
ks.test(clean_2kb$PW_ID,clean_4kb$PW_ID);
ks.test(clean_3kb$PW_ID,clean_4kb$PW_ID);

ks.test(clean_cds$PW_ID,clean_1kb$PW_ID);
ks.test(clean_cds$PW_ID,clean_2kb$PW_ID);
ks.test(clean_intron$PW_ID,clean_1kb$PW_ID);
ks.test(clean_intron$PW_ID,clean_2kb$PW_ID);
ks.test(clean_intron$PW_ID,clean_cds$PW_ID);


pdf("all_ID.pdf");
boxplot(clean_cds$PW_ID,clean_intron$PW_ID,clean_1kb$PW_ID,clean_2kb$PW_ID,
        clean_3kb$PW_ID,main="Boxplot Smac-Ncra Pairwise-ID");

hist(clean_cds$PW_ID,100,main="Pairwise % Identity for Ncra-Smac CDS");
hist(clean_intron$PW_ID,100,main="Pairwise % Identity for Ncra-Smac Intron");
hist(clean_1kb$PW_ID,100,main="Pairwise % Identity for Ncra-Smac 1kb 5'-Upstream");
hist(clean_2kb$PW_ID,100,main="Pairwise % Identity for Ncra-Smac 2kb 5'-Upstream");



