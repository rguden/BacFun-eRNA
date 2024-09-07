#################
#Load libraries
library(ggplot2)
library(plyr) #for join_all
library(dplyr) #for %>%
library(phyloseq)

#################
#CHANGE DATA DIRECTORY
setwd("/home/genomics/ahaegeman/RNA-seq_sediment/Taxonomy_rRNA_maximicrobial")

#################
#READ DATA

#READ ALL COUNT FILES FROM THE SAMPLES (DERIVED FROM KRAKEN2 OUTPUT)
counts.191206.07<-read.table("191206_07_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191206.07)<-c("R2_612","taxid")
counts.191206.08<-read.table("191206_08_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191206.08)<-c("WR2_612","taxid")
counts.191206.13<-read.table("191206_13_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191206.13)<-c("WR4_612","taxid")
counts.191206.15<-read.table("191206_15_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191206.15)<-c("R5_612","taxid")
counts.191206.17<-read.table("191206_17_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191206.17)<-c("R3_612","taxid")
counts.191206.18<-read.table("191206_18_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191206.18)<-c("R6_612","taxid")
counts.191206.19<-read.table("191206_19_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191206.19)<-c("R4_612","taxid")
counts.191206.21<-read.table("191206_21_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191206.21)<-c("WR3_612","taxid")
counts.191206.22<-read.table("191206_22_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191206.22)<-c("WR5_612","taxid")
counts.191206.23<-read.table("191206_23_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191206.23)<-c("WR6_612","taxid")
counts.191213.05<-read.table("191213_05_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191213.05)<-c("R2_1312","taxid")
counts.191213.06<-read.table("191213_06_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191213.06)<-c("WR2_1312","taxid")
counts.191213.09<-read.table("191213_09_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191213.09)<-c("R3_1312","taxid")
counts.191213.10<-read.table("191213_10_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191213.10)<-c("WR3_1312","taxid")
counts.191213.11<-read.table("191213_11_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191213.11)<-c("R4_1312","taxid")
counts.191213.12<-read.table("191213_12_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191213.12)<-c("WR4_1312","taxid")
counts.191213.14<-read.table("191213_14_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191213.14)<-c("WR5_1312","taxid")
counts.191213.16<-read.table("191213_16_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191213.16)<-c("R5_1312","taxid")
counts.191213.20<-read.table("191213_20_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191213.20)<-c("R6_1312","taxid")
counts.191213.24<-read.table("191213_24_rawcounts.txt",sep="\t",header=F,blank.lines.skip=T)
colnames(counts.191213.24)<-c("WR6_1312","taxid")

#MERGE ALL COUNT FILES
counts.all<-Reduce(merge, list(counts.191206.07,counts.191206.08,counts.191206.13,counts.191206.15,counts.191206.17,counts.191206.18,counts.191206.19,
                               counts.191206.21,counts.191206.22,counts.191206.23,counts.191213.05,counts.191213.06,counts.191213.09,counts.191213.10,
                               counts.191213.11,counts.191213.12,counts.191213.14,counts.191213.16,counts.191213.20,counts.191213.24))
dfs<-list(counts.191206.07,counts.191206.08,counts.191206.13,counts.191206.15,counts.191206.17,counts.191206.18,counts.191206.19,
          counts.191206.21,counts.191206.22,counts.191206.23,counts.191213.05,counts.191213.06,counts.191213.09,counts.191213.10,
          counts.191213.11,counts.191213.12,counts.191213.14,counts.191213.16,counts.191213.20,counts.191213.24)
counts.all<-join_all(dfs, by='taxid', type='full') #do a full join (include all taxids seen in all count tables)

#replace all NA values by 0
counts.all[is.na(counts.all)] <- 0

#some rows are duplicated for an unknown reason. The counts are the same for all samples over the duplicate rows, so we decide to remove these duplicated rows.
counts.all.unique<-counts.all %>% distinct(taxid, .keep_all = TRUE)

#write merged table to file
write.table(counts.all.unique, file="counts.txt",sep="\t",col.names=T,row.names=F,quote=F)

#READ TAXONOMY FILE
#Remark: in the originally retrieved taxonomy, there were some aberrant characters (',#, % -- or ,) in some names, these were removed since these give trouble when reading into R
#Remark2: also make sure there are no duplicate rows in the data
taxonomy<-read.table(file="taxonomy.txt",sep="\t",blank.lines.skip=T,header=T,na.strings="")

#READ METADATA FILE
metadata<-read.table(file="metadata_for_R.txt",sep="\t",header=T,colClasses = "character")
#order dataframe to extract the wanted sample order for the plots, first by date, then experiment, then replicate
wanted_sample_order<-metadata[order( metadata[,4], metadata[,6], metadata[,5] ),  ]$full_name
levels(metadata$full_name)<-wanted_sample_order
#set rownames
rownames(metadata)<-metadata$full_name

#################
#PHYLOSEQ OBJECT CREATION

#reshape counts to be used in phyloseq object
counts.all.for.ps<-as.matrix(counts.all.unique  %>% dplyr::select(-taxid)) #remove taxid column and convert to matrix
rownames(counts.all.for.ps)<-counts.all.unique$taxid  #set taxid as rownames

#reshape taxonomy to be used in phyloseq object
taxonomy.for.ps<-as.matrix(taxonomy  %>% dplyr::select(-taxid))
rownames(taxonomy.for.ps)<-taxonomy$taxid  #set taxid as rownames

#make phyloseq object
ps<-phyloseq(otu_table(counts.all.for.ps, taxa_are_rows = TRUE),tax_table(taxonomy.for.ps),sample_data(metadata))

#subset phyloseq object by selecting samples by timepoint
ps.1312<-subset_samples(ps, date=="191213")
ps.0612<-subset_samples(ps, date=="191206")

#filter taxa
ps.1312.filt = filter_taxa(ps.1312, function(x) sum(x) >= 1, TRUE)
ps.0612.filt = filter_taxa(ps.0612, function(x) sum(x) >= 1, TRUE)

#write to output
ps.1312.filt.tableforoutput<-cbind("ASV"=rownames(tax_table(ps.1312.filt)),tax_table(ps.1312.filt),otu_table(ps.1312.filt))
write.table(ps.1312.filt.tableforoutput,file="ASVtable_191213.txt",sep="\t",row.names=F,col.names=T,quote=F)

ps.0612.filt.tableforoutput<-cbind("ASV"=rownames(tax_table(ps.0612.filt)),tax_table(ps.0612.filt),otu_table(ps.0612.filt))
write.table(ps.0612.filt.tableforoutput,file="ASVtable_191206.txt",sep="\t",row.names=F,col.names=T,quote=F)


#################
#AGGLOMERATE TAXA ON CERTAIN TAXONOMIC LEVEL

ps.taxglom.superkingdom<-tax_glom(ps, taxrank="superkingdom") #superkingdom level
ps.taxglom.phylum<-tax_glom(ps, taxrank="phylum") #phylum level
ps.taxglom.class<-tax_glom(ps, taxrank="class") #class level
ps.taxglom.order<-tax_glom(ps, taxrank="order") #order level

#################
#SELECT ONLY BACTERIA

#subset of bacteria
ps.bacteria<-subset_taxa(ps, superkingdom=="Bacteria")

#filter low abundant taxa
ps.bacteria.filt = filter_taxa(ps.bacteria, function(x) sum(x) >= 20, TRUE)

#check count table
otu_table(ps.bacteria.filt)
#check tax table
tax_table(ps.bacteria.filt)

################
#MAKE BARPLOTS OF PHYLOSEQ OBJECTS AND SUBSETS

#SUPERKINGDOM
#calculate percentages
ps.taxglom.superkingdom.frac <- transform_sample_counts(ps.taxglom.superkingdom, function(OTU) OTU/sum(OTU)*100)
#barplot
barplot<-plot_bar(ps.taxglom.superkingdom.frac, fill="superkingdom") + 
  theme(axis.text=element_text(size=8),legend.text=element_text(size=8),legend.key.size=unit(0.4, "cm")) + 
  geom_bar(aes(color=superkingdom, fill=superkingdom), stat="identity", position="stack", color="white", width=0.4, linetype=0) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend(ncol = 1)) #add this line to force the legend in 1 column
barplot$data$Sample <- factor(barplot$data$Sample, levels = wanted_sample_order) #change the order of the samples in the plot
pdf(file="Barplot_superkingdom.pdf",width=11.7,height=8.3)
barplot 
trash<-dev.off()

#BACTERIA ON PHYLUM LEVEL
#agglomerate all counts on phylum level
ps.bacteria.filt.taxglom.phylum<-tax_glom(ps.bacteria.filt, taxrank="phylum") #phylum level
#calculate percentages
ps.bacteria.filt.taxglom.phylum.frac <- transform_sample_counts(ps.bacteria.filt.taxglom.phylum, function(OTU) OTU/sum(OTU)*100)
#barplot
barplot<-plot_bar(ps.bacteria.filt.taxglom.phylum.frac, fill="phylum") + 
  theme(axis.text=element_text(size=8),legend.text=element_text(size=8),legend.key.size=unit(0.4, "cm")) + 
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack", color="white", width=0.6, linetype=0) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend(ncol = 2)) #add this line to force the legend in 2 columns
barplot$data$Sample <- factor(barplot$data$Sample, levels = wanted_sample_order) #change the order of the samples in the plot
pdf(file="Barplot_Bacteria_phylum.pdf",width=11.7,height=8.3)
barplot 
trash<-dev.off()

#BACTERIA ON ORDER LEVEL
#agglomerate all counts on order level
ps.bacteria.filt.taxglom.order<-tax_glom(ps.bacteria.filt, taxrank="order") #order level
#calculate percentages
ps.bacteria.filt.taxglom.order.frac <- transform_sample_counts(ps.bacteria.filt.taxglom.order, function(OTU) OTU/sum(OTU)*100)
#barplot
barplot<-plot_bar(ps.bacteria.filt.taxglom.order.frac, fill="order") + 
  theme(axis.text=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.3, "cm")) + 
  geom_bar(aes(color=order, fill=order), stat="identity", position="stack", color="white", width=0.6, linetype=0) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend(ncol = 3)) #add this line to force the legend in 3 columns
barplot$data$Sample <- factor(barplot$data$Sample, levels = wanted_sample_order) #change the order of the samples in the plot
pdf(file="Barplot_Bacteria_order.pdf",width=11.7,height=8.3)
barplot 
trash<-dev.off()

#MAKE VENN DIAGRAM OF SHARED TAXA

#Install required package
#install.packages("remotes")
#remotes::install_github("adrientaudiere/MiscMetabar")
#install.packages("venneuler")
library(MiscMetabar)
library(venneuler)

#plot venn diagram. If print_values=T, the numbers are shown on the graph; if print_values=F, the resulting object is a venndiagram in a ggplot object
#the minimum number of sequences that need to be present can be specified with min_nb_seq=0; every value below this will be considered 0.
venn.13.12.samples.values<-venn_phyloseq(ps.1312.filt, "full_name", min_nb_seq = 0, print_values = T) 
venn.13.12.samples.plot<-venn_phyloseq(ps.1312.filt, "full_name", min_nb_seq = 0, print_values = F) 

venn.13.12.experiment.plot<-venn_phyloseq(ps.1312.filt, "experiment", min_nb_seq = 0, print_values = F)

