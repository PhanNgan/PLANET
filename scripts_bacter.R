library(ggplot2)
library(svglite)
library(cowplot)
library(dplyr)
library(plyr)
library(magrittr)
library(nlme)
library(emmeans)
library(phyloseq)
library(textshape)
library(tidyverse)
library(eulerr)
library(microbiome)
library(vegan)
library(MASS)
library(lme4)
library(car)
#library(lsmeans)
library(multcomp)
#library(multcompView)
library(stringr)
library(DAtest)
library(tidymodels)
#library(broom.mixed)
library(Hmisc)
library(corrplot)


# Colors and shapes for ggplot2
mycolors <- scale_color_manual(values = c("CA"="#1B9E77", "CT"="#D95F02", "IR64"="lightblue2", "IR504"="cornflowerblue", "Azucena"="khaki3", "Zhonghua"="gold"))
myshapes <- scale_shape_manual(values = c("CA"=3, "CT"=5, "IR64"=15, "IR504"=16, "Azucena"=17, "Zhonghua"=18))


#https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
setwd("D:/PLANET/Cambodia_2018_bacteria/in_R")


# Load abundance file from Qiime2 output
otu_mat<-read.table("Cambodia_2018_16S-table-FINAL-unrarefied.txt", header=TRUE, check.names=FALSE, row.names=NULL, sep="\t")


samples_df <-read.table("metadata.txt", header=TRUE, check.names=FALSE, row.names=NULL, sep="\t", colClasses = c('character','character', 'factor', 'factor',
                                                                                                                 'factor', 'factor')) 
tax_mat<-read.table(file="taxonomy_table.txt", header=TRUE, check.names=FALSE, row.names=NULL, sep="\t")

tax_mat<-read.table(file="taxonomy_table_Nematozoa_soil_highest_blast_identity.txt", header=TRUE, check.names=FALSE, row.names=NULL, sep="\t")


otu_mat <- otu_mat %>%
  tibble::column_to_rownames("OTU") 
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("OTU")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

data <- phyloseq(OTU, TAX, samples)
data
rank_names(data)
sample_variables(data)

sample_names(data)

#https://micca.readthedocs.io/en/latest/phyloseq.html
df1 <- data
df2 <- subset_taxa(df1, Domain=="Bacteria")
df2 <- subset_samples(df2, Practices !="NU")
###DESEQ2
library(microbiome)
pseq.filt <- core(df2, detection = 0, prevalence = .1)

install.packages("DESeq2")
library(DESeq2)
ds2 <- phyloseq_to_deseq2(pseq.filt, ~ Practices)

ds2$Practices <- relevel(ds2$Practices, ref = "CA")

dds = ds2[ rowSums(counts(ds2)) > 5, ]
cts = counts(dds)
geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))

dds = estimateSizeFactors(dds, geoMeans=geoMeans)

dds1 <- DESeq(dds)
res05 <- results(dds1, alpha=0.05)
summary(res05)
plotMA(res05, ylim=c(-2,2))
deseq.results <- as.data.frame(results(dds1))
deseq.results$taxon <- rownames(results(dds1))


library(knitr)
library(dplyr)
deseq.results <- deseq.results %>% arrange(pvalue, log2FoldChange)
knitr::kable(deseq.results %>% filter(padj < 0.05), digits = 5)
res = results(dds1, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
write.table(sigtab, "Deseq2_results_CA_CT.txt")


dataB <- read.delim("Deseq2_results_taxon_Root_CA_CT.txt", header=TRUE, check.names = FALSE)

mycolors <- scale_color_manual(values = c("WPS-2"="darkseagreen1", "Solibacterales"="deeppink3","Methylococcales"="#fabebe", "uncultured"="palevioletred1", "Burkholderiales"="black", "Acidobacteriales"="peachpuff", "Geobacterales"="green", "Rhizobiales"="#1f78b4", "Elsterales"="#000075", "Pedosphaerales"="#fdbf6f", "Haliangiales"="#ff7f00", "Myxococcales"="#ffff99", "Subgroup_2"="khaki", "Gammaproteobacteria_Incertae_Sedis"="#b15928", "Xanthomonadales"="red", "Sphingomonadales"="royalblue", "Enterobacterales"="#3cb44b", "Bryobacterales"="#ffe119", "Bacillales"="purple4", "Micropepsales"="seagreen4", "Obscuribacterales"="#911eb4", "Kapabacteriales"="skyblue4", "Anaerolineales"="#46f0f0", "Kryptoniales"="lightgoldenrod", "Exiguobacterales"="rosybrown1", "Nitrospirales"="saddlebrown", "Cytophagales"="lightcyan", "Vicinamibacterales"="mediumblue", "Subgroup_5"="aquamarine", "Chitinophagales"="#f032e6", "Pseudomonadales"="#bcf60c", "Ktedonobacterales"="yellow3"))

mycolors <- scale_color_manual(values = c("Burkholderiales"="#b2df8a","Xanthomonadales"="#e6194b","Enterobacterales"="#3cb44b","Rhizobiales"="#1f78b4","Acidobacteriales"="peachpuff","Sphingomonadales"="royalblue","Methylococcales"="#fabebe"))

mycolors <- scale_color_manual(values = c("Roots"="#b2df8a","Rhizosphere"="#e6194b"))


B <- ggplot(dataB, aes(x=reorder(TaxoN, log2FC), y=log2FC, color=Compartment)) + 
  geom_point(aes(size=2)) + theme_minimal() + xlab ("ESV's assignment at species level") +
  ylab("Log2(Fold Change of abundance)") + coord_flip() + mycolors + 
  theme(axis.text.y = element_text (face="italic"), axis.text = element_text(size=11), axis.title=element_text(size=14,face="bold")) + 
  geom_hline(yintercept=0, color="red") +
  guides(color = guide_legend(override.aes = list(size = 5) ) )

B


myshapes <- scale_shape_manual(values = c("IR64"=15, "IR504"=16, "Azucena"=17, "Zhonghua"=18)) #"CA"=3, "CT"=5,



B <- ggplot(dataB, aes(x=reorder(TaxoN, log2FC), y=log2FC, color=Order)) + 
  geom_point(aes(size=2)) + theme_minimal() + xlab ("ESV's assignment at species level") +
  ylab("Log2(Fold Change of abundance)") + coord_flip() + mycolors + 
  theme(axis.text.y = element_text (face="italic"), axis.text = element_text(size=11), axis.title=element_text(size=14,face="bold")) + 
  geom_hline(yintercept=0, color="red") +
  guides(color = guide_legend(override.aes = list(size = 5) ) )

B
####



######Differential taxonomic enrichments
# Subset the samples
df3 <- subset_samples(df2, Compartment=="Roots")
df3 <- subset_samples(df2, Compartment=="Rhizosphere")
# Subset the samples
compart <- subset_samples(df2, Compartment=="Roots")
compart <- subset_samples(df2, Compartment=="Rhizosphere")

ps_4_16S_Zhonghua <- subset_samples(df3, Cultivar=="IR504")
ps_4_16S_Zhonghua <- subset_samples(df3, Cultivar=="IR64")
ps_4_16S_Zhonghua <- subset_samples(df3, Cultivar=="Azucena")
ps_4_16S_Zhonghua <- subset_samples(df3, Cultivar=="Zhonghua")
# Trim low abundant features to retain features that are at least present in 25% of the samples (25% = present in 2 of over samples) and that have at least a total count of 10 across all samples

ps_4_16S_Zhonghua_f <- preDA(ps_4_16S_Zhonghua, min.samples = 2, min.reads = 10, min.abundance = 0)
compart_f <- preDA(compart, min.samples = 4, min.reads = 10, min.abundance = 0)
# Test from phyloseq object
test <- testDA(ps_4_16S_Zhonghua_f, predictor = "Practices", relative = FALSE) # Runs all methods and compares their performance # My data are absolute abundances/counts so relative = FALSE

test <- testDA(compart_f, predictor = "Practices", relative = FALSE) # Runs all methods and compares their performance # My data are absolute abundances/counts so relative = FALSE

summary(test)
plot(test)
saveRDS(test, "DAtest_16S_IR64.rds")

# Chose the best test according to the DAtest results:
## Power should be as high as possible
## Score should be as high as possible
## AUC should be as high as possible
## FPR should be lower than 0.05
## FDR should be as low as possible

# Make the test and save the results
IR504 <- DA.poi(ps_4_16S_Zhonghua_f, predictor = "Practices", relative = FALSE)
IR64 <- DA.poi(ps_4_16S_Zhonghua_f, predictor = "Practices", relative = FALSE)
Azucena <- DA.poi(ps_4_16S_Zhonghua_f, predictor = "Practices", relative = FALSE)
Zhonghua <- DA.poi(ps_4_16S_Zhonghua_f, predictor = "Practices", relative = FALSE)
Roots <- DA.poi(compart_f, predictor = "Practices", relative = FALSE)
Rhizosphere <- DA.poi(compart_f, predictor = "Practices", relative = FALSE)


# Import the results by cultivar
IR504<-read.table("DAtest_poi_IR504.txt", header=TRUE, check.names = FALSE, sep = "\t", row.names = 1)
IR64<-read.table("DAtest_poi_IR64.txt", header=TRUE, check.names = FALSE, sep = "\t", row.names = 1)
Azucena<-read.table("DAtest_poi_Azucena.txt", header=TRUE, check.names = FALSE, sep = "\t", row.names = 1)
Zhonghua<-read.table("DAtest_poi_Zhonghua.txt", header=TRUE, check.names = FALSE, sep = "\t", row.names = 1)
nrow(IR504)
nrow(IR64)
nrow(Azucena)
nrow(Zhonghua)

nrow(Roots)
nrow(Rhizosphere)

# Add a column to inform the cultivar
vec<-rep("IR504",194) # 629 for 16S, 139 for ITS, 21 for NEM
IR504["cultivar"] <- vec
vec<-rep("IR64",165) # 551 for 16S, 165 for ITS, 18 for NEM
IR64["cultivar"] <- vec
vec<-rep("Azucena",152) # 576 for 16S, 132 for ITS, 21 for NEM
Azucena["cultivar"] <- vec
vec<-rep("Zhonghua",186) # 574 for 16S, 141 for ITS, 18 for NEM
Zhonghua["cultivar"] <- vec


vec<-rep("Roots",277) # 576 for 16S, 132 for ITS, 21 for NEM
Roots["compartment"] <- vec
vec<-rep("Rhizosphere",485) # 574 for 16S, 141 for ITS, 18 for NEM
Rhizosphere["compartment"] <- vec
# Merge the tables
ALL_16S <- rbind(IR504,IR64, Azucena, Zhonghua)
ALL_16S <- rbind(Roots, Rhizosphere)
write.csv(ALL_16S,"DAtest_16S_allcultivars.csv")
# Filter the differentially abundant SVs based on p-value
SIG_16S <- ALL_16S[ALL_16S$pval.adj < 0.05, ]
write.csv(SIG_16S,"DAtest_16S_Rhizo_Roots_practices_sig.csv")
##Add species name base on blast to NCBI '16S RNA sequence bacteria and Archaea"
library(dplyr)
library(tidyverse)
setwd("D:/PLANET/Cambodia_2018_bacteria/in_R/DAtest/")
blast_hit <- read.table("Root_Rhizo_DAtest_CT_CA_1.txt", sep = "\t", header = T, quote="", stringsAsFactors = FALSE)
blast_hit <- blast_hit %>% rename("query_acc.ver"=V1, "accession.version"=V2, "identity"=V3, "alignment_length"=V4, "mismatches"=V5, "gap_opens"=V6, "q._start"=V7, "q._end"=V8, "s._start"=V9, "s._end"=V10, "evalue"=V11, "bit_score"=V12)
tax_id <- read.table("Root_Rhizo_DAtest_CT_CA_acc_taxid.txt", sep = "\t", header = T, quote="", stringsAsFactors = FALSE)
tax_id <- tax_id %>% rename("accession"=V1, "accession.version"=V2, "tax_id"=V3, "gi"=V4)
tax_lig <- read.csv("Root_Rhizo_DAtest_CT_CA_lig1.txt", sep = "\t", header = T, quote="", stringsAsFactors = FALSE) 
tax_table <- left_join(blast_hit, tax_id, by="accession.version")
tax_lig_table <- left_join(tax_table, tax_lig, by="tax_id")
write.table(tax_lig_table, "104_seqs_tax_lig_table.txt" , sep = "\t", quote=F, row.names = F, col.names = T)

tax_lig_table_2 <- tax_lig_table %>% 
  rowwise() %>% 
  group_by(query_acc.ver) %>% 
  dplyr::arrange(desc(identity)) %>% 
  slice_head(n=1) %>% 
  dplyr::select(query_acc.ver,tax_id,phylum,order,family,genus,species)

Datest_out <- read.csv("Root_Rhizo_DAtest_CT_CA_out.txt", sep = "\t", header = T, quote="", stringsAsFactors = FALSE) 

Final_table <- left_join(Datest_out, tax_lig_table_2, by="query_acc.ver")

write.table(Final_table, "Enrichment_DAtest_Roots_Rhizo_log2FC_with_NCBI_lig_table.txt" , sep = "\t", quote=F, row.names = F, col.names = T)

# Remove the NA and arrange the columns to show in the graph

END_16S <- Final_table

filter1 <- END_16S %>% rowwise() %>% 
  group_by(species,cultivar) %>% 
  dplyr::summarize(Prac_both_ordering=paste(ordering,collapse = "/"),species) %>% 
  slice_head(n=1)

write.csv(filter1,"DAtest_filter1.csv")
END_16S <- read.table("DAtest_16S_allcultivars_arranged_filter2.txt",sep='\t', header=TRUE, row.names = , check.names=FALSE)
# Draw the graph
END_16S$cultivar<-ordered(END_16S$cultivar, levels=c('IR504', 'IR64', 'Azucena', 'Zhonghua'))


mycolors <- scale_color_manual(values = c("IR64"="blue", "IR504"="red", "Azucena"="green", "Zhonghua"="black"))

myshapes <- scale_shape_manual(values = c("IR64"=15, "IR504"=16, "Azucena"=17, "Zhonghua"=18)) #"CA"=3, "CT"=5,
mycolors <- scale_color_manual(values = c("Burkholderiales"="#b2df8a","Xanthomonadales"="#e6194b","Sphingomonadales"="#3cb44b","Enterobacterales"="#1f78b4","Methylococcales"="#fabebe"))

## 16S and ITS
taxaenriched_16S <- ggplot(END_16S, aes(x=reorder(species, log2FC), y=log2FC, color=cultivar)) + 
  geom_point(aes(size=1)) + theme_minimal() + xlab ("Species") + ylab("Log2FC(absolute abundance)") + 
  coord_flip() + 
  theme(axis.text.y = element_text (face="italic"), axis.text = element_text(size=11), axis.title=element_text(size=14,face="bold")) + 
  geom_hline(yintercept=0, color="red") + mycolors +
  guides(color = guide_legend(override.aes = list(size = 5) ) )



taxaenriched_16S

# NEM
taxaenriched_NEM <- ggplot(END_NEM, aes(x=reorder(Family, log2FC), y=log2FC, shape = cultivar, color = functional_guild)) + geom_point() + theme_minimal() + xlab ("Family") + ylab("Log2FC(absolute abundance)") + coord_flip() + theme(axis.text.y = element_text (face="italic"), axis.text = element_text(size=11), axis.title=element_text(size=14,face="bold")) + geom_hline(yintercept=0, color="red") + myshapes + mycolors
taxaenriched_NEM
# Save individual plots
svglite(file="DAtest_NEM.svg", width = 6, height = 3, bg = "white", pointsize = 16)
taxaenriched_NEM
dev.off()
# width = 6, height = 3 for ITS













ps2 <- data

df2 <- subset_samples(df1, Compartment=="Rhizosphere")

df2 <- subset_samples(df2, Practices=="CT")

df3 = as(otu_table(df2), "matrix")
df3 = as.data.frame(df3)

write.table(df3,"OTU_Rhizo_CT",sep="\t")



ps2 <- prune_samples(sample_sums(df1)>=10, df1)

filter <- phyloseq::genefilter_sample(ps2, filterfun_sample(function(x) x >= 10))

ps2 <- prune_taxa(filter, ps2)


melted_data <- psmelt(ps2)
melted_data <- melted_data %>% filter(Abundance>=10)

write.table(melted_data,"melted_OTU_taxo_table2.txt", sep="\t")

length(unique(melted_data$Phylum))

table1 <-read.table("phylum_table_1", sep = "\t", header = F)
table2 <-read.table("phylum_table_2", sep = "\t", header = F, comment.char = "")
colnames(table1) <- c("OTU", "Practices", "Compartment", "Phylum")
colnames(table2) <- c("OTU", "Practices_col", "Compartment_col", "Phylum_col")
table3 <-  merge(table1, table2, by = "OTU")
write.table(table3,"phylum_table_3.txt",sep="\t", row.names = F, quote = F)

df2 <- subset_samples(df1, Compartment=="Rhizosphere")

ps1 <- subset_taxa(df2, Domain=="Bacteria")
ps2 <- subset_samples(ps1, Practices != "NU")
ps3 <- subset_samples(ps2, Practices == "CA")


ps3 <- filter_taxa(ps3, function(x) sum(x) > 0, TRUE) # Filter SVs in the taxo
check_taxo=ps3@tax_table@.Data


OTU_prac_CA = as(otu_table(ps3), "matrix")



OTU_prac_CA = as.data.frame(OTU_prac_CA)


write.table(OTU_prac_CA,"OTU_prac_CA_Rhizo.txt", sep="\t")

BiocManager::install(c("DECIPHER","phangorn","ape"))
install.packages("ape")
library(ape)
library(DECIPHER)
library(phangorn)
library(Biostrings)
n
seqs <- readDNAStringSet("OTU_prac_CT_Root.fasta")
# names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(seqs, anchor=NA,verbose=FALSE)

phangMat <- as(alignment, "matrix")
dm <- dist.ml(phyDat(phangMat, type = "DNA"))

treeNJ <- NJ(dm) 

myBoots <- boot.phylo(treeNJ, phangMat, function(xx) NJ(dist.ml(phyDat(xx, type = "DNA"))), 
                      B = 10, block = 1,
                      trees = T)
plot.phylo(treeNJ, "fan", cex = .5, show.tip.label = F)
nodelabels((myBoots$BP)/10*100, col = "blue", bg = NA, frame = "none", cex = .5)

total = median(sample_sums(ps2))
standf = function(x, t=total) round(t * (x / sum(x)))
ps2 = transform_sample_counts(ps2, standf)





plot_bar(ps2, x="Practices", fill = "Phylum", facet_grid = ~Cultivar) +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") 

plot_bar(ps2, x="Practices", fill = "Family") +
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") 



ps2 <- subset_taxa(ps2, !Family %in% c("Pratylenchidae"))
ps2 <- subset_taxa(ps2, Family %in% c("Anguinidae", "Meloidogynidae", "Pratylenchidae","Telotylenchidae"))
ps2 <- subset_taxa(ps2, Family %in% c("Anguinidae", "Meloidogynidae","Telotylenchidae"))



##merge by practice
physeq_merged_ITS_pr <- merge_samples(ps2, "Compartment")

physeq_merged_ITS_pr <- merge_samples(ps2, "Practices")

physeq_merged_ITS_pr <- merge_samples(ps2, "Cultivar")

#merge by practice and cultivar
var1 = as.character(get_variable(ps2, "Practices"))
var2 = as.character(get_variable(ps2, "Cultivar"))
sample_data(ps2)$PracCul <- mapply(paste0, var1, var2, collapse = "_")
physeq_merged_ITS_pr <- merge_samples(ps2, "PracCul")

var1 = as.character(get_variable(ps2, "Practices"))
var2 = as.character(get_variable(ps2, "Compartment"))
sample_data(ps2)$ComPrac <- mapply(paste0, var1, var2, collapse = "_")
physeq_merged_ITS_pr <- merge_samples(ps2, "ComPrac")



glom_ITS <- tax_glom(physeq_merged_ITS_pr, taxrank = 'Phylum') # 'Phylum' for 16S and ITS and 'Family' for NEM when TAXA # 'GROUP1' or 'GROUP2' for 16S, 'Trophic_Mode' for ITS and 'functional_guild' (or 'cp') for NEM when GUILDS


# Give abundance in %
glom2_ITS <- transform_sample_counts(glom_ITS, function(x) x / sum(x) )
glom2_ITS # should list taxa or guilds
# Create dataframe from phyloseq object
data_glom_ITS <- psmelt(glom2_ITS)
data_glom_ITS$Phylum <- as.character(data_glom_ITS$Phylum) # 'Phylum' for 16S and ITS and 'Family' for NEM when TAXA # 'GROUP1' or 'GROUP2' for 16S, 'Trophic_Mode' for ITS and 'functional_guild' (or 'cp') for NEM when GUILDS


# Recreate the "practice" and "cultivar" columns lost during sample merging # To draw graph # Don't need if you merge by practice
data_glom_ITS <- data_glom_ITS %>% 
  mutate(Practice2= stringr::str_extract(Sample, "^.{2}"))  %>% 
  mutate(Cultivar2=substring(Sample,3))


data_glom_ITS <- data_glom_ITS %>% 
  mutate(Practices2= stringr::str_extract(Sample, "^.{2}"))  %>% 
  mutate(Compartment2=substring(Sample,3))
write.table(data_glom_ITS, "data_glom.csv", quote = FALSE, sep = ",", row.names = FALSE)

##merge by practices_recreate_practices colum
data_glom_ITS <- data_glom_ITS %>% 
  mutate(Practice2= stringr::str_extract(Sample, "^.{2}"))  

# Rename phyla with < 1% abundance
data_glom2_ITS=data_glom_ITS

data_glom2_ITS$Phylum[data_glom2_ITS$Abundance < 0.01] <- "Other" # 'Phylum' for 16S and ITS, 'Family' for NEM


write.table(data_glom2_ITS, file = "Relabund_family_16S_merged_practice_glom2.csv", sep = "\t", row.names = FALSE)


Count = length(unique(data_glom2_ITS$Phylum)) # 'Phylum' for 16S and ITS, 'Family' for NEM when TAXA # 'GROUP1' or 'GROUP2' for 16S, 'Trophic_Mode' for ITS or 'cp' for NEM when GUILDS
Count

# Draw the figure for relabund by taxa

#Family_merge_by_practices
Relabund_taxa_NEM <- ggplot(data=data_glom2_ITS, aes(x=Abundance, fill=Phylum, y=Practice2)) + geom_bar(aes(), stat="identity", position="stack") + theme(legend.position="bottom") + theme_classic() + theme(axis.title = element_text(color="black", size=9, face="bold")) + theme(axis.text = element_text(color="black", size=7, face="bold")) + theme(legend.text = element_text(color="black", size=8, face="bold")) + theme(legend.title = element_text(color="black", size=8, face="bold")) + theme(strip.background = element_rect(fill = "#F0EEEC"), strip.text = element_text(colour = "black", face = "bold")) + theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype="solid")) + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(strip.text.y = element_text(size=8, angle=0, face = "bold")) + xlab("Relative Abundance") + scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + theme(legend.position = "bottom") + guides(fill=guide_legend(nrow=10,byrow=TRUE)) + theme(panel.spacing = unit(0.1, "lines")) # + scale_fill_manual(values=mycolors)
Relabund_taxa_NEM 

#Genus_merge_by_practices
Relabund_taxa_NEM <- ggplot(data=data_glom2_ITS, aes(x=Abundance, fill=Genus, y=Practice2)) + geom_bar(aes(), stat="identity", position="stack") + theme(legend.position="bottom") + theme_classic() + theme(axis.title = element_text(color="black", size=9, face="bold")) + theme(axis.text = element_text(color="black", size=7, face="bold")) + theme(legend.text = element_text(color="black", size=8, face="bold")) + theme(legend.title = element_text(color="black", size=8, face="bold")) + theme(strip.background = element_rect(fill = "#F0EEEC"), strip.text = element_text(colour = "black", face = "bold")) + theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype="solid")) + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(strip.text.y = element_text(size=8, angle=0, face = "bold")) + xlab("Relative Abundance") + scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + theme(legend.position = "bottom") + guides(fill=guide_legend(nrow=10,byrow=TRUE)) + theme(panel.spacing = unit(0.1, "lines")) # + scale_fill_manual(values=mycolors)
Relabund_taxa_NEM 

#Putative_feeding_merge_by_practices
Relabund_taxa_NEM <- ggplot(data=data_glom2_ITS, aes(x=Abundance, fill=Putative_feeding, y=Practice2)) + geom_bar(aes(), stat="identity", position="stack") + theme(legend.position="bottom") + theme_classic() + theme(axis.title = element_text(color="black", size=9, face="bold")) + theme(axis.text = element_text(color="black", size=7, face="bold")) + theme(legend.text = element_text(color="black", size=8, face="bold")) + theme(legend.title = element_text(color="black", size=8, face="bold")) + theme(strip.background = element_rect(fill = "#F0EEEC"), strip.text = element_text(colour = "black", face = "bold")) + theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype="solid")) + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(strip.text.y = element_text(size=8, angle=0, face = "bold")) + xlab("Relative Abundance") + scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + theme(legend.position = "bottom") + guides(fill=guide_legend(nrow=10,byrow=TRUE)) + theme(panel.spacing = unit(0.1, "lines")) # + scale_fill_manual(values=mycolors)
Relabund_taxa_NEM 

#+ facet_grid(rows=vars(Practice2))

#List <- c("Pratylenchidae", "Meloidogynidae", "Cephalobidae", "Aphelenchoididae", "Anguinidae", "Telotylenchidae", "Criconematidae")

glimpse(data_glom2_ITS)
Relabund_taxa_NEM <- ggplot(data_glom2_ITS, aes(x=Abundance,fill=Phylum , y=Compartment2)) + 
  geom_bar(aes(), stat="identity", position="stack") + theme(legend.position="bottom") + 
  theme_classic() + theme(axis.title = element_text(color="black", size=9, face="bold")) + 
  theme(axis.text = element_text(color="black", size=7, face="bold")) + 
  theme(legend.text = element_text(color="black", size=8, face="bold")) + 
  theme(legend.title = element_text(color="black", size=8, face="bold")) + 
  theme(strip.background = element_rect(fill = "#F0EEEC"), strip.text = element_text(colour = "black", face = "bold")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype="solid")) + 
  theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(strip.text.y = element_text(size=8, angle=0, face = "bold")) + 
  xlab("Relative Abundance") + ylab("Practices") + scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme(legend.position = "bottom") + guides(fill=guide_legend(nrow=10,byrow=TRUE)) + 
  theme(panel.spacing = unit(0.1, "lines"))  + facet_grid(rows=vars(Practices2)) # + scale_fill_manual(values=mycolors)

Relabund_taxa_NEM 

write.table(data_glom2_ITS, "data_glom2_ITS.csv", sep="\t")

##Total abundance 
#Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(ps2))
standf = function(x, t=total) round(t * (x / sum(x)))
ps2 = transform_sample_counts(ps2, standf)


ps2 <- subset_taxa(ps2, Putative_feeding != "NA")


plot_bar(ps2, x="Practices", fill = "Family", facet_grid = ~Cultivar) +
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") 

plot_bar(ps2, x="Practices", fill = "Family") +
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") 

#By genus
plot_bar(ps2, x="Practices", fill = "Genus", facet_grid = ~Cultivar) +
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") 

plot_bar(ps2, x="Practices", fill = "Genus") +
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") 

#by_feeder
plot_bar(ps2, x="Practices", fill = "Putative_feeding", facet_grid = ~Cultivar) +
  geom_bar(aes(color=Putative_feeding, fill=Putative_feeding), stat="identity", position="stack") 

plot_bar(ps2, x="Practices", fill = "Putative_feeding") +
  geom_bar(aes(color=Putative_feeding, fill=Putative_feeding), stat="identity", position="stack") 

# Calculate mean and standard deviation
mean_sd <- aggregate(data_glom2_ITS$Abundance, by=list(sample = data_glom2_ITS$Practice2), FUN = function(x) c(mean=mean(x), sd=sd(x), sum=sum(x), n=length(x))) # cultivar = relabund$cultivar OR practice = relabund$practice
mean_sd <- do.call(data.frame, mean_sd)
write.table(mean_sd, "mean_sd_sum_n_relabund_NEM_practice_roots.csv", quote = FALSE, sep = ",", row.names = FALSE)

# Make statistical tests
## ANOVA 2 factors
stats_df <- data_glom2_ITS %>% filter(Phylum == "Anguinidae") %>%  # start with data
  mutate(Practice2 = factor(Practice2, 
                            levels = c("CA", "CT")),
         Cultivar2 = factor(Cultivar2, 
                            labels = c("IR64", "IR504", "Azucena", "Zhonghua")))



glimpse(stats_df)

ad_aov <- aov(Abundance ~ Practice2, 
              data = stats_df)

ad_aov <- aov(Abundance ~ Cultivar2, 
              data = stats_df)
TukeyHSD(ad_aov)





AOV <- stats_df %>% group_by(Family) %>% aov(Abundance ~ Practice2)

write.table(AOV, "ANOVA_relabund_Family_roots_practices_cultivar_test.csv", quote = FALSE, sep = ",", row.names = FALSE)




AOV <- data_glom2_ITS %>% group_by(Phylum) %>% do(tidy(aov(Abundance ~ Compartment* Practices, .)))

TukeyHSD(ad_aov)



Turkey <- data_glom2_ITS %>% group_by(Family) %>% do(tidy(TukeyHSD(AOV)))

model_CEC <- lme(Abundance ~ Practices * Cultivar, random = ~ 1|replicate, data = meta_nem_and_soil_dt)
ANOVA_CEC <- anova(model_CEC)
CLD_CEC <- cld((emmeans(model_CEC, ~ practice * cultivar, adjust="tukey")), alpha=0.05, Letters=letters)




write.table(AOV, "ANOVA_relabund_phylum_cul_prac.csv", quote = FALSE, sep = ",", row.names = FALSE)
TukeyHSD(AOV)

## EMMEANS
CLD <- data_glom2_ITS %>% group_by(Family) %>% do(tidy(cld(emmeans(glm(Abundance ~ Practices, data = data_glom2_ITS), ~ Practices, adjust="tukey"), alpha=0.05, Letters=letters)))
write.table(CLD, "CLD_relabund_taxa_Family_roots_practice_cp.csv", quote = FALSE, sep = ",", row.names = FALSE)

AOV <- data_glom2_ITS %>% group_by(Family) %>% do(tidy(aov(Abundance ~ Practices, .)))
AOV

write.table(AOV, "ANOVA_relabund_Genusg_roots.csv", quote = FALSE, sep = ",", row.names = FALSE)
#by phylum
AOV <- data_glom2_ITS %>% group_by(Phylum) %>% do(tidy(aov(Abundance ~ Compartment* Practices, .)))
write.table(AOV, "ANOVA_relabund_Compartment_Practices.csv", quote = FALSE, sep = ",", row.names = FALSE)

##by genus
AOV <- data_glom2_ITS %>% group_by(Genus) %>% do(tidy(aov(Abundance ~ Practices, .)))
write.table(AOV, "ANOVA_relabund_Genusg_roots.csv", quote = FALSE, sep = ",", row.names = FALSE)

AOV <- data_glom2_ITS %>% group_by(Genus) %>% do(tidy(aov(Abundance ~ Practices* Cultivar, .)))
write.table(AOV, "ANOVA_relabund_genus_roots_pac_cul.csv", quote = FALSE, sep = ",", row.names = FALSE)


##by putative feeding
AOV <- data_glom2_ITS %>% group_by(Putative_feeding) %>% do(tidy(aov(Abundance ~ Practices, .)))
write.table(AOV, "ANOVA_relabund_Putative_feeding_roots_wo_NA.csv", quote = FALSE, sep = ",", row.names = FALSE)

AOV <- data_glom2_ITS %>% group_by(Putative_feeding) %>% do(tidy(aov(Abundance ~ Practices* Cultivar, .)))
write.table(AOV, "ANOVA_relabund_Putative_feeding_roots_pac_cul_wo_NA.csv", quote = FALSE, sep = ",", row.names = FALSE)

###OTHER GRAPH
# Filter phyloseq object
ps_3_NEM <- filter_taxa(df2, function(x) sum(x) > 0, TRUE) # Filter SVs in the taxo
check_taxo=ps_3_NEM@tax_table@.Data # Give the number of SVs based on above criteria
# Export phyloseq object
#saveRDS(ps_3_NEM, "physeq_NEM_roots.rds")
# Extract tables from the phyloseq object
OTU_NEM_ps3 = as(otu_table(ps_3_NEM), "matrix")
OTU_NEM_ps3_df = as.data.frame(OTU_NEM_ps3)
#write.csv(OTU_NEM_ps3_df,"OTU_NEM_onlyrhizo_woD-Z-R3_nemaplex.csv")

TAXO_NEM_ps3 = as(tax_table(ps_3_NEM), "matrix")
TAXO_NEM_ps3_df = as.data.frame(TAXO_NEM_ps3)
#write.csv(TAXO_NEM_ps3_df,"TAXO_NEM_onlyrhizo_woD-Z-R3_nemaplex.csv")

META_NEM = as(sample_data(ps_3_NEM), "matrix")
META_NEM_ps3_df = as.data.frame(META_NEM)
#write.csv(META_NEM_ps3_df,"META_NEM_onlyrhizo_woD-Z-R3_nemaplex.csv")
# Create SV_use (from OTU_ps + META_ps) and matrix_use for the analysis
SV_use_NEM_ps3=merge(META_NEM_ps3_df, t(OTU_NEM_ps3_df), by="row.names")
colnames(SV_use_NEM_ps3)[1] <- "ID"
ncol(SV_use_NEM_ps3)
write.csv(SV_use_NEM_ps3,"META_NEM_roots.csv")
matrix_NEM_ps3<-SV_use_NEM_ps3[c(7:2345)] # (number of metadata variables +1 : total number of variables in SV_use) #58 rhizophere, 108 roots
matrix_NEM_ps3<-mapply(matrix_NEM_ps3, FUN=as.numeric)
matrix_use_NEM_ps3<-matrix_NEM_ps3[,colSums(matrix_NEM_ps3)>=1]
#5.B. Community description
#For each community, replace "16S", "ITS" or "NEM" and run the whole part.

#5.B.Global (Venn diagrams)
# Make the lists of sample per practice
list_Root<-subset(META_NEM_ps3_df, Compartment=='Roots')

list_Rhizosphere<-subset(META_NEM_ps3_df, Compartment=='Rhizosphere')



list_Root_CA<-subset(list_Root, Practices=='CA')
list_Root_CT<-subset(list_Root, Practices=='CT')
list_Rhizosphere_CA<-subset(list_Rhizosphere, Practices=='CA')
list_Rhizosphere_CT<-subset(list_Rhizosphere, Practices=='CT')
  
  
list_Root<-rownames(list_Root)
list_Root_CA<-rownames(list_Root_CA)
list_Root_CT<-rownames(list_Root_CT)
list_Rhizosphere<-rownames(list_Rhizosphere)
list_Rhizosphere_CA<-rownames(list_Rhizosphere_CA)
list_Rhizosphere_CT<-rownames(list_Rhizosphere_CT)

list_CA<-subset(META_NEM_ps3_df, Practices=='CA')
list_CA<-rownames(list_CA)
list_CT<-subset(META_NEM_ps3_df, Practices=='CT')
list_CT<-rownames(list_CT)

# Create column to resume abundance in binary (0=absent, 1=present) to draw the diagram
OTU_NEM_ps3_df$CT <- ifelse(rowSums(OTU_NEM_ps3_df[,list_CT], na.rm=T)>0,1,0)
OTU_NEM_ps3_df$CA <- ifelse(rowSums(OTU_NEM_ps3_df[,list_CA], na.rm=T)>0,1,0)

OTU_NEM_ps3_df$Roots <- ifelse(rowSums(OTU_NEM_ps3_df[,list_Root], na.rm=T)>0,1,0)
OTU_NEM_ps3_df$Rhizosphere <- ifelse(rowSums(OTU_NEM_ps3_df[,list_Rhizosphere], na.rm=T)>0,1,0)

OTU_NEM_ps3_df$Roots_CA <- ifelse(rowSums(OTU_NEM_ps3_df[,list_Root_CA], na.rm=T)>0,1,0)
OTU_NEM_ps3_df$Roots_CT <- ifelse(rowSums(OTU_NEM_ps3_df[,list_Root_CT], na.rm=T)>0,1,0)
OTU_NEM_ps3_df$Rhizosphere_CA <- ifelse(rowSums(OTU_NEM_ps3_df[,list_Rhizosphere_CA], na.rm=T)>0,1,0)
OTU_NEM_ps3_df$Rhizosphere_CT <- ifelse(rowSums(OTU_NEM_ps3_df[,list_Rhizosphere_CT], na.rm=T)>0,1,0)

# Make the lists of sample per cultivar
list_IR504<-subset(META_NEM_ps3_df, Cultivar=='IR504')
list_IR504<-rownames(list_IR504)
list_IR64<-subset(META_NEM_ps3_df, Cultivar=='IR64')
list_IR64<-rownames(list_IR64)
list_Azucena<-subset(META_NEM_ps3_df, Cultivar=='Azucena')
list_Azucena<-rownames(list_Azucena)
list_Zhonghua<-subset(META_NEM_ps3_df, Cultivar=='Zhonghua')
list_Zhonghua<-rownames(list_Zhonghua)

# Create column to resume abundance in binary (0=absent, 1=present) to draw the diagram
OTU_NEM_ps3_df$IR504 <- ifelse(rowSums(OTU_NEM_ps3_df[,list_IR504], na.rm=T)>0,1,0)
OTU_NEM_ps3_df$IR64 <- ifelse(rowSums(OTU_NEM_ps3_df[,list_IR64], na.rm=T)>0,1,0)
OTU_NEM_ps3_df$Azucena <- ifelse(rowSums(OTU_NEM_ps3_df[,list_Azucena], na.rm=T)>0,1,0)
OTU_NEM_ps3_df$Zhonghua <- ifelse(rowSums(OTU_NEM_ps3_df[,list_Zhonghua], na.rm=T)>0,1,0)
set.seed=100

# Draw the diagrams

## Compartments
a<-euler(OTU_NEM_ps3_df[,c('Roots_CA','Rhizosphere_CA','Roots_CT','Rhizosphere_CT')],shape="circle", control=list(extraopt_threshold=0,extraopt=T))
venn_practice_NEM <-plot(a,quantities = list(cex = .8),labels = list(cex = 1), fills = list(fill = c("cornflowerblue","lightblue2","khaki3","gold"), alpha = 0.6))
venn_practice_NEM


## Practice
a<-euler(OTU_NEM_ps3_df[,c('CA','CT')],shape="circle", control=list(extraopt_threshold=0,extraopt=T))
venn_practice_NEM <-plot(a,quantities = list(cex = .8),labels = list(cex = 1), fills = list(fill = c("#1B9E77", "#D95F02"), alpha = 0.6))
venn_practice_NEM

## cultivar
a<-euler(OTU_NEM_ps3_df[,c('IR504', 'IR64', 'Azucena', 'Zhonghua')],shape="ellipse", control=list(extraopt_threshold=0,extraopt=T))
venn_cultivar_NEM <-plot(a,quantities = list(cex = .8),labels = list(cex = 1), fills = list(fill = c("cornflowerblue","lightblue2","khaki3","gold"), alpha = 0.6))
venn_cultivar_NEM


##5.B.Structure
# Characterize abundance table
totreads <- sum(matrix_use_NEM_ps3)

totSVs <- ncol(matrix_use_NEM_ps3)

raremin <- min(rowSums(matrix_use_NEM_ps3))
raremax <- max(rowSums(matrix_use_NEM_ps3))
raremedian <- median(rowSums(matrix_use_NEM_ps3))
totreads
totSVs
raremin
raremax
raremedian
#write.csv(matrix_use_NEM_ps3,"matrix_use_NEM_ps3_rhizo.csv")
# Check rarefaction curve
rar=rarecurve(matrix_use_NEM_ps3, step=20, cex = 0.6, label=TRUE, col="blue", xlab="Sample size", ylab="Taxa count") # ylab="SVs count" for 16S and ITS or ylab="Taxa count" for NEM
# Prepare and draw NMDS graph for structure
NMDS_NEM <- metaMDS(matrix_use_NEM_ps3, distance = "bray", trymax = 100)
stressplot(NMDS_NEM)
NMDS_NEM

NMDSsites=scores(NMDS_NEM, display="sites")
SV_use_NEM_NMDS=cbind(SV_use_NEM_ps3,NMDSsites)


SV_use_NEM_NMDS$Compartment<-ordered(SV_use_NEM_NMDS$Compartment, levels=c("Roots", "Rhizosphere"))

SV_use_NEM_NMDS$Practices<-ordered(SV_use_NEM_NMDS$Practices, levels=c("CT", "CA"))
SV_use_NEM_NMDS$Cultivar<-ordered(SV_use_NEM_NMDS$Cultivar, levels=c("IR504", "IR64", "Azucena", "Zhonghua"))
write.csv(SV_use_NEM_NMDS,"SV_use_NEM_NMDS.csv")

# Colors and shapes for ggplot2
mycolors <- scale_color_manual(values = c("Roots"="blue1", "Rhizosphere"="coral4")) # "IR64"="lightblue2", "IR504"="cornflowerblue", "Azucena"="khaki3", "Zhonghua"="gold"))

mycolors <- scale_color_manual(values = c("CA"="#1B9E77", "CT"="#D95F02"))
myshapes <- scale_shape_manual(values = c("IR64"=15, "IR504"=16, "Azucena"=17, "Zhonghua"=18)) #"CA"=3, "CT"=5,

NMDS_NEM_plot = ggplot(data=SV_use_NEM_NMDS, aes(NMDS1, NMDS2, color=Practices, shape=Cultivar)) + geom_point(size=4.5) + theme_bw() + xlab("NMDS1") + ylab("NMDS2") + mycolors + myshapes + theme(axis.text.y = element_text (), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
NMDS_NEM_plot

NMDS_NEM_plot = ggplot(data=SV_use_NEM_NMDS, aes(NMDS1, NMDS2, color=Compartment, shape=Cultivar)) + geom_point(size=4.5) + theme_bw() + xlab("NMDS1") + ylab("NMDS2") + mycolors + myshapes + theme(axis.text.y = element_text (), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
NMDS_NEM_plot


xlim(-316,-314.5) + ylim(-0.5,0.75)
# Export NMDS coordinates
coordinates=cbind(SV_use_NEM_ps3[,c(1:6)],NMDSsites)
write.csv(coordinates,"NMDScoordinates_NEM.csv")
# Make betadisper test on practice effect to test multivariate homogeneity of variance before adonis
dis <- vegdist(matrix_use_NEM_ps3)
Practices=SV_use_NEM_ps3$Practices
Cultivar=SV_use_NEM_ps3$Cultivar

mod_practice <- betadisper(dis, Practices)
permutest(mod_practice)

mod_cultivar <- betadisper(dis, Cultivar)
permutest(mod_cultivar)
# Make adonis: Permutational Multivariate Analysis of Variance Using Distance Matrices (Vegan)
Alltreatments=SV_use_NEM_ps3[c(3:5)]



adonis <- adonis(matrix_use_NEM_ps3 ~ Practices*Cultivar, random = ~ 1|Replicate, data=Alltreatments)
adonis


##5.B.Diversity
# Calculate diversity measures : Richness (SVs/taxa richness) and Shannon (Shannon index)
Richness <- specnumber(matrix_use_NEM_ps3) 
Shannon <- diversity(matrix_use_NEM_ps3)
divtable= cbind(Richness,Shannon)

SV_use_div_NEM=cbind(SV_use_NEM_ps3,divtable)
#SV_use_div_NEM <- SV_use_div_NEM %>% filter(!Practices %in% "NU")
# Extract indices
indices = cbind(SV_use_NEM_ps3[,c(1:6)],divtable)

write.csv(indices,"Diversity_NEM.csv")
# Calculate mean and sd
Icompartment_Richness_NEM=ddply(indices, "Compartment", summarise, mean=mean(Richness), sd=sd(Richness))
Icompartment_Shannon_NEM=ddply(indices, "Compartment", summarise, mean=mean(Shannon), sd=sd(Shannon))

Ipractice_Richness_NEM=ddply(indices, "Practices", summarise, mean=mean(Richness), sd=sd(Richness))
Ipractice_Shannon_NEM=ddply(indices, "Practices", summarise, mean=mean(Shannon), sd=sd(Shannon))
Icultivar_Richness_NEM=ddply(indices, "Cultivar", summarise, mean=mean(Richness), sd=sd(Richness))
Icultivar_Shannon_NEM=ddply(indices, "Cultivar", summarise, mean=mean(Shannon), sd=sd(Shannon))

# Export the results of indices
write.csv(Ipractice_Richness_NEM, "Ipractice_Richness_NEM.csv")
write.csv(Ipractice_Shannon_NEM, "Ipractice_Shannon_NEM.csv")
write.csv(Icultivar_Richness_NEM, "Icultivar_Richness_NEM.csv")
write.csv(Icultivar_Shannon_NEM, "Icultivar_Shannon_NEM.csv")


#Richness -> Poisson (count - ex=richness)
hist(Richness)
model_Richness <- glmer(Richness ~ Compartment*Practices + (1|Replicate), family="poisson", data=SV_use_div_NEM)
model_Richness <- glmer(Richness ~ Practices*Cultivar + (1|Replicate), family="poisson", data=SV_use_div_NEM)

model_Richness <- glmer(Richness ~ Compartment + (1|Replicate), family="poisson", data=SV_use_div_NEM)
model_Richness <- glmer(Richness ~ Practices + (1|Replicate), family="poisson", data=SV_use_div_NEM)
summary(model_Richness)


Anova(model_Richness)
plot(model_Richness)
E <- resid(model_Richness)
hist(E)
qqnorm(E)

CLDem <- cld((emmeans(model_Richness, ~ Practices*Cultivar, adjust="tukey")), alpha=0.05, Letters=letters)
CLDem <- cld((emmeans(model_Richness, ~ Compartment*Practices, adjust="tukey")), alpha=0.05, Letters=letters)

CLDem <- cld((emmeans(model_Richness, ~ Practices, adjust="tukey")), alpha=0.05, Letters=letters)

CLDem
# For Shannon -> Gaussian (identity)
hist(Shannon)


model_Shannon <- glmer(Shannon ~ Practices*Cultivar + (1|Replicate), family=gaussian(link="identity"), data=SV_use_div_NEM)
model_Shannon <- glmer(Shannon ~ Compartment*Practices + (1|Replicate), family=gaussian(link="identity"), data=SV_use_div_NEM)

model_Shannon <- glmer(Shannon ~ Practices + (1|Replicate), family=gaussian(link="identity"), data=SV_use_div_NEM)

summary(model_Shannon)
Anova(model_Shannon)
plot(model_Shannon)
E <- resid(model_Shannon)
hist(E)
qqnorm(E)

CLDem <- cld((emmeans(model_Shannon, ~ Practices*Cultivar, adjust="tukey")), alpha=0.05, Letters=letters)
CLDem <- cld((emmeans(model_Shannon, ~ Compartment*Practices, adjust="tukey")), alpha=0.05, Letters=letters)
CLDem <- cld((emmeans(model_Shannon, ~ Practices, adjust="tukey")), alpha=0.05, Letters=letters)

CLDem
# Draw the graphs of diversity
SV_use_div_NEM$Practices<-ordered(SV_use_div_NEM$Practices, levels=c("CT", "CA"))
SV_use_div_NEM$Compartment<-ordered(SV_use_div_NEM$Compartment, levels=c("Roots", "Rhizosphere"))
## Richness
### cultivar
Richness_cultivar_NEM <- ggplot(data=SV_use_div_NEM, aes(x=Cultivar, y=Richness, fill=Practices), size=6) + geom_boxplot(outlier.shape=NA) + theme_bw() + xlab("Rice cultivar") + ylab("Observed richness") + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + scale_fill_manual(values = c("CA"="#1B9E77", "CT"="#D95F02")) + scale_x_discrete(limits=c("IR504", "IR64", "Azucena", "Zhonghua"))
Richness_cultivar_NEM

### practice
Richness_practice_NEM <- ggplot(data=SV_use_div_NEM, aes(x=Practices, y=Richness, fill=Practices), size=6) + geom_boxplot(outlier.shape=NA) + theme_bw() + xlab("Practice") + ylab("Observed richness") + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + scale_fill_manual(values = c("CA"="#1B9E77", "CT"="#D95F02"))
Richness_practice_NEM

### compartment
Richness_practice_NEM <- ggplot(data=SV_use_div_NEM, aes(x=Compartment, y=Richness, fill=Compartment), size=6) + geom_boxplot(outlier.shape=NA) + theme_bw() + xlab("Compartment") + ylab("Observed richness") + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + scale_fill_manual(values = c("Roots"="cyan", "Rhizosphere"="darkviolet")) 
Richness_practice_NEM

Richness_practice_NEM <- ggplot(data=SV_use_div_NEM, aes(x=Practices, y=Richness, fill=Compartment), size=6) + 
  geom_boxplot(outlier.shape=NA) + theme_bw() + xlab("Practice") + ylab("Observed richness") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + 
  scale_fill_manual(values = c("Roots"="cyan", "Rhizosphere"="darkviolet")) + scale_x_discrete(limits=c("CA", "CT"))
Richness_practice_NEM


Richness_practice_NEM <- ggplot(data=SV_use_div_NEM, aes(x=Compartment, y=Richness, fill=Practices), size=6) + 
  geom_boxplot(outlier.shape=NA) + theme_bw() + xlab("Compartment") + ylab("Observed richness") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + 
  scale_fill_manual(values = c("CA"="#1B9E77", "CT"="#D95F02")) + scale_x_discrete(limits=c("Roots", "Rhizosphere"))
Richness_practice_NEM


Richness_practice_NEM <- ggplot(data=SV_use_div_NEM, aes(x=Compartment, y=Richness, fill=Practices), size=6) + 
  geom_boxplot(outlier.shape=NA) + theme_bw() + xlab("Compartment") + ylab("Observed richness") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + 
  scale_fill_manual(values = c("CA"="#1B9E77", "CT"="#D95F02")) + scale_x_discrete(limits=c("Roots", "Rhizosphere"))
Richness_practice_NEM

## Shannon
### cultivar
Shannon_cultivar_NEM <- ggplot(data=SV_use_div_NEM, aes(x=Cultivar, y=Shannon, fill=Practices), size=6) + geom_boxplot(outlier.shape=NA) + theme_bw() + xlab("Rice cultivar") + ylab("Shannon index") + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + scale_fill_manual(values = c("CA"="#1B9E77", "CT"="#D95F02")) + scale_x_discrete(limits=c("IR504", "IR64", "Azucena", "Zhonghua"))
Shannon_cultivar_NEM

### practice
Shannon_practice_NEM <- ggplot(data=SV_use_div_NEM, aes(x=Practices, y=Shannon, fill=Practices), size=6) + 
  geom_boxplot(outlier.shape=NA) + theme_bw() + xlab("Practice") + ylab("Shannon index") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + 
  scale_fill_manual(values = c("CA"="#1B9E77", "CT"="#D95F02"))
Shannon_practice_NEM

Shannon_practice_NEM <- ggplot(data=SV_use_div_NEM, aes(x=Compartment, y=Shannon, fill=Compartment), size=6) + 
  geom_boxplot(outlier.shape=NA) + theme_bw() + xlab("Compartment") + ylab("Shannon index") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + 
  scale_fill_manual(values = c("Roots"="cyan", "Rhizosphere"="darkviolet"))
Shannon_practice_NEM

Shannon_practice_NEM<- ggplot(data=SV_use_div_NEM, aes(x=Practices, y=Shannon, fill=Compartment), size=6) + 
  geom_boxplot(outlier.shape=NA) + theme_bw() + xlab("Practice") + ylab("Shannon index") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + 
  scale_fill_manual(values = c("Roots"="cyan", "Rhizosphere"="darkviolet")) + scale_x_discrete(limits=c("CA", "CT"))
Shannon_practice_NEM

Shannon_practice_NEM <- ggplot(data=SV_use_div_NEM, aes(x=Compartment, y=Shannon, fill=Practices), size=6) + 
  geom_boxplot(outlier.shape=NA) + theme_bw() + xlab("Compartment") + ylab("Shannon index") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + 
  scale_fill_manual(values = c("CA"="#1B9E77", "CT"="#D95F02")) + scale_x_discrete(limits=c("Roots", "Rhizosphere"))
Shannon_practice_NEM
# Export the graphs of the three communities in one figure
Richness_cultivar_16S_f = Richness_cultivar_16S + ylab(NULL) + xlab(NULL)
Richness_practice_16S_f = Richness_practice_16S + ylab(NULL) + xlab(NULL)
Shannon_cultivar_16S_f = Shannon_cultivar_16S + ylab(NULL) + xlab(NULL)
Shannon_practice_16S_f = Shannon_practice_16S + ylab(NULL) + xlab(NULL)
Richness_cultivar_ITS_f = Richness_cultivar_ITS + ylab(NULL) + xlab(NULL)
Richness_practice_ITS_f = Richness_practice_ITS + ylab(NULL) + xlab(NULL)
Shannon_cultivar_ITS_f = Shannon_cultivar_ITS + ylab(NULL) + xlab(NULL)
Shannon_practice_ITS_f = Shannon_practice_ITS + ylab(NULL) + xlab(NULL)
Richness_cultivar_NEM_f = Richness_cultivar_NEM + ylab(NULL) + xlab(NULL)
Richness_practice_NEM_f = Richness_practice_NEM + ylab(NULL) + xlab(NULL)
Shannon_cultivar_NEM_f = Shannon_cultivar_NEM + ylab(NULL) + xlab(NULL)
Shannon_practice_NEM_f = Shannon_practice_NEM + ylab(NULL) + xlab(NULL)

cultivar = plot_grid(Richness_cultivar_16S_f, Richness_cultivar_ITS_f, Richness_cultivar_NEM_f, Shannon_cultivar_16S_f, Shannon_cultivar_ITS_f, Shannon_cultivar_NEM_f, align = "hv", labels=c('A', 'C', 'E', 'G', 'I', 'K'), hjust=-3, nrow=2, ncol=3)

practice = plot_grid(Richness_practice_16S_f, Richness_practice_ITS_f, Richness_practice_NEM_f, Shannon_practice_16S_f, Shannon_practice_ITS_f, Shannon_practice_NEM_f, align = "hv", labels=c('B', 'D', 'F', 'H', 'J', 'L'), hjust = -3, nrow=2, ncol=3)

save_plot("Figure_X_diversity_cultivar.svg", cultivar, base_width=5, base_height=2.5, nrow=2, ncol=3)
save_plot("Figure_X_diversity_practice.svg", practice, base_width=2.5, base_height=2.5, nrow=2, ncol=3)

save_plot("Figure_X_diversity_richnessaxis.svg", Richness_cultivar_ITS, base_width=5, base_height=2.5)
save_plot("Figure_X_diversity_shannonaxis.svg", Shannon_cultivar_ITS, base_width=5, base_height=2.5)

##
