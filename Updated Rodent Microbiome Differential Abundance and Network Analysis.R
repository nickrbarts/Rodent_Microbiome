rm(list=ls())
library(SpiecEasi)
library(NetCoMi)
library(phyloseq)
library(qgraph)
library(igraph)
library(vegan)
library(MCL)
library(qiime2R)
library(ggplot2)
library(phyloseqCompanion)
library(CoDaSeq)
library(microbiome)
library(DESeq2)
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
library(scales)
library(dplyr)
library(grid)
library(reshape2)
library(ape)
library(ashr)
library(limma)
library(edgeR)
library(splines)
library(RColorBrewer)

setwd("~/Desktop/RodentMicrobiome2024")
###Herbivore - Vole######
#####create your phyloseq object##########
##16S##
physeq_16S_MOVO <- qza_to_phyloseq(features="table-bacteria-MOVO.qza", taxonomy="taxonomy.qza", metadata="rodent-metadata-ut.txt")
rank_names(physeq_16S_MOVO)
table(tax_table(physeq_16S_MOVO)[, "Genus"], exclude = NULL) #view low abundant or uncharacterized taxa
physeq_16S_MOVO_sub <- subset_taxa(physeq_16S_MOVO, !is.na(Genus) & !Genus %in% c("", "uncharacterized")) #remove any with Genus = NA
# Compute prevalence of each feature, store as data.frame
prevdf_MOVO = apply(X = otu_table(physeq_16S_MOVO_sub),
                    MARGIN = ifelse(taxa_are_rows(physeq_16S_MOVO_sub), yes = 1, no = 2),
                    FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_MOVO = data.frame(Prevalence = prevdf_MOVO,
                         TotalAbundance = taxa_sums(physeq_16S_MOVO_sub),
                         tax_table(physeq_16S_MOVO_sub))
plyr::ddply(prevdf_MOVO, "Genus", function(df){cbind(mean(df$Prevalence),sum(df$Prevalence))}) #Shows prevalence of phyla, remove low prevalnence(?)
prevdf1_MOVO = subset(prevdf_MOVO, Genus %in% get_taxa_unique(physeq_16S_MOVO_sub, "Genus"))
ggplot(prevdf1_MOVO, aes(TotalAbundance, Prevalence / nsamples(physeq_16S_MOVO_sub),color=Genus)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Genus) + theme(legend.position="none")
prevalenceThreshold = 0.05* nsamples(physeq_16S_MOVO_sub)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1_MOVO)[(prevdf1_MOVO$Prevalence >= prevalenceThreshold)]
ps_16S_MOVO = prune_taxa(keepTaxa, physeq_16S_MOVO_sub)
# How many genera would be present after filtering?
length(get_taxa_unique(ps_16S_MOVO, taxonomic.rank = "Genus"))


##ITS##
physeq_ITS_MOVO <- qza_to_phyloseq(features="~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/table-sans-no-wild-fecals-MOVO-only.qza", taxonomy="~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/taxonomy.qza", metadata="~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/rodent-metadata-ut-full.txt")
rank_names(physeq_ITS_MOVO)
table(tax_table(physeq_ITS_MOVO)[, "Genus"], exclude = NULL) #view low abundant or uncharacterized taxa
physeq_ITS_MOVO_sub <- subset_taxa(physeq_ITS_MOVO, !is.na(Genus) & !Genus %in% c("", "uncharacterized")) #remove any with Phyla = NA
# Compute prevalence of each feature, store as data.frame
prevdf_ITS_MOVO = apply(X = otu_table(physeq_ITS_MOVO_sub),
                        MARGIN = ifelse(taxa_are_rows(physeq_ITS_MOVO_sub), yes = 1, no = 2),
                        FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_ITS_MOVO = data.frame(Prevalence = prevdf_ITS_MOVO,
                             TotalAbundance = taxa_sums(physeq_ITS_MOVO_sub),
                             tax_table(physeq_ITS_MOVO_sub))
plyr::ddply(prevdf_ITS_MOVO, "Genus", function(df){cbind(mean(df$Prevalence),sum(df$Prevalence))}) #Shows prevalence of phyla, remove low prevalnence(?)
prevdf1_ITS_MOVO = subset(prevdf_ITS_MOVO, Genus %in% get_taxa_unique(physeq_ITS_MOVO_sub, "Genus"))
ggplot(prevdf1_ITS_MOVO, aes(TotalAbundance, Prevalence / nsamples(physeq_ITS_MOVO_sub),color=Genus)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
prevalenceThreshold = 0.05* nsamples(physeq_ITS_MOVO_sub)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1_ITS_MOVO)[(prevdf1_ITS_MOVO$Prevalence >= prevalenceThreshold)]
ps_ITS_MOVO = prune_taxa(keepTaxa, physeq_ITS_MOVO_sub)
# How many genera would be present after filtering?
length(get_taxa_unique(ps_ITS_MOVO, taxonomic.rank = "Genus"))


#Merge#
ps_ITS_16S_MOVO <- merge_phyloseq(ps_16S_MOVO, ps_ITS_MOVO)
ps_ITS_16S_MOVO_genus <- aggregate_taxa(ps_ITS_16S_MOVO, "Genus", verbose=FALSE)
otu <- otu_table(ps_ITS_16S_MOVO_genus)
otu <- as.data.frame(otu)
otu <- t(otu)
metadata <- read.table("rodent-metadata-ut.txt", header=TRUE)
otu <- cbind(otu, metadata[, 4:6])
write.csv(otu, "~/Desktop/MOVO_OTU_Table.csv")


###tax_glom to get to genus level###
ps_ITS_16S_MOVO_genus <- tax_glom(ps_ITS_16S_MOVO, taxrank="Genus", NArm=FALSE)
ps_ITS_16S_MOVO_genus_OTU <- as(otu_table(ps_ITS_16S_MOVO_genus), "matrix")


######DESeq2 Differential Abundance######
#DESeq Analysis
dds = phyloseq_to_deseq2(ps_ITS_16S_MOVO_genus, ~ Fiber*Protein)  ####NOTE: CANNOT USE ps_g_16S_MOVO due to an error, but can use phyloseq object####
dds$Fiber <- relevel(dds$Fiber, "Low")
dds$Protein <- relevel(dds$Protein, "Low")

#Calculate geometric means
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, type="poscounts")

#Run DESeq2
dds = DESeq(dds, fitType="local")
resultsNames(dds) # "Intercept" "Fiber_High_vs_Low" "Protein_High_vs_Low" "FiberHigh.ProteinHigh"
res = results(dds, name="FiberHigh.ProteinHigh")
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_ITS_16S_MOVO_genus)[rownames(sigtab), ], "matrix"))
head(sigtab)
summary(sigtab)
summary(res)

#Shrink
shrink <- lfcShrink(dds, coef=c("FiberHigh.ProteinHigh"), sigtab, type="ashr", lfcThreshold = 0)
shrink <- shrink[order(shrink$padj, na.last=NA), ]
sigtab1 = shrink[(shrink$padj < alpha), ]
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(ps_ITS_16S_MOVO_genus)[rownames(sigtab1), ], "matrix"))
summary(sigtab1)
write.csv(sigtab1, "MOVO DESeq Interaction.csv")

DESeq2::plotMA(shrink, ylim=c(-25,25)) #MAplot to look at 
shrink1 = as.data.frame(shrink)
shrink1 = mutate(shrink1, sig=ifelse(shrink1$padj<0.1, "FDR<0.1", "NotSig"))
shrink1[which(abs(shrink1$log2FoldChange)<1.0), "sig"] = "NotSig"
shrink1 = cbind(as(shrink1, "data.frame"), as(tax_table(ps_ITS_16S_MOVO_genus)[rownames(shrink1), ], "matrix"))


write.csv(shrink1, "MOVO_DEseq2_Interaction_Genus.csv")


####Taxon Filtering####
#get 16S genus#
# Agglomerate to genus level
ps_g_16S_MOVO <- phyloseq::tax_glom(ps_16S_MOVO, taxrank = "Genus")
taxtab <- ps_g_16S_MOVO@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Family)
miss_f <- which(taxtab[, "Family"] == "f__")
miss_g <- which(taxtab[, "Genus"] == "g__")

# Number unspecified genera
taxtab[miss_f, "Family"] <- paste0("f__", 1:length(miss_f))
taxtab[miss_g, "Genus"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g <- which(duplicated(taxtab[, "Genus"]) |
                  duplicated(taxtab[, "Genus"], fromLast = TRUE))

for(i in seq_along(taxtab)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Order"], ")")
  } else if(i %in% miss_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Family"], ")")
  }
  
  # Family names are added to duplicate genera
  if(i %in% dupl_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Family"], ")")
  }
}

ps_g_16S_MOVO@tax_table@.Data <- taxtab
rownames(ps_g_16S_MOVO@otu_table@.Data) <- taxtab[, "Genus"]

#get ITS genus#
ps_genus_ITS_MOVO <- phyloseq::tax_glom(ps_ITS_MOVO, taxrank = "Genus")
taxtab_ITS_MOVO <- ps_genus_ITS_MOVO@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Family)
miss_f <- which(taxtab_ITS_MOVO[, "Family"] == "f__")
miss_g <- which(taxtab_ITS_MOVO[, "Genus"] == "g__")

# Number unspecified genera
taxtab_ITS_MOVO[miss_f, "Family"] <- paste0("f__", 1:length(miss_f))
taxtab_ITS_MOVO[miss_g, "Genus"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g1 <- which(duplicated(taxtab_ITS_MOVO[, "Genus"]) |
                   duplicated(taxtab_ITS_MOVO[, "Genus"], fromLast = TRUE))

for(i in seq_along(taxtab_ITS_MOVO)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab_ITS_MOVO[i, "Genus"] <- paste0(taxtab_ITS_MOVO[i, "Genus"], "(", taxtab_ITS_MOVO[i, "Order"], ")")
  } else if(i %in% miss_g){
    taxtab_ITS_MOVO[i, "Genus"] <- paste0(taxtab_ITS_MOVO[i, "Genus"], "(", taxtab_ITS_MOVO[i, "Family"], ")")
  }
  
  # Family names are added to duplicate genera
  if(i %in% dupl_g1){
    taxtab_ITS_MOVO[i, "Genus"] <- paste0(taxtab_ITS_MOVO[i, "Genus"], "(", taxtab_ITS_MOVO[i, "Family"], ")")
  }
}

ps_genus_ITS_MOVO@tax_table@.Data <- taxtab_ITS_MOVO
rownames(ps_genus_ITS_MOVO@otu_table@.Data) <- taxtab_ITS_MOVO[, "Genus"]


####merged######
ps_ITS_16S_MOVO_g <- phyloseq::tax_glom(ps_ITS_16S_MOVO, taxrank = "Genus")
taxtab_ITS_16S_MOVO_g <- ps_ITS_16S_MOVO_g@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Family)
miss_f <- which(taxtab_ITS_16S_MOVO_g[, "Family"] == "f__")
miss_g <- which(taxtab_ITS_16S_MOVO_g[, "Genus"] == "g__")

# Number unspecified genera
taxtab_ITS_16S_MOVO_g[miss_f, "Family"] <- paste0("f__", 1:length(miss_f))
taxtab_ITS_16S_MOVO_g[miss_g, "Genus"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g1 <- which(duplicated(taxtab_ITS_16S_MOVO_g[, "Genus"]) |
                   duplicated(taxtab_ITS_16S_MOVO_g[, "Genus"], fromLast = TRUE))

for(i in seq_along(taxtab_ITS_16S_MOVO_g)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab_ITS_16S_MOVO_g[i, "Genus"] <- paste0(taxtab_ITS_16S_MOVO_g[i, "Genus"], "(", taxtab_ITS_16S_MOVO_g[i, "Order"], ")")
  } else if(i %in% miss_g){
    taxtab_ITS_16S_MOVO_g[i, "Genus"] <- paste0(taxtab_ITS_16S_MOVO_g[i, "Genus"], "(", taxtab_ITS_16S_MOVO_g[i, "Family"], ")")
  }
  
  # Family names are added to duplicate genera
  if(i %in% dupl_g1){
    taxtab_ITS_16S_MOVO_g[i, "Genus"] <- paste0(taxtab_ITS_16S_MOVO_g[i, "Genus"], "(", taxtab_ITS_16S_MOVO_g[i, "Family"], ")")
  }
}

ps_ITS_16S_MOVO_g@tax_table@.Data <- taxtab_ITS_16S_MOVO_g
rownames(ps_ITS_16S_MOVO_g@otu_table@.Data) <- taxtab_ITS_16S_MOVO_g[, "Genus"]
ITS_16S_MOVO_genus_table <- as(otu_table(ps_ITS_16S_MOVO_g), "matrix")


######glmmTMB########
library(broom)
library(car)
library(magicfor)
library(glmmTMB)
library(AICcmodavg)
library(bbmle)

data <- ITS_16S_MOVO_genus_table
df_transpose <- as.data.frame(t(as.matrix(data)))
metadata <- read.table("rodent-metadata-ut.txt", header=TRUE)
df_metadata <- cbind(df_transpose, metadata$Protein)
df_metadata <- cbind(df_metadata, metadata$Fiber)
df_metadata <- write.csv(df_metadata, "df_metadata.csv")
df_metadata <- read.csv("df_metadata.csv")
fiber <- df_metadata$Fiber
protein <- df_metadata$Protein


# names of lipids
microbe.names <- colnames(df_metadata)[4:120]
no.microbes <- length(microbe.names)

# create a named list to hold the fitted models
fitlistnb2 <- as.list(1:no.microbes)
names(fitlistnb2) <- microbe.names
review_nb2 <- vector("list", length=117)

# loop over lipid names
for(i in microbe.names){ 
  
  # print status
  print(paste("Hurry the fuck up:", i, "which is", which(microbe.names==i), "out of", no.microbes))
  
  # create temporary data matrix and model formula
  tmp <- df_metadata[, c(i,"Fiber","Protein")]
  fml <- as.formula(paste( i, "~", paste(c("Fiber*Protein"), collapse="+") ) )
  
  # assign fit to list by name
  fitlistnb2[[i]]<- glmmTMB(fml, family=nbinom2(link="log"), data=tmp)
  
  #print summary
  summariesnb2 <- print(summary(fitlistnb2[[i]])) 
  cofsnb2 <- coef(summary(fitlistnb2[[i]]))
  outputnb2 <- rep(do.call(rbind.data.frame, cofsnb2))
  review_nb2[[i]] <- outputnb2
}

LIST_nb2 <- do.call(rbind.data.frame, review_nb2)

write.csv(LIST_nb2, "MOVO glmmTMB Genus Summary.csv")



#########NETWORK CODE##############3
library(SpiecEasi)
library(NetCoMi)
library(phyloseq)
library(qgraph)
library(igraph)
library(vegan)
library(MCL)
library(qiime2R)
library(ggplot2)
library(phyloseqCompanion)
library(CoDaSeq)
library(microbiome)
library(DESeq2)
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
library(scales)
library(dplyr)
library(grid)
library(reshape2)
library(ape)
library(ashr)
library(limma)
library(edgeR)
library(splines)
library(RColorBrewer)

setwd("~/Desktop/RodentMicrobiome2024")

#####create your phyloseq object##########
##16S##
physeq_16S_MOVO <- qza_to_phyloseq(features="table-bacteria-MOVO.qza", taxonomy="taxonomy.qza", metadata="rodent-metadata-ut.txt")
rank_names(physeq_16S_MOVO)
table(tax_table(physeq_16S_MOVO)[, "Genus"], exclude = NULL) #view low abundant or uncharacterized taxa
physeq_16S_MOVO_sub <- subset_taxa(physeq_16S_MOVO, !is.na(Genus) & !Genus %in% c("", "uncharacterized")) #remove any with Genus = NA
# Compute prevalence of each feature, store as data.frame
prevdf_MOVO = apply(X = otu_table(physeq_16S_MOVO_sub),
                    MARGIN = ifelse(taxa_are_rows(physeq_16S_MOVO_sub), yes = 1, no = 2),
                    FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_MOVO = data.frame(Prevalence = prevdf_MOVO,
                         TotalAbundance = taxa_sums(physeq_16S_MOVO_sub),
                         tax_table(physeq_16S_MOVO_sub))
plyr::ddply(prevdf_MOVO, "Genus", function(df){cbind(mean(df$Prevalence),sum(df$Prevalence))}) #Shows prevalence of phyla, remove low prevalnence(?)
prevdf1_MOVO = subset(prevdf_MOVO, Genus %in% get_taxa_unique(physeq_16S_MOVO_sub, "Genus"))
ggplot(prevdf1_MOVO, aes(TotalAbundance, Prevalence / nsamples(physeq_16S_MOVO_sub),color=Genus)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Genus) + theme(legend.position="none")
prevalenceThreshold = 0.05* nsamples(physeq_16S_MOVO_sub)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1_MOVO)[(prevdf1_MOVO$Prevalence >= prevalenceThreshold)]
ps_16S_MOVO = prune_taxa(keepTaxa, physeq_16S_MOVO_sub)
# How many genera would be present after filtering?
length(get_taxa_unique(ps_16S_MOVO, taxonomic.rank = "Genus"))

# Agglomerate to genus level
ps_16S_MOVO_filter <- phyloseq::tax_glom(ps_16S_MOVO, taxrank = "Genus")
taxtab <- ps_16S_MOVO_filter@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Family)
miss_f <- which(taxtab[, "Family"] == "f__")
miss_g <- which(taxtab[, "Genus"] == "g__")

# Number unspecified genera
taxtab[miss_f, "Family"] <- paste0("f__", 1:length(miss_f))
taxtab[miss_g, "Genus"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g <- which(duplicated(taxtab[, "Genus"]) |
                  duplicated(taxtab[, "Genus"], fromLast = TRUE))

for(i in seq_along(taxtab)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Family"], ")")
  } else if(i %in% miss_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Order"], ")")
  }
  
  # Family names are added to duplicate genera
  if(i %in% dupl_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Order"], ")")
  }
}

ps_16S_MOVO_filter@tax_table@.Data <- taxtab
rownames(ps_16S_MOVO_filter@otu_table@.Data) <- taxtab[, "Genus"]
otu <- as.data.frame(otu_table(ps_16S_MOVO_filter))


##ITS##
physeq_ITS_MOVO <- qza_to_phyloseq(features="~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/table-sans-no-wild-fecals-MOVO-only.qza", taxonomy="~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/taxonomy.qza", metadata="~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/rodent-metadata-ut-full.txt")
rank_names(physeq_ITS_MOVO)
table(tax_table(physeq_ITS_MOVO)[, "Genus"], exclude = NULL) #view low abundant or uncharacterized taxa
physeq_ITS_MOVO_sub <- subset_taxa(physeq_ITS_MOVO, !is.na(Genus) & !Genus %in% c("", "uncharacterized")) #remove any with Phyla = NA
# Compute prevalence of each feature, store as data.frame
prevdf_ITS_MOVO = apply(X = otu_table(physeq_ITS_MOVO_sub),
                        MARGIN = ifelse(taxa_are_rows(physeq_ITS_MOVO_sub), yes = 1, no = 2),
                        FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_ITS_MOVO = data.frame(Prevalence = prevdf_ITS_MOVO,
                             TotalAbundance = taxa_sums(physeq_ITS_MOVO_sub),
                             tax_table(physeq_ITS_MOVO_sub))
plyr::ddply(prevdf_ITS_MOVO, "Genus", function(df){cbind(mean(df$Prevalence),sum(df$Prevalence))}) #Shows prevalence of phyla, remove low prevalnence(?)
prevdf1_ITS_MOVO = subset(prevdf_ITS_MOVO, Genus %in% get_taxa_unique(physeq_ITS_MOVO_sub, "Genus"))
ggplot(prevdf1_ITS_MOVO, aes(TotalAbundance, Prevalence / nsamples(physeq_ITS_MOVO_sub),color=Genus)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
prevalenceThreshold = 0.05* nsamples(physeq_ITS_MOVO_sub)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1_ITS_MOVO)[(prevdf1_ITS_MOVO$Prevalence >= prevalenceThreshold)]
ps_ITS_MOVO = prune_taxa(keepTaxa, physeq_ITS_MOVO_sub)
# How many genera would be present after filtering?
length(get_taxa_unique(ps_ITS_MOVO, taxonomic.rank = "Genus"))


ps_ITS_MOVO_filter <- phyloseq::tax_glom(ps_ITS_MOVO, taxrank = "Genus")
taxtab1 <- ps_ITS_MOVO_filter@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Family)
miss_f <- which(taxtab1[, "Family"] == "f__")
miss_g <- which(taxtab1[, "Genus"] == "g__")

# Number unspecified genera
taxtab1[miss_f, "Family"] <- paste0("f__", 1:length(miss_f))
taxtab1[miss_g, "Genus"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g1 <- which(duplicated(taxtab1[, "Genus"]) |
                   duplicated(taxtab1[, "Genus"], fromLast = TRUE))

for(i in seq_along(taxtab1)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab1[i, "Genus"] <- paste0(taxtab1[i, "Genus"], "(", taxtab1[i, "Family"], ")")
  } else if(i %in% miss_g){
    taxtab1[i, "Genus"] <- paste0(taxtab1[i, "Genus"], "(", taxtab1[i, "Order"], ")")
  }
  
  # Family names are added to duplicate genera
  if(i %in% dupl_g1){
    taxtab1[i,"Genus"] <- paste0(taxtab1[i,"Genus"],"(", taxtab1[i, "Family"], ")")
  }
}

ps_ITS_MOVO_filter@tax_table@.Data <- taxtab1
rownames(ps_ITS_MOVO_filter@otu_table@.Data) <- taxtab1[, "Genus"]
otu <- as.data.frame(otu_table(ps_ITS_MOVO_filter))



counts_16S <- as.matrix(t(phyloseq::otu_table(ps_16S_MOVO_filter)@.Data))
counts_16S <- t(counts_16S)
counts_ITS <- as.matrix(t(phyloseq::otu_table(ps_ITS_MOVO_filter)@.Data))
counts_ITS <- t(counts_ITS)

####FIBER#####
high_fiber_16S <- counts_16S[,c(1:3, 7, 9, 11:14,18:19,23,25:26,30:32,34,37,38,40)]
high_fiber_16S <- t(high_fiber_16S)
low_fiber_16S <- counts_16S[, c(4:6,8,10,15:17,20:22,24,27:29,33,35,36,39,41)]
low_fiber_16S <-t(low_fiber_16S)

high_fiber_ITS <- counts_ITS[,c(1:3, 7, 9, 11:14,18:19,23,25:26,30:32,34,37,38,40)]
high_fiber_ITS <- t(high_fiber_ITS)
low_fiber_ITS <- counts_ITS[, c(4:6,8,10,15:17,20:22,24,27:29,33,35,36,39,41)]
low_fiber_ITS <- t(low_fiber_ITS)


set.seed(123456)


# Run SpiecEasi and create association matrix for group 1
spiec_result_gr1 <- multi.spiec.easi(list(high_fiber_16S, high_fiber_ITS), 
                                     method='mb', nlambda=100, 
                                     lambda.min.ratio=1e-2, 
                                     pulsar.params = list(thresh = 0.05))
assoMat1 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec_result_gr1), mode = "ave")
assoMat1 <- as.matrix(assoMat1)

spiec_result_gr2 <- multi.spiec.easi(list(low_fiber_16S, low_fiber_ITS), 
                                     method='mb', nlambda=100, 
                                     lambda.min.ratio=1e-2, 
                                     pulsar.params = list(thresh = 0.05))

assoMat2 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec_result_gr2), mode = "ave")

assoMat2 <- as.matrix(assoMat2)


#Get taxa names
taxnames <- c(taxa_names(ps_16S_MOVO_filter), taxa_names(ps_ITS_MOVO_filter))

colnames(assoMat1) <- rownames(assoMat1) <- taxnames
diag(assoMat1) <- 1

colnames(assoMat2) <- rownames(assoMat2) <- taxnames
diag(assoMat2) <- 1


library(NetCoMi)
##
# Network construction (pass association matrices to netConstruct)
# - sparsMethod must be set to "none" because sparsification is already included in SpiecEasi
net_fiber_16S_ITS <- netConstruct(data = assoMat1, data2 = assoMat2, 
                                  dataType = "condDependence", 
                                  sparsMethod = "none")

# Network analysis
netprops_fiber_16S_ITS <- netAnalyze(net_fiber_16S_ITS, hubPar = "eigenvector")
nodeCols <- c(rep("lightblue", ntaxa(ps_16S_MOVO_filter)), rep("orange", ntaxa(ps_ITS_MOVO_filter)))
names(nodeCols) <- taxnames

plot(netprops_fiber_16S_ITS, 
     sameLayout = FALSE, 
     groupsChanged = TRUE,
     layoutGroup = "union",
     shortenLabels = "none",
     labelLength = 6L,
     nodeColor = "colorVec", 
     rmSingles = "all",
     colorVec = nodeCols,
     nodeSize = "eigenvector", 
     nodeSizeSpread = 4,
     labelScale = FALSE,
     cexNodes = 5, 
     cexLabels = .2,
     cexHubLabels = .6,
     hubBorderCol = "black",
     hubTransp = 0,
     labels= FALSE,
     cexTitle = 2,
     groupNames = c("High Fiber", "Low Fiber"))


netcomp_fiber_16S_ITS <- netCompare(netprops_fiber_16S_ITS, permTest = FALSE, nPerm = 100, storeAssoPerm = TRUE, fileStoreAssoPerm = "assoPerm_comp", storeCountsPerm = FALSE, seed = 123456)
summary(netcomp_fiber_16S_ITS, groupNames = c("High Fiber", "Low Fiber"))

net_single<- netConstruct(data = assoMat2,
                          dataType = "condDependence", 
                          sparsMethod = "none",
                          measure = "spieceasi")

props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

summary(props_single, numbNodes = 5L)


#####PROTEIN######
high_protein_16S <- counts_16S[,c(2,5,9,11,13:15,19,21,22,24:26,31:33,35,36,38,39,41)]
high_protein_16S <- t(high_protein_16S)
low_protein_16S <- counts_16S[, c(1,3,4,6:8,10,12,16:18,20,23,27:30,34,37,40)]
low_protein_16S <-t(low_protein_16S)

high_protein_ITS <- counts_ITS[,c(2,5,9,11,13:15,19,21,22,24:26,31:33,35,36,38,39,41)]
high_protein_ITS <- t(high_protein_ITS)
low_protein_ITS <- counts_ITS[, c(1,3,4,6:8,10,12,16:18,20,23,27:30,34,37,40)]
low_protein_ITS <- t(low_protein_ITS)


set.seed(123456)


# Run SpiecEasi and create association matrix for group 1
spiec_result_gr1 <- multi.spiec.easi(list(high_protein_16S, high_protein_ITS), 
                                     method='mb', nlambda=100, 
                                     lambda.min.ratio=1e-2, 
                                     pulsar.params = list(thresh = 0.05))
assoMat1 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec_result_gr1), mode = "ave")
assoMat1 <- as.matrix(assoMat1)

spiec_result_gr2 <- multi.spiec.easi(list(low_protein_16S, low_protein_ITS), 
                                     method='mb', nlambda=40, 
                                     lambda.min.ratio=1e-2, 
                                     pulsar.params = list(thresh = 0.05))

assoMat2 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec_result_gr2), mode = "ave")

assoMat2 <- as.matrix(assoMat2)


#Get taxa names
taxnames <- c(taxa_names(ps_16S_MOVO_filter), taxa_names(ps_ITS_MOVO_filter))

colnames(assoMat1) <- rownames(assoMat1) <- taxnames
diag(assoMat1) <- 1

colnames(assoMat2) <- rownames(assoMat2) <- taxnames
diag(assoMat2) <- 1


library(NetCoMi)
##
# Network construction (pass association matrices to netConstruct)
# - sparsMethod must be set to "none" because sparsification is already included in SpiecEasi
net_protein_16S_ITS <- netConstruct(data = assoMat1, data2 = assoMat2, 
                                    dataType = "condDependence", 
                                    sparsMethod = "none")

# Network analysis
netprops_protein_16S_ITS <- netAnalyze(net_protein_16S_ITS, hubPar = "eigenvector")

nodeCols <- c(rep("lightblue", ntaxa(ps_16S_MOVO_filter)), rep("orange", ntaxa(ps_ITS_MOVO_filter)))
names(nodeCols) <- taxnames


plot(netprops_protein_16S_ITS, 
     sameLayout = FALSE, 
     layoutGroup = "union",
     shortenLabels = "none",
     nodeColor = "colorVec", 
     colorVec = nodeCols,
     rmSingles=TRUE,
     nodeSize = "eigen", 
     nodeSizeSpread = 2,
     labelScale = FALSE,
     cexNodes = 4, 
     cexLabels = .5,
     cexHubLabels = .75,
     cexTitle = 3,
     groupNames = c("High Protein", "Low Protein"))


netcomp_protein_16S_ITS <- netCompare(netprops_protein_16S_ITS, permTest = FALSE)
summary(netcomp_protein_16S_ITS, groupNames = c("High Protein", "Low Protein"))

net_single<- netConstruct(data = assoMat1,
                          dataType = "condDependence", 
                          sparsMethod = "none",
                          measure = "spieceasi")

props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

summary(props_single, numbNodes = 5L)



###WFMO#######
#####create your phyloseq object##########
##16S##
physeq_16S_WFMO <- qza_to_phyloseq(features="table-bacteria-WFMO.qza", taxonomy="taxonomy.qza", metadata="rodent-metadata-ky.txt")
rank_names(physeq_16S_WFMO)
table(tax_table(physeq_16S_WFMO)[, "Genus"], exclude = NULL) #view low abundant or uncharacterized taxa
physeq_16S_WFMO_sub <- subset_taxa(physeq_16S_WFMO, !is.na(Genus) & !Genus %in% c("", "uncharacterized")) #remove any with Phyla = NA
# Compute prevalence of each feature, store as data.frame
prevdf_WFMO = apply(X = otu_table(physeq_16S_WFMO_sub),
                    MARGIN = ifelse(taxa_are_rows(physeq_16S_WFMO_sub), yes = 1, no = 2),
                    FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_WFMO = data.frame(Prevalence = prevdf_WFMO,
                         TotalAbundance = taxa_sums(physeq_16S_WFMO_sub),
                         tax_table(physeq_16S_WFMO_sub))
plyr::ddply(prevdf_WFMO, "Genus", function(df){cbind(mean(df$Prevalence),sum(df$Prevalence))}) #Shows prevalence of phyla, remove low prevalnence(?)
prevdf1_WFMO = subset(prevdf_WFMO, Genus %in% get_taxa_unique(physeq_16S_WFMO_sub, "Genus"))
ggplot(prevdf1_WFMO, aes(TotalAbundance, Prevalence / nsamples(physeq_16S_WFMO_sub),color=Genus)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Genus) + theme(legend.position="none")
prevalenceThreshold = 0.1* nsamples(physeq_16S_WFMO_sub)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1_WFMO)[(prevdf1_WFMO$Prevalence >= prevalenceThreshold)]
ps_16S_WFMO = prune_taxa(keepTaxa, physeq_16S_WFMO_sub)
# How many genera would be present after filtering?
length(get_taxa_unique(ps_16S_WFMO, taxonomic.rank = "Genus"))




##ITS##
physeq_ITS_WFMO <- qza_to_phyloseq(features="~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/table-sans-no-wild-fecals-WFMO-only.qza", taxonomy= "~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/taxonomy.qza", metadata="~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/rodent-metadata-ky-full.txt")
rank_names(physeq_ITS_WFMO)
table(tax_table(physeq_ITS_WFMO)[, "Genus"], exclude = NULL) #view low abundant or uncharacterized taxa
physeq_ITS_WFMO_sub <- subset_taxa(physeq_ITS_WFMO, !is.na(Genus) & !Genus %in% c("", "uncharacterized")) #remove any with Phyla = NA
# Compute prevalence of each feature, store as data.frame
prevdf_ITS_WFMO = apply(X = otu_table(physeq_ITS_WFMO_sub),
                        MARGIN = ifelse(taxa_are_rows(physeq_ITS_WFMO_sub), yes = 1, no = 2),
                        FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_ITS_WFMO = data.frame(Prevalence = prevdf_ITS_WFMO,
                             TotalAbundance = taxa_sums(physeq_ITS_WFMO_sub),
                             tax_table(physeq_ITS_WFMO_sub))
plyr::ddply(prevdf_ITS_WFMO, "Genus", function(df){cbind(mean(df$Prevalence),sum(df$Prevalence))}) #Shows prevalence of phyla, remove low prevalnence(?)
prevdf1_ITS_WFMO = subset(prevdf_ITS_WFMO, Genus %in% get_taxa_unique(physeq_ITS_WFMO_sub, "Genus"))
ggplot(prevdf1_ITS_WFMO, aes(TotalAbundance, Prevalence / nsamples(physeq_ITS_WFMO_sub),color=Genus)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Genus) + theme(legend.position="none")
prevalenceThreshold = 0.1* nsamples(physeq_ITS_WFMO_sub)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1_ITS_WFMO)[(prevdf1_ITS_WFMO$Prevalence >= prevalenceThreshold)]
ps_ITS_WFMO = prune_taxa(keepTaxa, physeq_ITS_WFMO_sub)
# How many genera would be present after filtering?
length(get_taxa_unique(ps_ITS_WFMO, taxonomic.rank = "Genus"))


#Merge#
ps_ITS_16S_WFMO <- merge_phyloseq(ps_16S_WFMO, ps_ITS_WFMO)
ps_ITS_16S_WFMO_genus <- aggregate_taxa(ps_ITS_16S_WFMO, "Genus", verbose=FALSE)

otu <- otu_table(ps_ITS_16S_WFMO_genus)
otu <- as.data.frame(otu)
otu <- t(otu)
metadata <- read.table("rodent-metadata-ky.txt", header=TRUE)
otu <- cbind(otu, metadata[, 5:6])
write.csv(otu, "~/Desktop/otu_wfmo.csv")

plot <- ggplot(otu, aes(x = Fiber, y = `Bacteroides`, fill = Protein)) +
  geom_boxplot() +
  theme_bw() 
plot

p <- ggplot(otu, aes(x=Protein, y=`Bacteroides`, fill =Fiber)) + geom_violin(trim=FALSE)
p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +  geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.9)) + scale_fill_manual(values=c('dodgerblue','blue')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Protein", y="Observed ASVs") + theme(plot.background = element_rect(fill = 'black'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill='black'), axis.line = element_line(colour = "white"))+ theme(axis.text.x=element_text(colour="white"), axis.text.y=element_text(colour="white"), axis.title.x=element_text(colour="white"), axis.title.y=element_text(colour="white")) 



###tax_glom to get to genus level###
ps_ITS_16S_WFMO_genus <- tax_glom(ps_ITS_16S_WFMO, taxrank="Genus", NArm=FALSE)
ps_ITS_16S_WFMO_genus_OTU <- as(otu_table(ps_ITS_16S_WFMO_genus), "matrix")


######DESeq2 Differential Abundance######
#DESeq Analysis
dds = phyloseq_to_deseq2(ps_ITS_16S_WFMO_genus, ~ Fiber*Protein)  ####NOTE: CANNOT USE ps_g_16S_MOVO due to an error, but can use phyloseq object####
dds$Fiber <- relevel(dds$Fiber, "Low")
dds$Protein <- relevel(dds$Protein, "Low")

#Calculate geometric means
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, type="poscounts")


#Run DESeq2
dds = DESeq(dds, fitType="local")
resultsNames(dds) # "Intercept" "Fiber_High_vs_Low" "Protein_High_vs_Low" "FiberHigh.ProteinHigh"
res = results(dds, name="FiberHigh.ProteinHigh")
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_ITS_16S_WFMO_genus)[rownames(sigtab), ], "matrix"))
head(sigtab)
summary(sigtab)
summary(res)

#Shrink
shrink <- lfcShrink(dds, coef=c("FiberHigh.ProteinHigh"), sigtab, type="ashr", lfcThreshold = 0)
shrink <- shrink[order(shrink$padj, na.last=NA), ]
sigtab1 = shrink[(shrink$padj < alpha), ]
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(ps_ITS_16S_WFMO_genus)[rownames(sigtab1), ], "matrix"))
summary(sigtab1)
write.csv(sigtab1, "DESEq GRMO Interaction.csv")

DESeq2::plotMA(shrink, ylim=c(-30,30)) #MAplot to look at 
shrink1 = as.data.frame(shrink)
shrink1 = mutate(shrink1, sig=ifelse(shrink1$padj<0.1, "FDR<0.1", "NotSig"))
shrink1[which(abs(shrink1$log2FoldChange)<1.0), "sig"] = "NotSig"
shrink1 = cbind(as(shrink1, "data.frame"), as(tax_table(ps_ITS_16S_WFMO_genus)[rownames(shrink1), ], "matrix"))

write.csv(shrink1, "~/Desktop/WFMO_DEseq2_Interaction-SILVA.csv")


####Taxon Filtering####
#get 16S genus#
# Agglomerate to genus level
ps_g_16S_WFMO <- phyloseq::tax_glom(ps_16S_WFMO, taxrank = "Genus")
taxtab <- ps_g_16S_WFMO@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Family)
miss_f <- which(taxtab[, "Family"] == "f__")
miss_g <- which(taxtab[, "Genus"] == "g__")

# Number unspecified genera
taxtab[miss_f, "Family"] <- paste0("f__", 1:length(miss_f))
taxtab[miss_g, "Genus"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g <- which(duplicated(taxtab[, "Genus"]) |
                  duplicated(taxtab[, "Genus"], fromLast = TRUE))

for(i in seq_along(taxtab)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Order"], ")")
  } else if(i %in% miss_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Family"], ")")
  }
  
  # Family names are added to duplicate genera
  if(i %in% dupl_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Family"], ")")
  }
}

ps_g_16S_WFMO@tax_table@.Data <- taxtab
rownames(ps_g_16S_WFMO@otu_table@.Data) <- taxtab[, "Genus"]

#get ITS genus#
ps_genus_ITS_WFMO <- phyloseq::tax_glom(ps_ITS_WFMO, taxrank = "Genus")
taxtab_ITS_WFMO <- ps_genus_ITS_WFMO@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Family)
miss_f <- which(taxtab_ITS_WFMO[, "Family"] == "f__")
miss_g <- which(taxtab_ITS_WFMO[, "Genus"] == "g__")

# Number unspecified genera
taxtab_ITS_WFMO[miss_f, "Family"] <- paste0("f__", 1:length(miss_f))
taxtab_ITS_WFMO[miss_g, "Genus"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g1 <- which(duplicated(taxtab_ITS_WFMO[, "Genus"]) |
                   duplicated(taxtab_ITS_WFMO[, "Genus"], fromLast = TRUE))

for(i in seq_along(taxtab_ITS_WFMO)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab_ITS_WFMO[i, "Genus"] <- paste0(taxtab_ITS_WFMO[i, "Genus"], "(", taxtab_ITS_WFMO[i, "Order"], ")")
  } else if(i %in% miss_g){
    taxtab_ITS_WFMO[i, "Genus"] <- paste0(taxtab_ITS_WFMO[i, "Genus"], "(", taxtab_ITS_WFMO[i, "Family"], ")")
  }
  
  # Family names are added to duplicate genera
  if(i %in% dupl_g1){
    taxtab_ITS_WFMO[i, "Genus"] <- paste0(taxtab_ITS_WFMO[i, "Genus"], "(", taxtab_ITS_WFMO[i, "Family"], ")")
  }
}

ps_genus_ITS_WFMO@tax_table@.Data <- taxtab_ITS_WFMO
rownames(ps_genus_ITS_WFMO@otu_table@.Data) <- taxtab_ITS_WFMO[, "Genus"]


####merged######
ps_ITS_16S_WFMO_g <- phyloseq::tax_glom(ps_ITS_16S_WFMO, taxrank = "Genus")
taxtab_ITS_16S_WFMO_g <- ps_ITS_16S_WFMO_g@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Family)
miss_f <- which(taxtab_ITS_16S_WFMO_g[, "Family"] == "f__")
miss_g <- which(taxtab_ITS_16S_WFMO_g[, "Genus"] == "g__")

# Number unspecified genera
taxtab_ITS_16S_WFMO_g[miss_f, "Family"] <- paste0("f__", 1:length(miss_f))
taxtab_ITS_16S_WFMO_g[miss_g, "Genus"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g1 <- which(duplicated(taxtab_ITS_16S_WFMO_g[, "Genus"]) |
                   duplicated(taxtab_ITS_16S_WFMO_g[, "Genus"], fromLast = TRUE))

for(i in seq_along(taxtab_ITS_16S_WFMO_g)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab_ITS_16S_WFMO_g[i, "Genus"] <- paste0(taxtab_ITS_16S_WFMO_g[i, "Genus"], "(", taxtab_ITS_16S_WFMO_g[i, "Order"], ")")
  } else if(i %in% miss_g){
    taxtab_ITS_16S_WFMO_g[i, "Genus"] <- paste0(taxtab_ITS_16S_WFMO_g[i, "Genus"], "(", taxtab_ITS_16S_WFMO_g[i, "Family"], ")")
  }
  
  # Family names are added to duplicate genera
  if(i %in% dupl_g1){
    taxtab_ITS_16S_WFMO_g[i, "Genus"] <- paste0(taxtab_ITS_16S_WFMO_g[i, "Genus"], "(", taxtab_ITS_16S_WFMO_g[i, "Family"], ")")
  }
}

ps_ITS_16S_WFMO_g@tax_table@.Data <- taxtab_ITS_16S_WFMO_g
rownames(ps_ITS_16S_WFMO_g@otu_table@.Data) <- taxtab_ITS_16S_WFMO_g[, "Genus"]
ITS_16S_WFMO_genus_table <- as(otu_table(ps_ITS_16S_WFMO_g), "matrix")


######glmmTMB########
library(broom)
library(car)
library(magicfor)
library(glmmTMB)
library(AICcmodavg)
library(bbmle)

data <- ITS_16S_WFMO_genus_table
df_transpose <- as.data.frame(t(as.matrix(data)))
metadata <- read.table("rodent-metadata-ky.txt", header=TRUE)
df_metadata <- cbind(df_transpose, metadata$Protein)
df_metadata <- cbind(df_metadata, metadata$Fiber)
write.csv(df_metadata, "~/Desktop/RodentMicrobiome2024/df_metdata.csv")
df_metadata <- read.csv("df_metdata.csv")
fiber <- df_metadata$Fiber
protein <- df_metadata$Protein


# names of lipids
microbe.names <- colnames(df_metadata)[4:114]
no.microbes <- length(microbe.names)

# create a named list to hold the fitted models
fitlistnb2 <- as.list(1:no.microbes)
names(fitlistnb2) <- microbe.names
review_nb2 <- vector("list", length=110)

# loop over lipid names
for(i in microbe.names){ 
  
  # print status
  print(paste("Hurry the fuck up:", i, "which is", which(microbe.names==i), "out of", no.microbes))
  
  # create temporary data matrix and model formula
  tmp <- df_metadata[, c(i,"Fiber","Protein")]
  fml <- as.formula(paste( i, "~", paste(c("Fiber*Protein"), collapse="+") ) )
  
  # assign fit to list by name
  fitlistnb2[[i]]<- glmmTMB(fml,family=nbinom2(link="log"), data=tmp)
  
  #print summary
  summariesnb2 <- print(summary(fitlistnb2[[i]])) 
  cofsnb2 <- coef(summary(fitlistnb2[[i]]))
  outputnb2 <- rep(do.call(rbind.data.frame, cofsnb2))
  review_nb2[[i]] <- outputnb2
}

LIST_nb2 <- do.call(rbind.data.frame, review_nb2)

write.csv(LIST_nb2, "WFMO SILVA glmmTMB Genus Summary.csv")



#########NETWORK CODE##############3
library(SpiecEasi)
library(NetCoMi)
library(phyloseq)
library(qgraph)
library(igraph)
library(vegan)
library(MCL)
library(qiime2R)
library(ggplot2)
library(phyloseqCompanion)
library(CoDaSeq)
library(microbiome)
library(DESeq2)
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
library(scales)
library(dplyr)
library(grid)
library(reshape2)
library(ape)
library(ashr)
library(limma)
library(edgeR)
library(splines)
library(RColorBrewer)

setwd("~/Desktop/RodentMicrobiome2024")

#####create your phyloseq object##########
##16S##
physeq_16S_WFMO <- qza_to_phyloseq(features="table-bacteria-WFMO.qza", taxonomy="taxonomy.qza", metadata="rodent-metadata-ky.txt")
rank_names(physeq_16S_WFMO)
table(tax_table(physeq_16S_WFMO)[, "Genus"], exclude = NULL) #view low abundant or uncharacterized taxa
physeq_16S_WFMO_sub <- subset_taxa(physeq_16S_WFMO, !is.na(Genus) & !Genus %in% c("", "uncharacterized")) #remove any with Phyla = NA
# Compute prevalence of each feature, store as data.frame
prevdf_WFMO = apply(X = otu_table(physeq_16S_WFMO_sub),
                    MARGIN = ifelse(taxa_are_rows(physeq_16S_WFMO_sub), yes = 1, no = 2),
                    FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_WFMO = data.frame(Prevalence = prevdf_WFMO,
                         TotalAbundance = taxa_sums(physeq_16S_WFMO_sub),
                         tax_table(physeq_16S_WFMO_sub))
plyr::ddply(prevdf_WFMO, "Genus", function(df){cbind(mean(df$Prevalence),sum(df$Prevalence))}) #Shows prevalence of phyla, remove low prevalnence(?)
prevdf1_WFMO = subset(prevdf_WFMO, Genus %in% get_taxa_unique(physeq_16S_WFMO_sub, "Genus"))
ggplot(prevdf1_WFMO, aes(TotalAbundance, Prevalence / nsamples(physeq_16S_WFMO_sub),color=Genus)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Genus) + theme(legend.position="none")
prevalenceThreshold = 0.1* nsamples(physeq_16S_WFMO_sub)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1_WFMO)[(prevdf1_WFMO$Prevalence >= prevalenceThreshold)]
ps_16S_WFMO = prune_taxa(keepTaxa, physeq_16S_WFMO_sub)
# How many genera would be present after filtering?
length(get_taxa_unique(ps_16S_WFMO, taxonomic.rank = "Genus"))

# Agglomerate to genus level
ps_16S_WFMO_filter <- phyloseq::tax_glom(ps_16S_WFMO, taxrank = "Genus")
taxtab <- ps_16S_WFMO_filter@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Family)
miss_f <- which(taxtab[, "Family"] == "f__")
miss_g <- which(taxtab[, "Genus"] == "g__")

# Number unspecified genera
taxtab[miss_f, "Family"] <- paste0("f__", 1:length(miss_f))
taxtab[miss_g, "Genus"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g <- which(duplicated(taxtab[, "Genus"]) |
                  duplicated(taxtab[, "Genus"], fromLast = TRUE))

for(i in seq_along(taxtab)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Family"], ")")
  } else if(i %in% miss_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Order"], ")")
  }
  
  # Family names are added to duplicate genera
  if(i %in% dupl_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Order"], ")")
  }
}

ps_16S_WFMO_filter@tax_table@.Data <- taxtab
rownames(ps_16S_WFMO_filter@otu_table@.Data) <- taxtab[, "Genus"]
otu <- as.data.frame(otu_table(ps_16S_WFMO_filter))

##ITS##
physeq_ITS_WFMO <- qza_to_phyloseq(features="~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/table-sans-no-wild-fecals-WFMO-only.qza", taxonomy= "~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/taxonomy.qza", metadata="~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/rodent-metadata-ky-full.txt")
rank_names(physeq_ITS_WFMO)
table(tax_table(physeq_ITS_WFMO)[, "Genus"], exclude = NULL) #view low abundant or uncharacterized taxaRodent Microbiome Differential Abundance and Network Analysis.R
physeq_ITS_WFMO_sub <- subset_taxa(physeq_ITS_WFMO, !is.na(Genus) & !Genus %in% c("", "uncharacterized")) #remove any with Phyla = NA
# Compute prevalence of each feature, store as data.frame
prevdf_ITS_WFMO = apply(X = otu_table(physeq_ITS_WFMO_sub),
                        MARGIN = ifelse(taxa_are_rows(physeq_ITS_WFMO_sub), yes = 1, no = 2),
                        FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_ITS_WFMO = data.frame(Prevalence = prevdf_ITS_WFMO,
                             TotalAbundance = taxa_sums(physeq_ITS_WFMO_sub),
                             tax_table(physeq_ITS_WFMO_sub))
plyr::ddply(prevdf_ITS_WFMO, "Genus", function(df){cbind(mean(df$Prevalence),sum(df$Prevalence))}) #Shows prevalence of phyla, remove low prevalnence(?)
prevdf1_ITS_WFMO = subset(prevdf_ITS_WFMO, Genus %in% get_taxa_unique(physeq_ITS_WFMO_sub, "Genus"))
ggplot(prevdf1_ITS_WFMO, aes(TotalAbundance, Prevalence / nsamples(physeq_ITS_WFMO_sub),color=Genus)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Genus) + theme(legend.position="none")
prevalenceThreshold = 0.1* nsamples(physeq_ITS_WFMO_sub)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1_ITS_WFMO)[(prevdf1_ITS_WFMO$Prevalence >= prevalenceThreshold)]
ps_ITS_WFMO = prune_taxa(keepTaxa, physeq_ITS_WFMO_sub)
# How many genera would be present after filtering?
length(get_taxa_unique(ps_ITS_WFMO, taxonomic.rank = "Genus"))


ps_ITS_WFMO_filter <- phyloseq::tax_glom(ps_ITS_WFMO, taxrank = "Genus")
taxtab1 <- ps_ITS_WFMO_filter@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Family)
miss_f <- which(taxtab1[, "Family"] == "f__")
miss_g <- which(taxtab1[, "Genus"] == "g__")

# Number unspecified genera
taxtab1[miss_f, "Family"] <- paste0("f__", 1:length(miss_f))
taxtab1[miss_g, "Genus"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g1 <- which(duplicated(taxtab1[, "Genus"]) |
                   duplicated(taxtab1[, "Genus"], fromLast = TRUE))

for(i in seq_along(taxtab1)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab1[i, "Genus"] <- paste0(taxtab1[i, "Genus"], "(", taxtab1[i, "Family"], ")")
  } else if(i %in% miss_g){
    taxtab1[i, "Genus"] <- paste0(taxtab1[i, "Genus"], "(", taxtab1[i, "Order"], ")")
  }
  
  # Family names are added to duplicate genera
  if(i %in% dupl_g1){
    taxtab1[i,"Genus"] <- paste0(taxtab1[i,"Genus"],"(", taxtab1[i, "Family"], ")")
  }
}

ps_ITS_WFMO_filter@tax_table@.Data <- taxtab1
rownames(ps_ITS_WFMO_filter@otu_table@.Data) <- taxtab1[, "Genus"]
otu <- as.data.frame(otu_table(ps_ITS_WFMO_filter))



counts_16S <- as.matrix(t(phyloseq::otu_table(ps_16S_WFMO_filter)@.Data))
counts_16S <- t(counts_16S)
counts_ITS <- as.matrix(t(phyloseq::otu_table(ps_ITS_WFMO_filter)@.Data))
counts_ITS <- t(counts_ITS)
write.csv(counts_16S, "~/Desktop/counts_16S_WFMO.csv")
write.csv(counts_ITS, "~/Desktop/counts_ITS_WFMO.csv")

####FIBER#####
high_fiber_16S <- counts_16S[,c(2, 5:7, 9, 11:15,20,22,23,25,26,30,31,35,36,37)]
high_fiber_16S <- t(high_fiber_16S)
low_fiber_16S <- counts_16S[, c(1,3,4,8,10,16:19,21,24,27,28,29,32:34,38:40)]
low_fiber_16S <-t(low_fiber_16S)

high_fiber_ITS <- counts_ITS[,c(2, 5:7, 9, 11:15,20,22,23,25,26,30,31,35,36,37)]
high_fiber_ITS <- t(high_fiber_ITS)
low_fiber_ITS <- counts_ITS[, c(1,3,4,8,10,16:19,21,24,27,28,29,32:34,38:40)]
low_fiber_ITS <- t(low_fiber_ITS)


set.seed(123456)


# Run SpiecEasi and create association matrix for group 1
spiec_result_gr1 <- multi.spiec.easi(list(high_fiber_16S, high_fiber_ITS), 
                                     method='mb', nlambda=100, 
                                     lambda.min.ratio=1e-2, 
                                     pulsar.params = list(thresh = 0.05))
assoMat1 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec_result_gr1), mode = "ave")
assoMat1 <- as.matrix(assoMat1)

spiec_result_gr2 <- multi.spiec.easi(list(low_fiber_16S, low_fiber_ITS), 
                                     method='mb', nlambda=100, 
                                     lambda.min.ratio=1e-2, 
                                     pulsar.params = list(thresh = 0.05))

assoMat2 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec_result_gr2), mode = "ave")

assoMat2 <- as.matrix(assoMat2)


#Get taxa names
taxnames <- c(taxa_names(ps_16S_WFMO_filter), taxa_names(ps_ITS_WFMO_filter))

colnames(assoMat1) <- rownames(assoMat1) <- taxnames
diag(assoMat1) <- 1

colnames(assoMat2) <- rownames(assoMat2) <- taxnames
diag(assoMat2) <- 1


library(NetCoMi)
##
# Network construction (pass association matrices to netConstruct)
# - sparsMethod must be set to "none" because sparsification is already included in SpiecEasi
net_fiber_16S_ITS <- netConstruct(data = assoMat1, data2 = assoMat2, 
                                  dataType = "condDependence", 
                                  sparsMethod = "none")


# Network analysis
netprops_fiber_16S_ITS <- netAnalyze(net_fiber_16S_ITS, hubPar = "eigenvector")
nodeCols <- c(rep("lightblue", ntaxa(ps_16S_WFMO_filter)), rep("orange", ntaxa(ps_ITS_WFMO_filter)))
names(nodeCols) <- taxnames

plot(netprops_fiber_16S_ITS, 
     sameLayout = FALSE, 
     layoutGroup = "union",
     shortenLabels = "none",
     nodeColor = "colorVec", 
     colorVec = nodeCols,
     rmSingles=TRUE,
     nodeSize = "eigen", 
     nodeSizeSpread = 10,
     labelScale = FALSE,
     cexNodes = 2, 
     cexLabels = .5,
     cexHubLabels = .6,
     cexTitle = 3,
     groupNames = c("High Fiber", "Low Fiber"))


netcomp_fiber_16S_ITS <- netCompare(netprops_fiber_16S_ITS, permTest = FALSE)
summary(netcomp_fiber_16S_ITS, groupNames = c("High Fiber", "Low Fiber"))


net_single<- netConstruct(data = assoMat2,
                          dataType = "condDependence", 
                          sparsMethod = "none",
                          measure = "spieceasi")

props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

summary(props_single, numbNodes = 5L)


#####PROTEIN######
high_protein_16S <- counts_16S[,c(1,2,3,5,7,10,12,14,18,20,22,24,27,29:32,36,38,40)]
high_protein_16S <- t(high_protein_16S)
low_protein_16S <- counts_16S[, c(4,6,8,9,11,13,15:17,19,21,23,25,26,28,33:35,37,39)]
low_protein_16S <-t(low_protein_16S)

high_protein_ITS <- counts_ITS[,c(1,2,3,5,7,10,12,14,18,20,22,24,27,29:32,36,38,40)]
high_protein_ITS <- t(high_protein_ITS)
low_protein_ITS <- counts_ITS[, c(4,6,8,9,11,13,15:17,19,21,23,25,26,28,33:35,37,39)]
low_protein_ITS <- t(low_protein_ITS)


set.seed(123456)


# Run SpiecEasi and create association matrix for group 1
spiec_result_gr1 <- multi.spiec.easi(list(high_protein_16S, high_protein_ITS), 
                                     method='mb', nlambda=100, 
                                     lambda.min.ratio=1e-2, 
                                     pulsar.params = list(thresh = 0.05))
assoMat1 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec_result_gr1), mode = "ave")
assoMat1 <- as.matrix(assoMat1)

spiec_result_gr2 <- multi.spiec.easi(list(low_protein_16S, low_protein_ITS), 
                                     method='mb', nlambda=100, 
                                     lambda.min.ratio=1e-2, 
                                     pulsar.params = list(thresh = 0.05))

assoMat2 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec_result_gr2), mode = "ave")

assoMat2 <- as.matrix(assoMat2)


#Get taxa names
taxnames <- c(taxa_names(ps_16S_WFMO_filter), taxa_names(ps_ITS_WFMO_filter))

colnames(assoMat1) <- rownames(assoMat1) <- taxnames
diag(assoMat1) <- 1

colnames(assoMat2) <- rownames(assoMat2) <- taxnames
diag(assoMat2) <- 1


library(NetCoMi)
##
# Network construction (pass association matrices to netConstruct)
# - sparsMethod must be set to "none" because sparsification is already included in SpiecEasi
net_protein_16S_ITS <- netConstruct(data = assoMat1, data2 = assoMat2, 
                                    dataType = "condDependence", 
                                    sparsMethod = "none")

# Network analysis
netprops_protein_16S_ITS <- netAnalyze(net_protein_16S_ITS, hubPar = "eigenvector")
nodeCols <- c(rep("lightblue", ntaxa(ps_16S_WFMO_filter)), rep("orange", ntaxa(ps_ITS_WFMO_filter)))
names(nodeCols) <- taxnames


plot(netprops_protein_16S_ITS, 
     sameLayout = FALSE,
     rmSingles = TRUE,
     layoutGroup = "union",
     shortenLabels = "none",
     nodeColor = "colorVec", 
     colorVec = nodeCols,
     nodeSize = "eigen", 
     nodeSizeSpread = 10,
     labelScale = FALSE,
     cexNodes = 4, 
     cexLabels = .5,
     cexHubLabels = .5,
     cexTitle = 3,
     groupNames = c("High Protein", "Low Protein"))


netcomp_protein_16S_ITS <- netCompare(netprops_protein_16S_ITS, permTest = FALSE)
summary(netcomp_protein_16S_ITS, groupNames = c("High Protein", "Low Protein"))

net_single<- netConstruct(data = assoMat2,
                          dataType = "condDependence", 
                          sparsMethod = "none",
                          measure = "spieceasi")

props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

summary(props_single, numbNodes = 5L)


####GRMO#####
####create your phyloseq object##########
##16S##
physeq_16S_GRMO <- qza_to_phyloseq(features="table-bacteria-lab-GRMO.qza", taxonomy="taxonomy.qza", metadata="rodent-metadata-az.txt")
rank_names(physeq_16S_GRMO)
table(tax_table(physeq_16S_GRMO)[, "Genus"], exclude = NULL) #view low abundant or uncharacterized taxa
physeq_16S_GRMO_sub <- subset_taxa(physeq_16S_GRMO, !is.na(Genus) & !Genus %in% c("", "uncharacterized")) #remove any with Phyla = NA
# Compute prevalence of each feature, store as data.frame
prevdf_GRMO = apply(X = otu_table(physeq_16S_GRMO_sub),
                    MARGIN = ifelse(taxa_are_rows(physeq_16S_GRMO_sub), yes = 1, no = 2),
                    FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_GRMO = data.frame(Prevalence = prevdf_GRMO,
                         TotalAbundance = taxa_sums(physeq_16S_GRMO_sub),
                         tax_table(physeq_16S_GRMO_sub))
plyr::ddply(prevdf_GRMO, "Genus", function(df){cbind(mean(df$Prevalence),sum(df$Prevalence))}) #Shows prevalence of phyla, remove low prevalnence(?)
prevdf1_GRMO = subset(prevdf_GRMO, Genus %in% get_taxa_unique(physeq_16S_GRMO_sub, "Genus"))
ggplot(prevdf1_GRMO, aes(TotalAbundance, Prevalence / nsamples(physeq_16S_GRMO_sub),color=Genus)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Genus) + theme(legend.position="none")
prevalenceThreshold = 0.05* nsamples(physeq_16S_GRMO_sub)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1_GRMO)[(prevdf1_GRMO$Prevalence >= prevalenceThreshold)]
ps_16S_GRMO = prune_taxa(keepTaxa, physeq_16S_GRMO_sub)
# How many genera would be present after filtering?
length(get_taxa_unique(ps_16S_GRMO, taxonomic.rank = "Genus"))


##ITS##
physeq_ITS_GRMO <- qza_to_phyloseq(features="~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/table-sans-no-wild-fecals-GRMO-only.qza", taxonomy= "~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/taxonomy.qza", metadata="~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/rodent-metadata-az-full.txt")
rank_names(physeq_ITS_GRMO)
table(tax_table(physeq_ITS_GRMO)[, "Genus"], exclude = NULL) #view low abundant or uncharacterized taxa
physeq_ITS_GRMO_sub <- subset_taxa(physeq_ITS_GRMO, !is.na(Genus) & !Genus %in% c("", "uncharacterized")) #remove any with Phyla = NA
# Compute prevalence of each feature, store as data.frame
prevdf_ITS_GRMO = apply(X = otu_table(physeq_ITS_GRMO_sub),
                        MARGIN = ifelse(taxa_are_rows(physeq_ITS_GRMO_sub), yes = 1, no = 2),
                        FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_ITS_GRMO = data.frame(Prevalence = prevdf_ITS_GRMO,
                             TotalAbundance = taxa_sums(physeq_ITS_GRMO_sub),
                             tax_table(physeq_ITS_GRMO_sub))
plyr::ddply(prevdf_ITS_GRMO, "Genus", function(df){cbind(mean(df$Prevalence),sum(df$Prevalence))}) #Shows prevalence of phyla, remove low prevalnence(?)
prevdf1_ITS_GRMO = subset(prevdf_ITS_GRMO, Genus %in% get_taxa_unique(physeq_ITS_GRMO_sub, "Genus"))
ggplot(prevdf1_ITS_GRMO, aes(TotalAbundance, Prevalence / nsamples(physeq_ITS_GRMO_sub),color=Genus)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Genus) + theme(legend.position="none")
prevalenceThreshold = 0.05* nsamples(physeq_ITS_GRMO_sub)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1_ITS_GRMO)[(prevdf1_ITS_GRMO$Prevalence >= prevalenceThreshold)]
ps_ITS_GRMO = prune_taxa(keepTaxa, physeq_ITS_GRMO_sub)
# How many genera would be present after filtering?
length(get_taxa_unique(ps_ITS_GRMO, taxonomic.rank = "Genus"))


#Merge#
ps_ITS_16S_GRMO <- merge_phyloseq(ps_16S_GRMO, ps_ITS_GRMO)


###tax_glom to get to genus level###
ps_ITS_16S_GRMO_genus <- tax_glom(ps_ITS_16S_GRMO, taxrank="Genus", NArm=FALSE)
ps_ITS_16S_GRMO_genus_OTU <- as(otu_table(ps_ITS_16S_GRMO_genus), "matrix")

otu <- otu_table(ps_ITS_16S_GRMO_genus)
otu <- as.data.frame(otu)


######DESeq2 Differential Abundance######
#DESeq Analysis
dds = phyloseq_to_deseq2(ps_ITS_16S_GRMO_genus, ~ Fiber*Protein)  ####NOTE: CANNOT USE ps_g_16S_GRMO due to an error, but can use phyloseq object####
dds$Fiber <- relevel(dds$Fiber, "Low")
dds$Protein <- relevel(dds$Protein, "Low")

#Calculate geometric means
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, type="poscounts")
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts <- as.data.frame(normalized_counts)
metadata <- read.table("rodent-metadata-az.txt", header=TRUE)
normalized_counts <- t(normalized_counts)
normalized_counts <- cbind(normalized_counts, metadata[, 4:6])
write.csv(normalized_counts, "GRMO normalized 16S_ITS counts.csv")

#Run DESeq2
dds = DESeq(dds, fitType="local")
resultsNames(dds) # "Intercept" "Fiber_High_vs_Low" "Protein_High_vs_Low" "FiberHigh.ProteinHigh"
res = results(dds, name="Protein_High_vs_Low")
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_ITS_16S_GRMO_genus)[rownames(sigtab), ], "matrix"))
head(sigtab)
summary(sigtab)
summary(res)

#Shrink
shrink <- lfcShrink(dds, coef=c("Protein_High_vs_Low"), sigtab, type="ashr", lfcThreshold = 0)
shrink <- shrink[order(shrink$padj, na.last=NA), ]
sigtab1 = shrink[(shrink$padj < alpha), ]
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(ps_ITS_16S_GRMO_genus)[rownames(sigtab1), ], "matrix"))
summary(sigtab1)

write.csv(sigtab1, "GRMO DESeq Fiber.csv")
DESeq2::plotMA(shrink, ylim=c(-25,25)) #MAplot to look at 
shrink1 = as.data.frame(shrink)
shrink1 = mutate(shrink1, sig=ifelse(shrink1$padj<0.05, "FDR<0.05", "NotSig"))
shrink1[which(abs(shrink1$log2FoldChange)<1.0), "sig"] = "NotSig"
shrink1 = cbind(as(shrink1, "data.frame"), as(tax_table(ps_ITS_16S_GRMO_genus)[rownames(shrink1), ], "matrix"))


write.csv(shrink1, "GRMO_DEseq2_Interaction_SILVA_Genus.csv")


####Taxon Filtering####
#get 16S genus#
# Agglomerate to genus level
ps_g_16S_GRMO <- phyloseq::tax_glom(ps_16S_GRMO, taxrank = "Genus")
taxtab <- ps_g_16S_GRMO@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Family)
miss_f <- which(taxtab[, "Family"] == "f__")
miss_g <- which(taxtab[, "Genus"] == "g__")

# Number unspecified genera
taxtab[miss_f, "Family"] <- paste0("f__", 1:length(miss_f))
taxtab[miss_g, "Genus"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g <- which(duplicated(taxtab[, "Genus"]) |
                  duplicated(taxtab[, "Genus"], fromLast = TRUE))

for(i in seq_along(taxtab)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Order"], ")")
  } else if(i %in% miss_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Family"], ")")
  }
  
  # Family names are added to duplicate genera
  if(i %in% dupl_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Family"], ")")
  }
}

ps_g_16S_GRMO@tax_table@.Data <- taxtab
rownames(ps_g_16S_GRMO@otu_table@.Data) <- taxtab[, "Genus"]

#get ITS genus#
ps_genus_ITS_GRMO <- phyloseq::tax_glom(ps_ITS_GRMO, taxrank = "Genus")
taxtab_ITS_GRMO <- ps_genus_ITS_GRMO@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Family)
miss_f <- which(taxtab_ITS_GRMO[, "Family"] == "f__")
miss_g <- which(taxtab_ITS_GRMO[, "Genus"] == "g__")

# Number unspecified genera
taxtab_ITS_GRMO[miss_f, "Family"] <- paste0("f__", 1:length(miss_f))
taxtab_ITS_GRMO[miss_g, "Genus"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g1 <- which(duplicated(taxtab_ITS_GRMO[, "Genus"]) |
                   duplicated(taxtab_ITS_GRMO[, "Genus"], fromLast = TRUE))

for(i in seq_along(taxtab_ITS_GRMO)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab_ITS_GRMO[i, "Genus"] <- paste0(taxtab_ITS_GRMO[i, "Genus"], "(", taxtab_ITS_GRMO[i, "Order"], ")")
  } else if(i %in% miss_g){
    taxtab_ITS_GRMO[i, "Genus"] <- paste0(taxtab_ITS_GRMO[i, "Genus"], "(", taxtab_ITS_GRMO[i, "Family"], ")")
  }
  
  # Family names are added to duplicate genera
  if(i %in% dupl_g1){
    taxtab_ITS_GRMO[i, "Genus"] <- paste0(taxtab_ITS_GRMO[i, "Genus"], "(", taxtab_ITS_GRMO[i, "Family"], ")")
  }
}

ps_genus_ITS_GRMO@tax_table@.Data <- taxtab_ITS_GRMO
rownames(ps_genus_ITS_GRMO@otu_table@.Data) <- taxtab_ITS_GRMO[, "Genus"]


####merged######
ps_ITS_16S_GRMO_g <- phyloseq::tax_glom(ps_ITS_16S_GRMO, taxrank = "Genus")
taxtab_ITS_16S_GRMO_g <- ps_ITS_16S_GRMO_g@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Family)
miss_f <- which(taxtab_ITS_16S_GRMO_g[, "Family"] == "f__")
miss_g <- which(taxtab_ITS_16S_GRMO_g[, "Genus"] == "g__")

# Number unspecified genera
taxtab_ITS_16S_GRMO_g[miss_f, "Family"] <- paste0("f__", 1:length(miss_f))
taxtab_ITS_16S_GRMO_g[miss_g, "Genus"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g1 <- which(duplicated(taxtab_ITS_16S_GRMO_g[, "Genus"]) |
                   duplicated(taxtab_ITS_16S_GRMO_g[, "Genus"], fromLast = TRUE))

for(i in seq_along(taxtab_ITS_16S_GRMO_g)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab_ITS_16S_GRMO_g[i, "Genus"] <- paste0(taxtab_ITS_16S_GRMO_g[i, "Genus"], "(", taxtab_ITS_16S_GRMO_g[i, "Order"], ")")
  } else if(i %in% miss_g){
    taxtab_ITS_16S_GRMO_g[i, "Genus"] <- paste0(taxtab_ITS_16S_GRMO_g[i, "Genus"], "(", taxtab_ITS_16S_GRMO_g[i, "Family"], ")")
  }
  
  # Family names are added to duplicate genera
  if(i %in% dupl_g1){
    taxtab_ITS_16S_GRMO_g[i, "Genus"] <- paste0(taxtab_ITS_16S_GRMO_g[i, "Genus"], "(", taxtab_ITS_16S_GRMO_g[i, "Family"], ")")
  }
}

ps_ITS_16S_GRMO_g@tax_table@.Data <- taxtab_ITS_16S_GRMO_g
rownames(ps_ITS_16S_GRMO_g@otu_table@.Data) <- taxtab_ITS_16S_GRMO_g[, "Genus"]
ITS_16S_GRMO_genus_table <- as(otu_table(ps_ITS_16S_GRMO_g), "matrix")


######glmmTMB########
library(broom)
library(car)
library(magicfor)
library(glmmTMB)
library(AICcmodavg)
library(bbmle)

data <- ITS_16S_GRMO_genus_table
df_transpose <- as.data.frame(t(as.matrix(data)))
metadata <- read.table("rodent-metadata-az.txt", header=TRUE)
df_metadata <- cbind(df_transpose, metadata$Protein)
df_metadata <- cbind(df_metadata, metadata$Fiber)
df_metadata <- write.csv(df_metadata, "~/Desktop/RodentMicrobiome2024/glmm_metadata.csv")
df_metadata <- read.csv("~/Desktop/RodentMicrobiome2024/glmm_metadata.csv")
Fiber <- df_metadata$Fiber
Protein <- df_metadata$Protein


# names of lipids
microbe.names <- colnames(df_metadata)[4:132]
no.microbes <- length(microbe.names)

# create a named list to hold the fitted models
fitlistnb2 <- as.list(1:no.microbes)
names(fitlistnb2) <- microbe.names
review_nb2 <- vector("list", length=129)

# loop over lipid names
for(i in microbe.names){ 
  
  # print status
  print(paste("Hurry the fuck up:", i, "which is", which(microbe.names==i), "out of", no.microbes))
  
  # create temporary data matrix and model formula
  tmp <- df_metadata[, c(i,"Fiber","Protein")]
  fml <- as.formula(paste( i, "~", paste(c("Fiber*Protein"), collapse="+") ) )
  
  # assign fit to list by name
  fitlistnb2[[i]]<- glmmTMB(fml, family=nbinom2(link="log"), data=tmp)
  
  #print summary
  summariesnb2 <- print(summary(fitlistnb2[[i]])) 
  cofsnb2 <- coef(summary(fitlistnb2[[i]]))
  outputnb2 <- rep(do.call(rbind.data.frame, cofsnb2))
  review_nb2[[i]] <- outputnb2
}

LIST_nb2 <- do.call(rbind.data.frame, review_nb2)

write.csv(LIST_nb2, "GRMO glmmTMB Genus Summary SILVA.csv")

#########NETWORK CODE##############
library(SpiecEasi)
library(NetCoMi)
library(phyloseq)
library(qgraph)
library(igraph)
library(vegan)
library(MCL)
library(qiime2R)
library(ggplot2)
library(phyloseqCompanion)
library(CoDaSeq)
library(microbiome)
library(DESeq2)
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
library(scales)
library(dplyr)
library(grid)
library(reshape2)
library(ape)
library(ashr)
library(limma)
library(edgeR)
library(splines)
library(RColorBrewer)

setwd("~/Desktop/RodentMicrobiome2024")

#####create your phyloseq object##########
##16S##
physeq_16S_GRMO <- qza_to_phyloseq(features="table-bacteria-lab-GRMO.qza", taxonomy="taxonomy.qza", metadata="rodent-metadata-az.txt")
rank_names(physeq_16S_GRMO)
table(tax_table(physeq_16S_GRMO)[, "Genus"], exclude = NULL) #view low abundant or uncharacterized taxa
physeq_16S_GRMO_sub <- subset_taxa(physeq_16S_GRMO, !is.na(Genus) & !Genus %in% c("", "uncharacterized")) #remove any with Phyla = NA
# Compute prevalence of each feature, store as data.frame
prevdf_GRMO = apply(X = otu_table(physeq_16S_GRMO_sub),
                    MARGIN = ifelse(taxa_are_rows(physeq_16S_GRMO_sub), yes = 1, no = 2),
                    FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_GRMO = data.frame(Prevalence = prevdf_GRMO,
                         TotalAbundance = taxa_sums(physeq_16S_GRMO_sub),
                         tax_table(physeq_16S_GRMO_sub))
plyr::ddply(prevdf_GRMO, "Genus", function(df){cbind(mean(df$Prevalence),sum(df$Prevalence))}) #Shows prevalence of phyla, remove low prevalnence(?)
prevdf1_GRMO = subset(prevdf_GRMO, Genus %in% get_taxa_unique(physeq_16S_GRMO_sub, "Genus"))
ggplot(prevdf1_GRMO, aes(TotalAbundance, Prevalence / nsamples(physeq_16S_GRMO_sub),color=Genus)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Genus) + theme(legend.position="none")
prevalenceThreshold = 0.05* nsamples(physeq_16S_GRMO_sub)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1_GRMO)[(prevdf1_GRMO$Prevalence >= prevalenceThreshold)]
ps_16S_GRMO = prune_taxa(keepTaxa, physeq_16S_GRMO_sub)
# How many genera would be present after filtering?
length(get_taxa_unique(ps_16S_GRMO, taxonomic.rank = "Genus"))

# Agglomerate to genus level
ps_16S_GRMO_filter <- phyloseq::tax_glom(ps_16S_GRMO, taxrank = "Genus")
taxtab <- ps_16S_GRMO_filter@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Family)
miss_f <- which(taxtab[, "Family"] == "f__")
miss_g <- which(taxtab[, "Genus"] == "g__")

# Number unspecified genera
taxtab[miss_f, "Family"] <- paste0("f__", 1:length(miss_f))
taxtab[miss_g, "Genus"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g <- which(duplicated(taxtab[, "Genus"]) |
                  duplicated(taxtab[, "Genus"], fromLast = TRUE))

for(i in seq_along(taxtab)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Family"], ")")
  } else if(i %in% miss_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Order"], ")")
  }
  
  # Family names are added to duplicate genera
  if(i %in% dupl_g){
    taxtab[i, "Genus"] <- paste0(taxtab[i, "Genus"], "(", taxtab[i, "Order"], ")")
  }
}

ps_16S_GRMO_filter@tax_table@.Data <- taxtab
rownames(ps_16S_GRMO_filter@otu_table@.Data) <- taxtab[, "Genus"]
otu <- as.data.frame(otu_table(ps_16S_GRMO_filter))


##ITS##
physeq_ITS_GRMO <- qza_to_phyloseq(features="~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/table-sans-no-wild-fecals-GRMO-only.qza", taxonomy= "~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/taxonomy.qza", metadata="~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/rodent-metadata-az-full.txt")
rank_names(physeq_ITS_GRMO)
table(tax_table(physeq_ITS_GRMO)[, "Genus"], exclude = NULL) #view low abundant or uncharacterized taxa
physeq_ITS_GRMO_sub <- subset_taxa(physeq_ITS_GRMO, !is.na(Genus) & !Genus %in% c("", "uncharacterized")) #remove any with Phyla = NA
# Compute prevalence of each feature, store as data.frame
prevdf_ITS_GRMO = apply(X = otu_table(physeq_ITS_GRMO_sub),
                        MARGIN = ifelse(taxa_are_rows(physeq_ITS_GRMO_sub), yes = 1, no = 2),
                        FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_ITS_GRMO = data.frame(Prevalence = prevdf_ITS_GRMO,
                             TotalAbundance = taxa_sums(physeq_ITS_GRMO_sub),
                             tax_table(physeq_ITS_GRMO_sub))
plyr::ddply(prevdf_ITS_GRMO, "Genus", function(df){cbind(mean(df$Prevalence),sum(df$Prevalence))}) #Shows prevalence of phyla, remove low prevalnence(?)
prevdf1_ITS_GRMO = subset(prevdf_ITS_GRMO, Genus %in% get_taxa_unique(physeq_ITS_GRMO_sub, "Genus"))
ggplot(prevdf1_ITS_GRMO, aes(TotalAbundance, Prevalence / nsamples(physeq_ITS_GRMO_sub),color=Genus)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Genus) + theme(legend.position="none")
prevalenceThreshold = 0.05* nsamples(physeq_ITS_GRMO_sub)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1_ITS_GRMO)[(prevdf1_ITS_GRMO$Prevalence >= prevalenceThreshold)]
ps_ITS_GRMO = prune_taxa(keepTaxa, physeq_ITS_GRMO_sub)
# How many genera would be present after filtering?
length(get_taxa_unique(ps_ITS_GRMO, taxonomic.rank = "Genus"))


ps_ITS_GRMO_filter <- phyloseq::tax_glom(ps_ITS_GRMO, taxrank = "Genus")
taxtab1 <- ps_ITS_GRMO_filter@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Family)
miss_f <- which(taxtab1[, "Family"] == "f__")
miss_g <- which(taxtab1[, "Genus"] == "g__")

# Number unspecified genera
taxtab1[miss_f, "Family"] <- paste0("f__", 1:length(miss_f))
taxtab1[miss_g, "Genus"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g1 <- which(duplicated(taxtab1[, "Genus"]) |
                   duplicated(taxtab1[, "Genus"], fromLast = TRUE))

for(i in seq_along(taxtab1)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab1[i, "Genus"] <- paste0(taxtab1[i, "Genus"], "(", taxtab1[i, "Family"], ")")
  } else if(i %in% miss_g){
    taxtab1[i, "Genus"] <- paste0(taxtab1[i, "Genus"], "(", taxtab1[i, "Order"], ")")
  }
  
  # Family names are added to duplicate genera
  if(i %in% dupl_g1){
    taxtab1[i,"Genus"] <- paste0(taxtab1[i,"Genus"],"(", taxtab1[i, "Family"], ")")
  }
}

ps_ITS_GRMO_filter@tax_table@.Data <- taxtab1
rownames(ps_ITS_GRMO_filter@otu_table@.Data) <- taxtab1[, "Genus"]
otu <- as.data.frame(otu_table(ps_ITS_GRMO_filter))



counts_16S <- as.matrix(t(phyloseq::otu_table(ps_16S_GRMO_filter)@.Data))
counts_16S <- t(counts_16S)
counts_ITS <- as.matrix(t(phyloseq::otu_table(ps_ITS_GRMO_filter)@.Data))
counts_ITS <- t(counts_ITS)
write.csv(counts_16S, "~/Desktop/counts_16S_GRMO.csv")
write.csv(counts_ITS, "~/Desktop/counts_ITS_GRMO.csv")

####FIBER#####
high_fiber_16S <- counts_16S[,c(1,4,6,7,11,12,14:16,18,20,22,24,27:29,34:36)]
high_fiber_16S <- t(high_fiber_16S)
low_fiber_16S <- counts_16S[, c(2,3,5,8:10,13,17,19,21,23,25,26,30:33,37,38)]
low_fiber_16S <-t(low_fiber_16S)

high_fiber_ITS <- counts_ITS[,c(1,4,6,7,11,12,14:16,18,20,22,24,27:29,34:36)]
high_fiber_ITS <- t(high_fiber_ITS)
low_fiber_ITS <- counts_ITS[, c(2,3,5,8:10,13,17,19,21,23,25,26,30:33,37,38)]
low_fiber_ITS <- t(low_fiber_ITS)


set.seed(123456)


# Run SpiecEasi and create association matrix for group 1
spiec_result_gr1 <- multi.spiec.easi(list(high_fiber_16S, high_fiber_ITS), 
                                     method='mb', nlambda=100, 
                                     lambda.min.ratio=0.01, 
                                     pulsar.params = list(thresh = 0.05))
assoMat1 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec_result_gr1), mode = "ave")
assoMat1 <- as.matrix(assoMat1)

spiec_result_gr2 <- multi.spiec.easi(list(low_fiber_16S, low_fiber_ITS), 
                                     method='mb', nlambda=100, 
                                     lambda.min.ratio=0.01, 
                                     pulsar.params = list(thresh = 0.05))

assoMat2 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec_result_gr2), mode = "ave")

assoMat2 <- as.matrix(assoMat2)


#Get taxa names
taxnames <- c(taxa_names(ps_16S_GRMO_filter), taxa_names(ps_ITS_GRMO_filter))

colnames(assoMat1) <- rownames(assoMat1) <- taxnames
diag(assoMat1) <- 1

colnames(assoMat2) <- rownames(assoMat2) <- taxnames
diag(assoMat2) <- 1


library(NetCoMi)
##
# Network construction (pass association matrices to netConstruct)
# - sparsMethod must be set to "none" because sparsification is already included in SpiecEasi
net_fiber_16S_ITS <- netConstruct(data = assoMat1, data2 = assoMat2, 
                                  dataType = "condDependence", 
                                  sparsMethod = "none")

# Network analysis
netprops_fiber_16S_ITS <- netAnalyze(net_fiber_16S_ITS, hubPar = "eigenvector")
nodeCols <- c(rep("lightblue", ntaxa(ps_16S_GRMO_filter)), rep("orange", ntaxa(ps_ITS_GRMO_filter)))
names(nodeCols) <- taxnames

plot(netprops_fiber_16S_ITS, 
     sameLayout = FALSE, 
     layoutGroup = "union",
     shortenLabels = "none",
     nodeColor = "colorVec", 
     colorVec = nodeCols,
     nodeSize = "eigen", 
     nodeSizeSpread = 2,
     labelScale = FALSE,
     cexNodes = 4, 
     cexLabels = .2,
     cexHubLabels = .5,
     cexTitle = 3,
     rmSingles = "inboth",
     groupNames = c("High Fiber", "Low Fiber"))


netcomp_fiber_16S_ITS <- netCompare(netprops_fiber_16S_ITS, permTest = FALSE)
summary(netcomp_fiber_16S_ITS, groupNames = c("High Fiber", "Low Fiber"))

net_single1 <- netConstruct(assoMat2,
                            filtTax = "highestFreq",
                            filtTaxPar = list(totalReads = 3),
                            measure = "spieceasi",
                            normMethod = "none", 
                            zeroMethod = "none",
                            sparsMethod = "none", 
                            dissFunc = "signed",
                            verbose = 3,
                            seed = 123456)

net_single<- netConstruct(data = assoMat1,
                          dataType = "condDependence", 
                          sparsMethod = "none",
                          measure = "spieceasi")

props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

summary(props_single, numbNodes = 5L)




#####PROTEIN######
high_protein_16S <- counts_16S[,c(1,2,4,6,11,13,17,19,24,26:29,31:34,37)]
high_protein_16S <- t(high_protein_16S)
low_protein_16S <- counts_16S[, c(3,5,7:10,12,14:16,18,20:23,25,30,35,36,38)]
low_protein_16S <-t(low_protein_16S)

high_protein_ITS <- counts_ITS[,c(1,2,4,6,11,13,17,19,24,26:29,31:34,37)]
high_protein_ITS <- t(high_protein_ITS)
low_protein_ITS <- counts_ITS[, c(3,5,7:10,12,14:16,18,20:23,25,30,35,36,38)]
low_protein_ITS <- t(low_protein_ITS)


set.seed(123456)


# Run SpiecEasi and create association matrix for group 1
spiec_result_gr1 <- multi.spiec.easi(list(high_protein_16S, high_protein_ITS), 
                                     method='mb', nlambda=100, 
                                     lambda.min.ratio=1e-2, 
                                     pulsar.params = list(thresh = 0.05))
assoMat1 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec_result_gr1), mode = "ave")
assoMat1 <- as.matrix(assoMat1)

spiec_result_gr2 <- multi.spiec.easi(list(low_protein_16S, low_protein_ITS), 
                                     method='mb', nlambda=100, 
                                     lambda.min.ratio=1e-2, 
                                     pulsar.params = list(thresh = 0.05))

assoMat2 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec_result_gr2), mode = "ave")

assoMat2 <- as.matrix(assoMat2)


#Get taxa names
taxnames <- c(taxa_names(ps_16S_GRMO_filter), taxa_names(ps_ITS_GRMO_filter))

colnames(assoMat1) <- rownames(assoMat1) <- taxnames
diag(assoMat1) <- 1

colnames(assoMat2) <- rownames(assoMat2) <- taxnames
diag(assoMat2) <- 1


library(NetCoMi)
##
# Network construction (pass association matrices to netConstruct)
# - sparsMethod must be set to "none" because sparsification is already included in SpiecEasi
net_protein_16S_ITS <- netConstruct(data = assoMat1, data2 = assoMat2, 
                                    dataType = "condDependence", 
                                    sparsMethod = "none")

# Network analysis
netprops_protein_16S_ITS <- netAnalyze(net_protein_16S_ITS, hubPar = "eigenvector")
nodeCols <- c(rep("lightblue", ntaxa(ps_16S_GRMO_filter)), rep("orange", ntaxa(ps_ITS_GRMO_filter)))
names(nodeCols) <- taxnames


plot(netprops_protein_16S_ITS, 
     sameLayout = TRUE, 
     layoutGroup = "union",
     shortenLabels = "none",
     nodeColor = "colorVec", 
     colorVec = nodeCols,
     nodeSize = "eigen", 
     nodeSizeSpread = 2,
     labelScale = FALSE,
     cexNodes = 4, 
     cexLabels = 0,
     cexHubLabels = .5,
     cexTitle = 3,
     groupNames = c("High Protein", "Low Protein"))


netcomp_protein_16S_ITS <- netCompare(netprops_protein_16S_ITS, permTest = FALSE)
summary(netcomp_protein_16S_ITS, groupNames = c("High Protein", "Low Protein"))


net_single<- netConstruct(data = assoMat2,
                          dataType = "condDependence", 
                          sparsMethod = "none",
                          measure = "spieceasi")

props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

summary(props_single, numbNodes = 5L)