library(qiime2R)
library(ape)
library(Biostrings)
library(biomformat)
library(phyloseq)
library(Hmisc)
library(yaml)
library(tidyr)
library(dplyr)
library(stats)
library(utils)
library(gplots)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(glmmTMB)
library(MicEco)
library(cowplot)
library(ggpubr)


setwd("~/Desktop/RodentMicrobiome2024")

ps <- qza_to_phyloseq(features="table-sans-archaea-mito-chloro-no-wild-fecals.qza", taxonomy="taxonomy.qza", metadata="rodent-metadata.txt")
dist = phyloseq::distance(ps, method="jaccard", binary = TRUE)
ordination = ordinate(ps, method="NMDS", distance=dist)
plot_ordination(ps, ordination, color="Species") + 
  theme_classic() +
  theme(strip.background = element_blank())
ps2 <- subset_samples(ps, SampleType == "Fecal")
GP.ord <- ordinate(ps2, "NMDS", "bray")
p2 = plot_ordination(ps2, GP.ord, type="samples", color="Species", shape="Diet") 
p2 + geom_polygon(aes(fill=SampleType)) + ggtitle("samples")
p2 + theme(plot.background = element_rect(fill = 'white'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill='white'), axis.line = element_line(colour = "black")) + geom_point(size = 2) + theme(legend.position='Species') + scale_color_manual(values=c("coral", "dark grey", "dodger blue")) + theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), axis.title.x=element_text(colour="black"), axis.title.y=element_text(colour="black")) 

#For all samples
SVs <- read_qza("table.qza")
taxonomy <- read_qza("taxonomy.qza")
head(taxonomy)
taxonomy <- parse_taxonomy(taxonomy$data)
head(taxonomy)

#For MOVO
metadata <- read.delim("rodent-metadata-ut.txt")

#MOVO alpha diversity
shannon <-read_qza("~/Desktop/RodentMicrobiome2024/core-metrics-results-MOVO-SILVA/shannon_vector.qza")
shannon <- shannon$data %>% rownames_to_column("SampleID")

gplots::venn(list(metadata=metadata$SampleID, shannon=shannon$SampleID))

metadata <- metadata %>% left_join(shannon)
head(metadata)
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels =c("Low", "High"))

metadata

p <- ggplot(metadata, aes(x=Fiber, y=shannon_entropy, fill =Protein)) + geom_boxplot(trim=FALSE)
MOVO_shannon <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Shannon Diversity")

MOVO_shannon



faith <- read_qza("~/Desktop/RodentMicrobiome2024/core-metrics-results-MOVO-SILVA/faith_pd_vector.qza")
faith <- faith$data %>% rownames_to_column("SampleID")
faith <- faith[, 2:3]
colnames(faith) <- c("SampleID", "faith_pd")
metadata <- metadata %>% left_join(faith)
head(metadata)

metadata$Protein <- ordered(metadata$Protein, levels =c("Low", "High"))
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))


p1 <- ggplot(metadata, aes(x=Fiber, y=faith_pd, fill =Protein)) + geom_boxplot(trim=FALSE)
MOVO_faith <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +  geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548'))  + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Faith's Phylogenetic Diversity")

MOVO_faith


  
ASVs <- read_qza("~/Desktop/RodentMicrobiome2024/core-metrics-results-MOVO-SILVA/observed_features_vector.qza")
ASVs<- ASVs$data %>% rownames_to_column("SampleID")
metadata <- metadata %>% left_join(ASVs)
head(metadata)
metadata$Protein <- ordered(metadata$Protein, levels =c("Low", "High"))
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))


p2 <- ggplot(metadata, aes(x=Fiber, y=observed_features, fill =Protein)) + geom_boxplot(trim=FALSE)
p2
MOVO_asv <-p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +  geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548'))  + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Observed ASVs") 
MOVO_asv


write.csv(metadata, "~/Desktop/MOVO_alpha_metrics.csv")


anova <- aov(shannon_entropy ~ Fiber + Protein, metadata)
summary(anova)
TukeyHSD(anova)


model <- glmmTMB(metadata$shannon_entropy ~ metadata$Fiber*metadata$Protein, family = gaussian)
summary(anova)
summary(model)



#Fiber is significant .0143
  
  
anova1 <- aov(faith_pd ~ Fiber*Protein, data = metadata)
model1 <- glmmTMB(metadata$faith_pd ~ metadata$Fiber*metadata$Protein, family = gaussian)
summary(anova1) 
TukeyHSD(anova1)
summary(model1)
#Interaction is significant (.0268)


anova2 <- aov(observed_features ~ Fiber*Protein, data=metadata)
summary(anova2)
model2 <- glmmTMB(metadata$observed_features ~ metadata$Fiber*metadata$Protein, family=gaussian)
summary(model2)
# Fiber and interaction are significant
#              Df Sum Sq Mean Sq F value  Pr(>F)   
#Fiber          1  18232   18232  10.030 0.00333 **
#Protein        1   3169    3169   1.743 0.19484   
#Fiber:Protein  1  11637   11637   6.402 0.01517 * 
#Residuals     37  67257    1818      

#MOVO beta diversity
bray_curtis<-read_qza("~/Desktop/RodentMicrobiome2024/core-metrics-results-MOVO-SILVA/bray_curtis_pcoa_results.qza")

bray_curtis2 <- bray_curtis$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata)

 
MOVO_BC <- ggplot(data=bray_curtis2, aes(x=PC1, y=PC2, color=`Protein`, shape=`Fiber`,size=1)) +
  geom_point(alpha=1) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() + scale_color_manual(values=c("Low" = "#FACBC4", "High"="#A65548")) + 
  theme(axis.title=element_text(size=18), text=element_text(family="serif"), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) + 
  xlab("PC1") +
  ylab("PC2") +
  scale_shape_manual(values=c(15, 16 ))
MOVO_BC

bcmatrix <- read.delim("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/Microbiome_9-28-21/exported-BC-matrix-MOVO-SILVA-9-28-21/distance-matrix.tsv")
bcadonis <- adonis2(bcmatrix ~ metadata$Protein*metadata$Fiber, method='bray', perm=999)
print(bcadonis)
adonis_OmegaSq(bcadonis, partial = TRUE)
#all significant (interaction .021)



#For GRMO
metadata1 <- read.delim("rodent-metadata-az.txt")

#GRMO alpha diversity
shannon <-read_qza("~/Desktop/RodentMicrobiome2024/core-metrics-results-GRMO-SILVA/shannon_vector.qza")
shannon <- shannon$data %>% rownames_to_column("SampleID")

gplots::venn(list(metadata=metadata1$SampleID, shannon=shannon$SampleID))

metadata <- metadata1 %>% left_join(shannon)
head(metadata)
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels = c("Low", "High"))

p <- ggplot(metadata, aes(x=Fiber, y=shannon_entropy, fill =Protein)) +geom_boxplot(trim=FALSE)
GRMO_shannon <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +  geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548'))  + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Shannon Diversity")
GRMO_shannon


faith <- read_qza("~/Desktop/RodentMicrobiome2024/core-metrics-results-GRMO-SILVA/faith_pd_vector.qza")
faith <- faith$data %>% rownames_to_column("SampleID")
faith <- faith[, 2:3]
colnames(faith) <- c("SampleID", "faith_pd")
metadata <- metadata %>% left_join(faith)
head(metadata)
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels = c("Low", "High"))

p <- ggplot(metadata, aes(x=Fiber, y=faith_pd, fill =Protein)) + geom_boxplot(trim=FALSE)
GRMO_faith <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +  geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548'))  + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Faith's Phylogenetic Diversity")
GRMO_faith

ASVs <- read_qza("~/Desktop/RodentMicrobiome2024/core-metrics-results-GRMO-SILVA/observed_features_vector.qza")
ASVs<- ASVs$data %>% rownames_to_column("SampleID")
metadata <- metadata %>% left_join(ASVs)
head(metadata)
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels = c("Low", "High"))

p <- ggplot(metadata, aes(x=Fiber, y=observed_features, fill =Protein)) + geom_boxplot(trim=FALSE)
GRMO_asv <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +  geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548'))  + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Observed ASVs")
GRMO_asv
write.csv(metadata, "~/Desktop/GRMO_alpha_metrics.csv")

anova <- aov(shannon_entropy ~ Protein + Fiber, data=metadata)
summary(anova)

model <- glmmTMB(shannon_entropy ~ Fiber*Protein, data=metadata, family=gaussian)
summary(model)
#Not significant    

anova1 <- aov(faith_pd ~ Protein*Fiber, data = metadata)
summary(anova1) 
model1 <- glmmTMB(faith_pd ~ Fiber*Protein, data=metadata, family = gaussian)
summary(model1)
#Not significant    

anova2 <- aov(observed_features ~ Protein*Fiber, data=metadata)
summary(anova2)
model2 <- glmmTMB(observed_features ~ Fiber*Protein, data = metadata, family=gaussian)
summary(model2)
# Not significant     

#GRMO beta diversity
bray_curtis<-read_qza("~/Desktop/RodentMicrobiome2024/core-metrics-results-GRMO-SILVA/bray_curtis_pcoa_results.qza")

bray_curtis2 <- bray_curtis$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata1)

GRMO_BC <- ggplot(data=bray_curtis2, aes(x=PC1, y=PC2, color=`Protein`, shape=`Fiber`,size=1)) +
  geom_point(alpha=1) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() + scale_color_manual(values=c("Low" = "#FACBC4", "High"="#A65548")) + 
  theme(axis.title=element_text(size=18), text=element_text(family="serif"), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) + 
  xlab("PC1") +
  ylab("PC2") +
  scale_shape_manual(values=c(15, 16 ))
GRMO_BC




bcmatrix <- read.delim("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/Microbiome_9-28-21/exported-BC-matrix-GRMO-SILVA-9-28-21/distance-matrix.tsv")
bcadonis <- adonis2(bcmatrix ~ metadata1$Fiber*metadata1$Protein, method='bray', perm=999)
print(bcadonis)
#all significant
#                                Df SumOfSqs      R2      F Pr(>F)    
#metadata$Protein                 1   0.5207 0.04304 1.6767  0.004 ** 
#metadata$Fiber                   1   0.5969 0.04934 1.9222  0.003 ***
#metadata$Protein:metadata$Fiber  1   0.4224 0.03492 1.3603  0.039 *  
#Residual                        34  10.5584 0.87271                  
#Total                           37  12.0984 1.00000  
adonis_OmegaSq(bcadonis, partial = TRUE)


#For WFMO
metadata2 <- read.delim("rodent-metadata-ky.txt")

#WFMO alpha diversity
shannon <-read_qza("~/Desktop/RodentMicrobiome2024/core-metrics-results-WFMO-SILVA/shannon_vector.qza")
shannon <- shannon$data %>% rownames_to_column("SampleID")

gplots::venn(list(metadata=metadata2$SampleID, shannon=shannon$SampleID))

metadata <- metadata2 %>% left_join(shannon)
head(metadata)
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels = c("Low", "High"))

p <- ggplot(metadata, aes(x=Fiber, y=shannon_entropy, fill =Protein)) + geom_boxplot(trim=FALSE)
WFMO_shannon <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +  geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548'))  + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Shannon Diversity")



faith <- read_qza("~/Desktop/RodentMicrobiome2024/core-metrics-results-WFMO-SILVA/faith_pd_vector.qza")
faith <- faith$data %>% rownames_to_column("SampleID")
faith <- faith[, 2:3]
colnames(faith) <- c("SampleID", "faith_pd")
metadata <- metadata %>% left_join(faith)
head(metadata)
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels = c("Low", "High"))

p <- ggplot(metadata, aes(x=Fiber, y=faith_pd, fill =Protein)) + geom_boxplot()
WFMO_faith <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +  geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548'))  + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Faith's Phylogenetic Diversity")



ASVs <- read_qza("~/Desktop/RodentMicrobiome2024/core-metrics-results-WFMO-SILVA/observed_features_vector.qza")
ASVs<- ASVs$data %>% rownames_to_column("SampleID")
metadata <- metadata %>% left_join(ASVs)
head(metadata)
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels = c("Low", "High"))

p <- ggplot(metadata, aes(x=Fiber, y=observed_features, fill =Protein)) + geom_boxplot()
WFMO_asv <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +  geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Observed ASVs")

write.csv(metadata, "~/Desktop/WFMO_alpha_metrics.csv")

anova <- aov(shannon_entropy ~ Protein*Fiber, data=metadata)
summary(anova)
TukeyHSD(anova)
model <- glmmTMB(shannon_entropy ~ Fiber*Protein, data=metadata, family=gaussian)
summary(model)
#Interaction is significant(.0124)

anova1 <- aov(faith_pd ~ Protein*Fiber, data = metadata)
summary(anova1) 
TukeyHSD(anova1)
model1 <- glmmTMB(faith_pd ~ Fiber*Protein, data=metadata, family=gaussian)
summary(model1)
#Interaction significant (.0116)


anova2 <- aov(observed_features ~ Protein*Fiber, data=metadata)
summary(anova2)
TukeyHSD(anova2)
model2 <- glmmTMB(metadata$observed_features ~ Fiber*Protein, data=metadata, family=gaussian)
summary(model2)
# Protein is significant (.0227)

#WFMO beta diversity
bray_curtis<-read_qza("~/Desktop/RodentMicrobiome2024/core-metrics-results-WFMO-SILVA/bray_curtis_pcoa_results.qza")

bray_curtis2 <- bray_curtis$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata2)

WFMO_BC <-ggplot(data=bray_curtis2, aes(x=PC1, y=PC2, color=`Protein`, shape=`Fiber`,size=1)) +
  geom_point(alpha=1) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() + scale_color_manual(values=c("Low" = "#FACBC4", "High"="#A65548")) + 
  theme(axis.title=element_text(size=18), text=element_text(family="serif"), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) + 
  xlab("PC1") +
  ylab("PC2") +
  scale_shape_manual(values=c(15, 16 ))
WFMO_BC

bcmatrix <- read.delim("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/Microbiome_9-28-21/exported-BC-matrix-WFMO-SILVA-9-28-21/distance-matrix.tsv")
bcadonis <- adonis2(bcmatrix ~ metadata2$Fiber*metadata2$Protein, method='bray', perm=999)
print(bcadonis)
#Protein is significant (.001)
adonis_OmegaSq(bcadonis, partial = TRUE)

###Group Figures#####
library(ggpubr)
alpha_16S <- ggarrange(MOVO_shannon, WFMO_shannon, GRMO_shannon, MOVO_faith, WFMO_faith, GRMO_faith, MOVO_asv, WFMO_asv, GRMO_asv, 
                    common.legend=TRUE, legend = "bottom", ncol = 3, nrow = 3)
alpha_16S

beta_16S <- ggarrange(MOVO_BC, WFMO_BC, GRMO_BC, 
                    common.legend=TRUE, legend = "bottom", ncol = 3, nrow = 1)
beta_16S











####ITS####
setwd("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21")

ps <- qza_to_phyloseq(features="table-sans-no-wild-fecals.qza", taxonomy="taxonomy.qza", metadata="rodent-metadata_subset.txt")
dist = phyloseq::distance(ps2, method="jaccard", binary = TRUE)
ordination = ordinate(ps2, method="NMDS", distance=dist)
p1 <- plot_ordination(ps, ordination, color="Species") + 
  theme_classic() +
  theme(strip.background = element_blank())
p1 + geom_polygon(aes(fill=SampleType)) + ggtitle("samples")
p1 + theme(plot.background = element_rect(fill = 'white'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill='white'), axis.line = element_line(colour = "black")) + geom_point(size=5) + theme(legend.position='none') + scale_color_manual(values=c("coral", "dark grey", "dodger blue")) + theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), axis.title.x=element_text(colour="black"), axis.title.y=element_text(colour="black")) 
p1 + theme(plot.background = element_rect(fill = 'white'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill='white'), axis.line = element_line(colour = "black")) + geom_point(size = 5) + theme(legend.position='Species') + scale_color_manual(values=c("coral", "dark grey", "dodger blue")) + theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), axis.title.x=element_text(colour="black"), axis.title.y=element_text(colour="black")) 


ps2 <- subset_samples(ps, SampleType == "Fecal")
GP.ord <- ordinate(ps2, "NMDS", "bray")
p2 = plot_ordination(ps2, GP.ord, type="samples", color="Diet", shape="SampleType") 
p2 + geom_polygon(aes(fill=SampleType)) + ggtitle("samples")
p2 + theme(plot.background = element_rect(fill = 'white'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill='white'), axis.line = element_line(colour = "black")) + geom_point(size=5) + theme(legend.position='none') + scale_color_manual(values=c("coral", "dark grey", "dodger blue")) + theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), axis.title.x=element_text(colour="black"), axis.title.y=element_text(colour="black")) 
p2 + theme(plot.background = element_rect(fill = 'white'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill='white'), axis.line = element_line(colour = "black")) + geom_point(size = 5) + theme(legend.position='Species') + scale_color_manual(values=c("coral", "dark grey", "dodger blue")) + theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), axis.title.x=element_text(colour="black"), axis.title.y=element_text(colour="black")) 




#For all samples
SVs <- read_qza("table.qza")
taxonomy <- read_qza("taxonomy.qza")
head(taxonomy)
taxonomy <- parse_taxonomy(taxonomy$data)
head(taxonomy)

#For MOVO
metadata <- read.delim("rodent-metadata-ut_ITS.txt")

#MOVO alpha diversity
shannon <-read_qza("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/core-metrics-results-MOVO-only-9-30-21/shannon_vector.qza")
shannon <- shannon$data %>% rownames_to_column("SampleID")

gplots::venn(list(metadata=metadata$SampleID, shannon=shannon$SampleID))

metadata <- metadata %>% left_join(shannon)
head(metadata)
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels =c("Low", "High"))


p <- ggplot(metadata, aes(x=Protein, y=shannon_entropy, fill =Fiber)) + geom_boxplot(trim=FALSE)
MOVO_shannon_ITS <-  p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                              panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Shannon Diversity")
MOVO_shannon_ITS

faith <- read_qza("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/core-metrics-results-MOVO-only-9-30-21/faith_pd_vector.qza")
faith <- faith$data %>% rownames_to_column("SampleID")
metadata <- metadata %>% left_join(faith)
head(metadata)
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels =c("Low", "High"))

metadata
p <- ggplot(metadata, aes(x=Fiber, y=faith_pd, fill =Protein)) + geom_boxplot(trim=FALSE)
MOVO_FPD_ITS <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Faith's Phylogenetic Diversity")

MOVO_FPD_ITS

ASVs <- read_qza("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/core-metrics-results-MOVO-only-9-30-21/observed_features_vector.qza")
ASVs<- ASVs$data %>% rownames_to_column("SampleID")
metadata <- metadata %>% left_join(ASVs)
head(metadata)
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels =c("Low", "High"))

metadata
p <- ggplot(metadata, aes(x=Fiber, y=observed_features, fill =Protein)) + geom_boxplot(trim=FALSE)
MOVO_Richness_ITS <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Observed ASVs")

MOVO_Richness_ITS


anova <- aov(metadata$shannon_entropy ~ metadata$Fiber + metadata$Protein)
summary(anova)
TukeyHSD(anova)
model <- glmmTMB(shannon_entropy ~ Fiber*Protein, data=metadata, family=gaussian)
summary(model)

anova1 <- aov(faith_pd ~ Fiber + Protein, data = metadata)
summary(anova1) 
TukeyHSD(anova1)

model1 <- glmmTMB(faith_pd ~ Fiber + Protein, data = metadata, family = gaussian)
summary(model1)

anova2 <- aov(observed_features ~ Fiber*Protein, data=metadata)
summary(anova2)
TukeyHSD(anova2)
model2 <- glmmTMB(observed_features ~ Fiber*Protein, data = metadata, family = gaussian)
summary(model2)

#MOVO beta diversity
bray_curtis<-read_qza("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/core-metrics-results-MOVO-only-9-30-21/bray_curtis_pcoa_results.qza")
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels =c("Low", "High"))
bray_curtis2 <- bray_curtis$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata)


MOVO_BC_ITS <- ggplot(data=bray_curtis2, aes(x=PC1, y=PC2, color=`Protein`, shape=`Fiber`,size=1)) +
  geom_point(alpha=1) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() + scale_color_manual(values=c("Low" = "#FACBC4", "High"="#A65548")) + 
  theme(axis.title=element_text(size=18), text=element_text(family="serif"), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) + 
  xlab("PC1") +
  ylab("PC2") +
  scale_shape_manual(values=c(15, 16 ))
MOVO_BC_ITS


bcmatrix <- read.delim("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/core-metrics-results-MOVO-only-9-30-21/distance-matrix.tsv")
bcmatrix <- bcmatrix[, -1]
bcadonis <- adonis2(bcmatrix ~ metadata$Protein*metadata$Fiber, method='bray', perm=999)
print(bcadonis)
#all significant (interaction .001)
adonis_OmegaSq(bcadonis, partial=TRUE)


#For GRMO
metadata1 <- read.delim("rodent-metadata-az_ITS.txt")

#GRMO alpha diversity
shannon <-read_qza("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome//ITS_9-30-21/core-metrics-results-GRMO-only-9-30-21/shannon_vector.qza")
shannon <- shannon$data %>% rownames_to_column("SampleID")

gplots::venn(list(metadata=metadata1$SampleID, shannon=shannon$SampleID))

metadata <- metadata1 %>% left_join(shannon)
head(metadata)
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels =c("Low", "High"))
p <- ggplot(metadata, aes(x=Fiber, y=shannon_entropy, fill =Protein)) + geom_boxplot(trim=FALSE)
GRMO_shannon_ITS <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Shannon Diversity")

GRMO_shannon_ITS

faith <- read_qza("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/core-metrics-results-GRMO-only-9-30-21/faith_pd_vector.qza")
faith <- faith$data %>% rownames_to_column("SampleID")
metadata <- metadata %>% left_join(faith)
head(metadata)
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels =c("Low", "High"))
p <- ggplot(metadata, aes(x=Fiber, y=faith_pd, fill =Protein)) + geom_boxplot(trim=FALSE)
GRMO_FPD_ITS <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Faith's Phylogenetic Diversity")

GRMO_FPD_ITS


ASVs <- read_qza("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/core-metrics-results-GRMO-only-9-30-21/observed_features_vector.qza")
ASVs<- ASVs$data %>% rownames_to_column("SampleID")
metadata <- metadata %>% left_join(ASVs)
head(metadata)
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels =c("Low", "High"))
p <- ggplot(metadata, aes(x=Fiber, y=observed_features, fill =Protein)) + geom_boxplot(trim=FALSE)
GRMO_ASV_ITS <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Observed ASVs")

GRMO_ASV_ITS


anova <- aov(shannon_entropy ~ Protein*Fiber, data=metadata)
summary(anova)
TukeyHSD(anova)
model <- glmmTMB(shannon_entropy ~ Fiber*Protein, data = metadata, family = gaussian)
summary(model)
#Not significant    

anova1 <- aov(faith_pd ~ Protein*Fiber, data = metadata)
summary(anova1) 
TukeyHSD(anova1)
model1 <- glmmTMB(faith_pd ~ Fiber + Protein, data = metadata, family = gaussian)
summary(model1)
#Not significant    

anova2 <- aov(observed_features ~ Protein*Fiber, data=metadata)
summary(anova2)
model2 <- glmmTMB(observed_features ~ Fiber*Protein, data = metadata, family = gaussian)
summary(model2)
# Not significant     

#GRMO beta diversity
bray_curtis<-read_qza("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/core-metrics-results-GRMO-only-9-30-21/bray_curtis_pcoa_results.qza")

bray_curtis2 <- bray_curtis$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata)

GRMO_BC_ITS <- ggplot(data=bray_curtis2, aes(x=PC1, y=PC2, color=`Protein`, shape=`Fiber`,size=1)) +
  geom_point(alpha=1) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() + scale_color_manual(values=c("Low" = "#FACBC4", "High"="#A65548")) + 
  theme(axis.title=element_text(size=18), text=element_text(family="serif"), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) + 
  xlab("PC1") +
  ylab("PC2") +
  scale_shape_manual(values=c(15, 16 ))
GRMO_BC_ITS

bcmatrix <- read.delim("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/exported-BC-matrix-GRMO-9-30-21/distance-matrix.tsv")
bcadonis <- adonis2(bcmatrix ~ metadata1$Protein*metadata1$Fiber, method='bray', perm=999)
print(bcadonis)
adonis_OmegaSq(bcadonis, partial=TRUE)
#interaction significant (.001)


#For WFMO
metadata2 <- read.delim("rodent-metadata-ky_ITS.txt")

#WFMO alpha diversity
shannon <-read_qza("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/core-metrics-results-WFMO-only-9-30-21/shannon_vector.qza")
shannon <- shannon$data %>% rownames_to_column("SampleID")

gplots::venn(list(metadata=metadata2$SampleID, shannon=shannon$SampleID))

metadata <- metadata2 %>% left_join(shannon)
head(metadata)
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels =c("Low", "High"))
p <- ggplot(metadata, aes(x=Fiber, y=shannon_entropy, fill =Protein)) + geom_boxplot(trim=FALSE)
WFMO_shannon_ITS <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Shannon Diversity")
WFMO_shannon_ITS


faith <- read_qza("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/core-metrics-results-WFMO-only-9-30-21/faith_pd_vector.qza")
faith <- faith$data %>% rownames_to_column("SampleID")
metadata <- metadata %>% left_join(faith)
head(metadata)
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels =c("Low", "High"))
p <- ggplot(metadata, aes(x=Fiber, y=faith_pd, fill =Protein)) + geom_boxplot(trim=FALSE)
WFMO_FPD_ITS <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Faith's Phylogenetic Diversity")

WFMO_FPD_ITS


ASVs <- read_qza("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/core-metrics-results-WFMO-only-9-30-21/observed_features_vector.qza")
ASVs<- ASVs$data %>% rownames_to_column("SampleID")
metadata <- metadata %>% left_join(ASVs)
head(metadata)
metadata$Fiber <- ordered(metadata$Fiber, levels = c("Low", "High"))
metadata$Protein <- ordered(metadata$Protein, levels =c("Low", "High"))
p <- ggplot(metadata, aes(x=Fiber, y=observed_features, fill =Protein)) + geom_boxplot(trim=FALSE)
WFMO_ASV_ITS <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c('#FACBC4','#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Observed ASVs")

WFMO_ASV_ITS


anova <- aov(shannon_entropy ~ Protein*Fiber, data=metadata)
summary(anova)
model <- glmmTMB(shannon_entropy ~ Fiber*Protein, data = metadata, family = gaussian)
summary(model)
#Protein is significant (.0331)

anova1 <- aov(faith_pd ~ Protein*Fiber, data = metadata)
summary(anova1) 
model1 <- glmmTMB(faith_pd ~ Fiber*Protein, data = metadata, family = gaussian)
summary(model1)
#Not significant


anova2 <- aov(observed_features ~ Protein*Fiber, data=metadata)
summary(anova2)
model2 <- glmmTMB(observed_features ~ Fiber*Protein, data = metadata, family = gaussian)
summary(model2)
#Not significant

#WFMO beta diversity
bray_curtis<-read_qza("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/core-metrics-results-WFMO-only-9-30-21/bray_curtis_pcoa_results.qza")

bray_curtis2 <- bray_curtis$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata)

WFMO_BC_ITS <- ggplot(data=bray_curtis2, aes(x=PC1, y=PC2, color=`Protein`, shape=`Fiber`,size=1)) +
  geom_point(alpha=1) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() + scale_color_manual(values=c("Low" = "#FACBC4", "High"="#A65548")) + 
  theme(axis.title=element_text(size=18), text=element_text(family="serif"), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) + 
  xlab("PC1") +
  ylab("PC2") +
  scale_shape_manual(values=c(15, 16 ))
WFMO_BC_ITS

bcmatrix <- read.delim("~/Desktop/NB Mac/Research/KohlLab/Kohl Lab Projects/Rodent Microbiome/ITS_9-30-21/core-metrics-results-WFMO-only-9-30-21/distance-matrix.tsv")
bcmatrix <- bcmatrix[,-1]
bcadonis <- adonis2(bcmatrix ~ metadata2$Fiber*metadata2$Protein, method='bray', perm=999)
print(bcadonis)
adonis_OmegaSq(bcadonis, partial=TRUE)
#Protein is significant (.027)

####Group Figure####
library(ggpubr)
alpha_ITS <- ggarrange(MOVO_shannon_ITS, WFMO_shannon_ITS, GRMO_shannon_ITS, MOVO_FPD_ITS, WFMO_FPD_ITS, GRMO_FPD_ITS, MOVO_Richness_ITS, WFMO_ASV_ITS, GRMO_ASV_ITS, 
                    common.legend=TRUE, legend = "bottom", ncol = 3, nrow = 3)
alpha_ITS

beta_ITS <- ggarrange(MOVO_BC_ITS, WFMO_BC_ITS, GRMO_BC_ITS, 
                      common.legend=TRUE, legend = "bottom", ncol = 3, nrow = 1)
beta_ITS





####Cecum Plasticity#####
data <- read.csv("~/Desktop/Rodent Microbiome_Phenotype.csv")
species <- data$Species
cecum_mass <- data$EmptyCecumMass..g.
data$cecum_body_ratio <- data$CecumContents.g./data$BodyMass..g.
fiber <- ordered(data$Fiber, levels = c("Low", "High"))
protein <- ordered(data$Protein, levels =c("Low", "High"))

grmo <- data[1:38,]
grmo$Fiber <- ordered(grmo$Fiber, levels = c("Low", "High"))
grmo$Protein <- ordered(grmo$Protein, levels =c("Low", "High"))
movo <- data[39:79,]
movo$Fiber <- ordered(movo$Fiber, levels = c("Low", "High"))
movo$Protein <- ordered(movo$Protein, levels =c("Low", "High"))
wfmo <- data[80:119,]
wfmo$Fiber <- ordered(wfmo$Fiber, levels = c("Low", "High"))
wfmo$Protein <- ordered(wfmo$Protein, levels =c("Low", "High"))

grmo_cecum <- glm(grmo$EmptyCecumMass..g. ~ grmo$Fiber*grmo$Protein)
summary(grmo_cecum)

grmo_ratio <- grmo[-34,]
grmo_ratio_model <- glm(grmo_ratio$cecum_body_ratio ~ grmo_ratio$Fiber*grmo_ratio$Protein)
summary(grmo_ratio_model)

movo_cecum <- glm(movo$EmptyCecumMass..g. ~ movo$Fiber + movo$Protein)
summary(movo_cecum)

movo_cecum_ratio <- glm(movo$cecum_body_ratio ~ movo$Fiber * movo$Protein)
summary(movo_cecum_ratio)

wfmo_cecum <- glm(wfmo$EmptyCecumMass..g. ~ wfmo$Fiber + wfmo$Protein, data=wfmo)
summary(wfmo_cecum)

wfmo_cecum_ratio <- glm(wfmo$cecum_body_ratio ~ wfmo$Fiber + wfmo$Protein)
summary(wfmo_cecum_ratio)

p <- ggplot(grmo, aes(x=Fiber, y=EmptyCecumMass..g., fill =Protein)) + geom_boxplot(trim=FALSE)
p
p_grmo_cecum <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c("Low" = '#FACBC4', "High"='#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Empty Cecum Mass (g)")
p_grmo_cecum

p <- ggplot(wfmo, aes(x=Fiber, y=wfmo$EmptyCecumMass..g., fill =Protein)) + geom_boxplot(trim=FALSE)
p
p_wfmo_cecum <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c("Low" = '#FACBC4', "High"='#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Empty Cecum Mass (g)")
p_wfmo_cecum

p <- ggplot(movo, aes(x=Fiber, y=movo$EmptyCecumMass..g., fill =Protein)) + geom_boxplot(trim=FALSE)
p
p_movo_cecum <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c("Low" = '#FACBC4', "High"='#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Empty Cecum Mass (g)")
p_movo_cecum


p <- ggplot(grmo, aes(x=Fiber, y=cecum_body_ratio, fill =Protein)) + geom_boxplot(trim=FALSE)
p
p_grmo_ratio <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c("Low" = '#FACBC4',"High" = '#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Cecum Contents/Body Mass (g)")
p_grmo_ratio


p <- ggplot(wfmo, aes(x=Fiber, y=wfmo$cecum_body_ratio, fill =Protein)) + geom_boxplot(trim=FALSE)
p
p_wfmo_ratio <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c("Low" = '#FACBC4',"High" = '#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Cecum Contents/Body Mass (g)")
p_wfmo_ratio


p <- ggplot(movo, aes(x=Fiber, y=movo$cecum_body_ratio, fill =Protein)) + geom_boxplot(trim=FALSE)
p
p_movo_ratio <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_dotplot(binaxis='y', stackdir = 'center', position=position_dodge(.75)) + scale_fill_manual(values=c("Low" = '#FACBC4',"High" = '#A65548')) + scale_x_discrete(limits=c("Low","High"), labels=c("Low", "High")) + labs(x="Fiber", y="Cecum Contents/Body Mass (g)")
p_movo_ratio

####Group Figure####
library(ggpubr)
cecum_all <- ggarrange(p_movo_cecum, p_wfmo_cecum, p_grmo_cecum, 
                       common.legend=TRUE, legend = "bottom", ncol = 3, nrow = 1)
cecum_all

cecum_ratio_all <- ggarrange(p_movo_ratio, p_wfmo_ratio, p_grmo_ratio, 
                                         common.legend=TRUE, legend = "bottom", ncol = 3, nrow = 1)
cecum_ratio_all

cecum_plots <- ggarrange(p_movo_cecum, p_wfmo_cecum, p_grmo_cecum, p_movo_ratio, p_wfmo_ratio, p_grmo_ratio,
                         common.legend=TRUE, legend = "bottom", ncol = 3, nrow = 2)
cecum_plots

ggsave(cecum_plots, device='tiff', dpi=700)

