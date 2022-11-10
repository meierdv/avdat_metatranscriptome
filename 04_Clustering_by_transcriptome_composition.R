require(data.table)
require(vegan)

Rewetting_metadata <- read.table("Metadata.csv", 
                                 sep="\t", 
                                 header=T, 
                                 row.names=1, 
                                 quote="")[c(1:24),]

#Clustering timepoints by expression pattern across all MAGs
Feature_counts <- read.csv("Avdat_featureCounts.tsv", 
                           sep="\t", 
                           header=T,
                           row.names = "Geneid",
                           check.names = F, 
                           quote="")[,c(3:26)]

Feature_counts_dist <- vegdist(t(decostand(Feature_counts, 
                                           method="total", 
                                           MARGIN=2)),
                               method="jaccard")

MAG_Composition <- read.csv("Avdat_TPMs_per_bin.tsv",
                            sep="\t",
                            header=T,
                            check.names = F,
                            row.names="Bin"
                            )

Composition_dist <- vegdist(t(MAG_Composition), method="jaccard")

TPMs <- read.csv("Avdat_per_bin_normalized_TPMs_DeSeq2.tsv", 
                 sep="\t", 
                 header=T,
                 row.names = "Geneid",
                 check.names = F, 
                 quote="")[,c(14:37)]

TPM_dist <- vegdist(t(TPMs), method="jaccard")

TPM_cluster <- hclust(TPM_dist, method="average")
plot(TPM_cluster, hang="-1", ylab="Dissimilarity")

#---------------------------------------------------------------------------------------------------
# Ordination plot, NMDS
#--------------------------------------------------------------------------------------------------
require(ggplot2)

comMDS <- metaMDS(TPM_dist, engine = c("monoMDS"))  # calculates NMDS ordination for the distance matrix

#Checking the goodness of fit and the stress of the calculated NMDS
goodness_comMDS <- goodness(comMDS)
stress_comMDS <- stressplot(comMDS)

#Creating an object plottable by GGplot. It is basically points with coordinates, which are extracted from the NMDS object.
#It's a datframe with x-axis called "MDS1" and y-axis called "MDS2"
com_ggNMDS <- data.frame(MDS1=comMDS$points[,1], 
                         MDS2=comMDS$points[,2])

#Habit_colours <- c("#003333", "#CCFF00", "#669900", "#00CCCC", "#FFFF00")  # Defining a colour vector corresponding to the sample categories in alphabetic order
#names(Habit_colours) <- cat(com_logic$Habitat_category) # assigning the sample categories to the colour vector
#Habit_colScale <- scale_colour_manual(name = "Habitat category", values = Habit_colours) # creating a manual colour scale for the plot

#Generating basic plot object
ggNMDS_plot <- ggplot(data=com_ggNMDS, 
                      aes(x=MDS1, 
                          y=MDS2, 
                          colour=Rewetting_metadata$Phase2), 
                      label=rownames(com_ggNMDS))

#Adding various attributes and settings to the plot
ggNMDS_plot +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  scale_y_continuous(limits = c(-0.4, 0.4)) +
  geom_point(size=6.5, 
             col="black") +
  geom_point(size=6) + 
  #  Habit_colScale +
  geom_text(size=4.5, 
            col="black", 
            aes(label=rownames(com_ggNMDS)), 
            hjust=0.5, 
            vjust=-1.5) +
  theme(legend.text=element_text(size=14, 
                                 lineheight=15), 
        legend.title=element_text(size=15, 
                                  lineheight=20))


  
#Analysis of similarities to test if assignment of sample categories holds true
anosim(t(TPMs), Rewetting_metadata$Phase2, 
       permutations=9999,
       distance="jaccard")

#perMANOVA
adonis2(t(MAG_Composition) ~ Phase2 + Timepoint + Water_content + Replicate,
        data=Rewetting_metadata,
        permutations=9999,
        method="jaccard")


