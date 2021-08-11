# Only do these steps once. 
# If you run into "non-zero exit status" for fields (MACOSX problem), 
# you'll need to download and install this: https://cran.r-project.org/bin/macosx/tools/gfortran-6.1.pkg
# install.packages("~/Downloads/fields_10.0.tar.gz", repos = NULL, type="source")
# BiocManager::install(version="3.10")
# BiocManager::install('phyloseq')
# BiocManager::install('DESeq2')
# BiocManager::install('biomformat')
# remotes::install_github("jbisanz/qiime2R")
# install.packages("tidyverse")

# Once the above have been done once, only need to load the packages now in each new session.
library(biomformat)     # New version of biom
library(ade4)           # Required by phyloseq
library(phyloseq)       # Used to manipulate phylo data
library(tidyverse)      # Installs ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats
library(gdata)
library(data.table)
library(vegan)
library(zCompositions)  # cmultRepl 
library(selbal)
library(ape)
library(vegan)
library(scales)
library(matrixStats)

rm(list = ls())

## Read in qiime2 biom file
canola.biom <- read_biom("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/feature-table.biom")

## Extract abundance data (returns ASV by sample data frame)
abundance <- data.frame(as.matrix(biom_data(canola.biom)), stringsAsFactors = F)
can.asv.otu <- abundance
sum(colSums(abundance))
mean(colSums(abundance))

# ****************************
## Process sample information
can.asv.sam <- data.frame(names(abundance), stringsAsFactors = F)
fastq.gz <- can.asv.sam$names.abundance. # This differentiates the 2016 samples with the prefixed gibberish to deal with later
nchar(can.asv.sam$names.abundance.)[1000] - nchar(can.asv.sam$names.abundance.)[1] # 27 characters in the prefixed gibberish
can.asv.sam$names.abundance.[startsWith(can.asv.sam$names.abundance., "M")] <- substr(can.asv.sam$names.abundance.[startsWith(can.asv.sam$names.abundance., "M")],28,46)
abu.colSums <- data.frame(colSums(abundance))

# Extract all the columns to populate the sample df
Sample.Code <- can.asv.sam$names.abundance.                                # Sample code (P2IRC identifier)
SampleType <- substr(Sample.Code,1,4)                                      # Sample type (ORIG/CHCK)
Year <- substr(Sample.Code,5,5)                                            # Year (2016/2017)
Site <- rep(NA,2654)                                                       # Site (S/L/M)
RootSoil <- substr(Sample.Code,7,7)                                        # Root or soil 
Plot <- substr(Sample.Code,8,11)                                           # Plot number
Week <- substr(Sample.Code,14,15)                                          # Week
Field.Duplicate <- substr(Sample.Code,12,12)                               # Extract field dupes
CollectionType <- substr(Sample.Code,13,13)                                # Collection (C01,C02,C03 = week 3,6,9) or weekly (weeks 1-10)
Extraction.Duplicate <- substr(Sample.Code,16,16)                          # Extract extraction dupes
PCR.Duplicate <- substr(Sample.Code,17,17)                                 # Extract PCR dupes
Sequencer.Replicate <- substr(Sample.Code,18,18)                           # Extract sequencer reps
Sequencer.Plate.Or.Site <- substr(Sample.Code,19,19)                       # Extract sequencer plate (2016) or site (2017)
ReadCounts <- abu.colSums$colSums.abundance.                               # Total reads for each sample
Dupes <- rep(NA, 2654)                                                     # Dupes to populate later

# Make a new df
can.asv.sam2 <- cbind.data.frame(fastq.gz = fastq.gz,
                                 SampleID = Sample.Code,
                                 Year = Year,
                                 Site = Site,
                                 RootSoil = RootSoil,
                                 Plot = Plot,
                                 Week = Week,
                                 Field.Duplicate = Field.Duplicate,
                                 CollectionType = CollectionType,
                                 Extraction.Duplicate = Extraction.Duplicate,
                                 PCR.Duplicate = PCR.Duplicate,
                                 Sequencer.Replicate = Sequencer.Replicate,
                                 Sequencer.Plate.Or.Site = Sequencer.Plate.Or.Site,
                                 ReadCounts = ReadCounts,
                                 Dupes = Dupes,
                                 stringsAsFactors = F)

# Assign the sites and tidy that column up
can.asv.sam2$Site <- can.asv.sam2$Sequencer.Plate.Or.Site
can.asv.sam2$Site[can.asv.sam2$Site == "a" | can.asv.sam2$Site == "b" | can.asv.sam2$Site == "c"] <- "L"
can.asv.sam2 <- droplevels.data.frame(can.asv.sam2)

# Duplicates - take the sample that has the greatest sequencing depth and keep (16 = "Dupes" column)
# Note: we don't actually use the next three lines now as we sum the abundances
asv.sam <- can.asv.sam2[order(can.asv.sam2$SampleID, -abs(can.asv.sam2$ReadCounts) ), ] # Sort by id and reverse of abs(value)
asv.sam[ duplicated(asv.sam$SampleID), 15] <- 1                                         # Assign lower values as duplicates
asv.sam[ !duplicated(asv.sam$SampleID), 15] <- 0                                        # Assign higher values as keepers
asv.sam$Year[asv.sam$Year == 1] <- 2016
asv.sam$Year[asv.sam$Year == 2] <- 2017
asv.sam$Site[asv.sam$Year == 2016] <- "L"
asv.sam$Site[asv.sam$SampleID == "ORIG2CR1107C01000S"] <- "S"
asv.sam$RootSoil[asv.sam$SampleID == "ORIG2CR1107C01000S"] <- "R"
asv.sam$Week[asv.sam$SampleID == "ORIG2CR1107C01000S"] <- "03"
asv.sam$Field.Duplicate[asv.sam$SampleID == "ORIG2CR1107C01000S"] <- "N"
asv.sam$CollectionType[asv.sam$SampleID == "ORIG2CR1107C01000S"] <- "C"
asv.sam$Extraction.Duplicate[asv.sam$SampleID == "ORIG2CR1107C01000S"] <- 0
asv.sam$PCR.Duplicate[asv.sam$SampleID == "ORIG2CR1107C01000S"] <- 0
asv.sam$Sequencer.Replicate[asv.sam$SampleID == "ORIG2CR1107C01000S"] <- 0
asv.sam$Sequencer.Plate.Or.Site[asv.sam$SampleID == "ORIG2CR1107C01000S"] <- "S"

# Three 2017 Llewellyn samples were included that weren't meant to be.
# 1155: NAM-25, 1205: B.juncea, 1230: NAM-66, 1237: NAM-71. These will be removed later

# Make sure weeks are coded correctly (e.g., "C01" is actually week 3). Make sure you recode C03 first, otherwise all week 3s will be recoded to week 9
asv.sam$Week  <- with(asv.sam, ifelse(CollectionType == "C" & Week == "03", "09", Week))
asv.sam$Week  <- with(asv.sam, ifelse(CollectionType == "C" & Week == "02", "06", Week))
asv.sam$Week  <- with(asv.sam, ifelse(CollectionType == "C" & Week == "01", "03", Week))

# Add yield and line information to the sample df
can.metadata <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/canola_metadata.csv", header = T)
names(can.metadata) <- c("SampleID","SiteYearPlot","Plot","CanolaLine","Year","Site_rm","Site","Yield.kg.ha")
can.metadata$SiteYear <- as.factor(paste(can.metadata$Site, can.metadata$Year, sep = "."))
can.metadata$site.year.plot <- as.factor(paste(can.metadata$Site, can.metadata$Year, can.metadata$Plot, sep = "."))
can.metadata <- can.metadata[can.metadata$SiteYear %in% c("L.2016", "L.2017", "M.2017","S.2017"),]
asv.sam$site.year.plot <- NA
asv.sam$site.year.plot <- paste(asv.sam$Site, asv.sam$Year, asv.sam$Plot, sep = ".")

# Merge the dfs and check to see they've merged correctly
can.asv.sam3 <- merge(asv.sam, can.metadata, by = "site.year.plot", all.x = TRUE)
can.asv.sam3$Plot.y[is.na(can.asv.sam3$Plot.y)] <- as.numeric(can.asv.sam3$Plot.x[is.na(can.asv.sam3$Plot.y)]) # Four erroneous plots to be removed later 
identical(as.integer(can.asv.sam3$Plot.x), as.integer(can.asv.sam3$Plot.y)) # TRUE

# Make a new data frame
can.asv.sam4 <- data.frame(fastq.gz = can.asv.sam3$SampleID.y,
                           SampleID = can.asv.sam3$SampleID.x, 
                           Year = can.asv.sam3$Year.x, 
                           Site = can.asv.sam3$Site.x,
                           RootSoil = can.asv.sam3$RootSoil,
                           Plot = can.asv.sam3$Plot.x,
                           CanolaLine = can.asv.sam3$CanolaLine,
                           Week = can.asv.sam3$Week,
                           Field.Duplicate = can.asv.sam3$Field.Duplicate,
                           Extraction.Duplicate = can.asv.sam3$Extraction.Duplicate,
                           PCR.Duplicate = can.asv.sam3$PCR.Duplicate,
                           Sequencer.Replicate = can.asv.sam3$Sequencer.Replicate,
                           Sequencer.Plate.Or.Site = can.asv.sam3$Sequencer.Plate.Or.Site,
                           ReadCounts = can.asv.sam3$ReadCounts,
                           Dupes = can.asv.sam3$Dupes,
                           Yield.kg.ha = can.asv.sam3$Yield.kg.ha,
                           stringsAsFactors = F)
can.asv.sam4$SiteYearLinePlotWeekRootSoil <- paste(can.asv.sam4$Site, can.asv.sam4$Year, can.asv.sam4$CanolaLine, can.asv.sam4$Plot, can.asv.sam4$Week, can.asv.sam4$RootSoil, sep = ".")
nlevels(as.factor(can.asv.sam4$SiteYearLinePlotWeekRootSoil)); print("of"); length(as.factor(can.asv.sam4$SiteYearLinePlotWeekRootSoil)) # 2122 / 2654 are unique (532 duplicates)

# ****************************
## Process taxonomic information
taxa <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/taxonomy.csv", header = T)
taxa2 <- taxa %>% separate(taxonomy, into = c("FeatureID","kingdom","phylum","class","order","family","genus","species","subspecies"), sep = ";")

# Remove the superfluous characters in each name ("D_0__", "D_1__", etc.)
Kingdom <- gsub("D_0__", "", taxa2$kingdom)
Phylum <- gsub("D_1__", "", taxa2$phylum)
Class <- gsub("D_2__", "", taxa2$class)
Order <- gsub("D_3__", "", taxa2$order)
Family <- gsub("D_4__", "", taxa2$family)
Genus <- gsub("D_5__", "", taxa2$genus)
Species <- gsub("D_6__", "", taxa2$species)

# Combine each taxonomic level into a matrix
taxa3 <- as.matrix(cbind(kingdom = Kingdom, phylum = Phylum, class = Class, order = Order, family = Family, genus = Genus, species = Species))

## Process the taxonomy information
# The subset_taxa command will remove rows with the specified taxonomy, but also NAs. So recode the NAs to "unclassified" here to avoid removing extra rows.
# https://github.com/joey711/phyloseq/issues/683
can.asv.tax <- as.data.frame(taxa3)
can.asv.tax <- apply(can.asv.tax, 2, 
                     function(x) 
                       gsub("^$|^ $", "unclassified", x))
can.asv.tax[can.asv.tax == "__"] <- "unclassified"
can.asv.tax <- as.matrix(can.asv.tax)                                   # Convert to matrix
row.names(can.asv.tax) <- row.names(abundance)

# ****************************
## Process sample information
# Check the row/col sums to make sure things are still peachy - here the df is 2663 x 18686
# Calculate the column sums
canola.asv.colsums <- colSums(can.asv.otu)
canola.asv.rowsums <- rowSums(can.asv.otu)
# Count the non-zero sum columns
length(canola.asv.colsums[canola.asv.colsums != 0]); print("of"); length(canola.asv.colsums) # 2654 of 2654 samples
length(canola.asv.rowsums[canola.asv.rowsums != 0]); print("of"); length(canola.asv.rowsums) # 17345 of 17345 taxa

## Check that the row and column names match up
identical(rownames(can.asv.otu),rownames(can.asv.tax)) # Yass
rownames(can.asv.sam4) <- can.asv.sam4[,2]
identical(names(can.asv.otu),rownames(can.asv.sam4)) # Nay
or.list <- names(can.asv.otu)
can.asv.sam5 <- can.asv.sam4[match(or.list, can.asv.sam4$SampleID),]
identical(names(can.asv.otu),rownames(can.asv.sam5)) # Yass

# *************
## Make phyloseq class. Double check to see if the taxonomy columns are capitalized or not. This will affect the coding below.
can_all.asv <- phyloseq(otu_table(can.asv.otu, taxa_are_rows = TRUE), sample_data(can.asv.sam5), tax_table(can.asv.tax)) # 17345 taxa x 2654 samples

# Inventory the number of samples and duplicates for roots and soil in each year
length(can_all.asv@sam_data$RootSoil[can_all.asv@sam_data$RootSoil == "R"]) # 1330
length(can_all.asv@sam_data$RootSoil[can_all.asv@sam_data$RootSoil == "S"]) # 1324
length(can_all.asv@sam_data$RootSoil[can_all.asv@sam_data$RootSoil == "R" & can_all.asv@sam_data$Year == 2016]) # 576
length(can_all.asv@sam_data$RootSoil[can_all.asv@sam_data$RootSoil == "R" & can_all.asv@sam_data$Year == 2017]) # 754
length(can_all.asv@sam_data$RootSoil[can_all.asv@sam_data$RootSoil == "R" & can_all.asv@sam_data$Year == 2016 & can_all.asv@sam_data$Field.Duplicate == "Y"])   # 37 field duplicates in 2016
length(can_all.asv@sam_data$RootSoil[can_all.asv@sam_data$RootSoil == "R" & can_all.asv@sam_data$Year == 2017 & can_all.asv@sam_data$Field.Duplicate == "Y"])   # 50 field duplicates in 2017
length(can_all.asv@sam_data$RootSoil[can_all.asv@sam_data$RootSoil == "R" & can_all.asv@sam_data$Year == 2016 & can_all.asv@sam_data$Extraction.Duplicate > 0]) # 2 extraction duplicates in 2016
length(can_all.asv@sam_data$RootSoil[can_all.asv@sam_data$RootSoil == "R" & can_all.asv@sam_data$Year == 2017 & can_all.asv@sam_data$Extraction.Duplicate > 0]) # 12 extraction duplicates in 2016
length(can_all.asv@sam_data$RootSoil[can_all.asv@sam_data$RootSoil == "R" & can_all.asv@sam_data$Year == 2016 & can_all.asv@sam_data$PCR.Duplicate > 0])        # 7 extraction duplicates in 2016
length(can_all.asv@sam_data$RootSoil[can_all.asv@sam_data$RootSoil == "R" & can_all.asv@sam_data$Year == 2017 & can_all.asv@sam_data$PCR.Duplicate > 0])        # 28 extraction duplicates in 2016
length(can_all.asv@sam_data$RootSoil[can_all.asv@sam_data$RootSoil == "R" & can_all.asv@sam_data$Year == 2016 & can_all.asv@sam_data$Sequencer.Replicate > 0])  # 46 extraction duplicates in 2016
length(can_all.asv@sam_data$RootSoil[can_all.asv@sam_data$RootSoil == "R" & can_all.asv@sam_data$Year == 2017 & can_all.asv@sam_data$Sequencer.Replicate > 0])  # 80 extraction duplicates in 2016

# Subset to just root bacteria
can_root.asv <- subset_samples(can_all.asv, RootSoil == "R") # 17345 x 1330

# Remove archaea, chloroplasts, and mitochondria (plant DNA)
can_root.asv <- subset_taxa(can_root.asv, kingdom != "Archaea")                             # Remove archaea (17345 x 1330)
can_root.asv <- subset_taxa(can_root.asv, class != "Chloroplast" & order != "Chloroplast")  # Remove chloroplasts (16002 x 1330)
can_root.no_mito <- subset_taxa(can_root.asv, family != "Mitochondria")                     # Remove mitochondrial contaminant (15342 x 1330)

# Remove erroneuously-included plots and create Aliivibrio and non datasets
can_root.no_mito@sam_data$SiteYearPlot <- paste(can_root.no_mito@sam_data$Site, 
                                                can_root.no_mito@sam_data$Year, 
                                                can_root.no_mito@sam_data$Plot, sep = ".")  # Create a grouping factor to use to remove errorneous samples
can_root.no_mito <- subset_samples(can_root.no_mito, SiteYearPlot != "L.2017.1155")         # Remove the erroneously included Llewellyn plot (NAM-25)
can_root.no_mito <- subset_samples(can_root.no_mito, SiteYearPlot != "L.2017.1205")         # Remove the erroneously included Llewellyn plot (B.juncea)
can_root.no_mito <- subset_samples(can_root.no_mito, SiteYearPlot != "L.2017.1230")         # Remove the erroneously included Llewellyn plot (NAM-66)
can_root.no_mito <- subset_samples(can_root.no_mito, SiteYearPlot != "L.2017.1237")         # Remove the erroneously included Llewellyn plot (NAM-71)
no_mito.fischeri <- subset_taxa(can_root.no_mito, genus == "Aliivibrio" )                   # Collect Aliivibrio internal standard (25 x 1326)
no_mito.no_fischeri <- subset_taxa(can_root.no_mito, genus != "Aliivibrio" )                # Remove Aliivibrio internal standard (13230 x 1326)

# 5 x 4
step0.hist <- as.data.frame(otu_table(no_mito.no_fischeri@otu_table))
step0.hist.unlist <- unlist(step0.hist)
par(mar = c(4,4,1,1))
step0.hist.out <- hist(log(step0.hist.unlist+1), breaks = seq(0,20,0.5), freq = FALSE, xlim = c(0,20), ylim = c(0,0.025), main = "", ylab = "", xlab = "")

x <- rownames(step0.hist)
x[grepl("^[[:digit:]]", x)] <- paste("X", x[grepl("^[[:digit:]]", x)], sep = "")
rownames(step0.hist) <- x

write.csv(step0.hist, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/step0.taxa.csv")

## Filter nuisance ASVs:

######################################################################
######################################################################
## 1) This step removes imaginary ASVs, which we say were only present twice in the entire data set 
# Noise filter: removes 13230 - 5754 = 7476 ASVs
# Remove any taxa with total # reads < 2 (5754 x 1326) - col sums > 2
no_mito.no_fischeri <- prune_taxa(taxa_sums(no_mito.no_fischeri) > 2, no_mito.no_fischeri)            # 13230 x 1326 to 5754 x 1326
# Remove zero-sum samples 
no_mito.no_fischeri <- prune_samples(sample_sums(no_mito.no_fischeri) > 0, no_mito.no_fischeri)       # 5754 x 1326 to 5754 x 1319

# Output data for MAMETetal_SBB_figures.R
step1.hist <- as.data.frame(otu_table(no_mito.no_fischeri@otu_table))
step1.hist.unlist <- unlist(step1.hist)
par(mar = c(4,4,1,1))
hist(log(step0.hist.unlist+1), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000), main = "", ylab = "", xlab = "", col = "grey90")
par(new = T)
step1.hist.out <- hist(log(step1.hist.unlist+1), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000), main = "", ylab = "", xlab = "", col = "#b94663")
x <- rownames(step1.hist)
x[grepl("^[[:digit:]]", x)] <- paste("X", x[grepl("^[[:digit:]]", x)], sep = "")
rownames(step1.hist) <- x
write.csv(step1.hist, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/step1.taxa.csv")

######################################################################
######################################################################
## 2) This step removes ASVs that only occurred in 2 plots, which are considered too rare to assess
# Trace filter: removes 5754 - 3681 = 2073 ASVs; 2/1319 = 0.15% of samples
no_mito.no_fischeri <- filter_taxa(no_mito.no_fischeri, function (x) {sum(x > 0) > 2}, prune=TRUE)    # Remove prevalence singletons and doubletons (3681 x 1319)
no_mito.fischeri <- subset_samples(no_mito.fischeri, 
                                   colnames(otu_table(no_mito.fischeri)) %in% 
                                     colnames(otu_table(no_mito.no_fischeri)))                        # Ensure the Aliivibrio dataset contains the same samples as the non (25 x 1319)

no_mito.no_fischeri@sam_data$SiteYearLineWeek <- paste(no_mito.no_fischeri@sam_data$Site, 
                                                       no_mito.no_fischeri@sam_data$Year, 
                                                       no_mito.no_fischeri@sam_data$CanolaLine, 
                                                       no_mito.no_fischeri@sam_data$Week, sep = ".")  # Create a grouping factor to tally total and non-zero abundances in each site GxT
asv.root.0 <- as.data.frame(no_mito.no_fischeri@otu_table)                                            # Extract the abundance table
sam.root.0 <- as.data.frame(no_mito.no_fischeri@sam_data)                                             # Extract the sample table
asv.root.1 <- cbind.data.frame(SiteYearLineWeek = sam.root.0$SiteYearLineWeek, t(asv.root.0))         # Combine the grouping factor and abundance table into one df

# Output data for MAMETetal_SBB_figures.R
step2.hist <- asv.root.0
step2.hist.unlist <- unlist(step2.hist)
par(mar = c(4,4,1,1))
hist(log(step1.hist.unlist+1), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000), main = "", ylab = "", xlab = "", col = "#b94663")
par(new = T)
step2.hist.out <- hist(log(step2.hist.unlist+1), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000),  main = "", ylab = "", xlab = "", col = "#bc7d39")
x <- rownames(step2.hist)
x[grepl("^[[:digit:]]", x)] <- paste("X", x[grepl("^[[:digit:]]", x)], sep = "")
rownames(step2.hist) <- x

write.csv(step2.hist, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/step2.taxa.csv")

# Count the number of non-zero abundances in each site GxT. This may take a while.
root.non.zeros <- asv.root.1 %>%
  group_by(SiteYearLineWeek) %>% 
  summarise_all(list(~ sum(as.logical(.))))

# Count the total number of observations for each site GxT
root.totals <- asv.root.1 %>%
  group_by(SiteYearLineWeek) %>% 
  tally()

# Divide the non-zero count by the total count to get a proportion. 
# E.g., at L.2016, there are 3 reps of NAM-0. For an ASV to be included, it should be found in 2/3 samples (prop = 0.6).
# A problem arises when there are an even number of reps/dupes.
# E.g., if there are 4 "samples" (e.g., 3 reps and 1 dupe), we would still like to keep the ASV if it occurs in 2/4 samples.
# So we need to keep >= 0.5, rather than >0.5 or we'll lose the 2/4, 3/6, etc. 
# Also, we should reclass occurrences of 1 to 0 before the proportion step so we don't include the 1/2s.
root.non.zeros2 <- root.non.zeros
root.non.zeros2[root.non.zeros2 == 1] <- 0

root.props <- root.non.zeros2[,-1] / root.totals$n
rownames(root.props) <- root.totals$SiteYearLineWeek
# In this step, if set to >= 0.5, for a few samples there are only 2 reps (1136 ASVs). If > 0.5, down to 667 ASVs
# If present at >=0.3 = 3515; at >= 1 = 176
# 3) This is step that insures that the ASVs are present in the majority of all plots that make up a sample
# Majority filter: only removes 3681-980 = 2701 ASVs.
root.pa <- as.data.frame(ifelse(root.props >= 0.5, 1, 0)) # Convert proportions >= 0.5 to 1, everything elseo to 0
root.pa2 <- root.pa[,colSums(root.pa) > 0]                # Remove all taxa that have been codified as 0
sub.list <- names(root.pa2)                               # Extract the ASV id list to use to subset the phyloseq object later

######################################################################
## Prepare an Aliivibrio standardized dataset

# Prepare the Aliivibrio and non datasets
can_root.fischeri <- no_mito.fischeri
can_root.no_mito.no_fischeri <- subset_taxa(no_mito.no_fischeri, rownames(tax_table(no_mito.no_fischeri)) %in% sub.list) # 980 x 1319

# Using internal standard to correct abundances
# Formula is: X = Rs*(wi/gi*ci)/Ri
# Where:
# Ri (fischeri) is the number of reads assigned to the internal standard
# Rs (asv.table) is the number of reads assigned to other taxa
# wi is the weight of internal standard gDNA add
# gi is the weight of the genome of internal standard
# ci is the 16S copy number of the internal standard
# X is the (initirooty unknown) number of 16S rRNA genes per sample
fischeri <- as.vector(sample_sums(can_root.fischeri))                               # Create vector of internal standard
wi <- 1e-10                                                                         # Weight of internal standard gDNA added to the samples
gi <- 4.49e-15                                                                      # Weight of genome of the internal standard
ci <- 8                                                                             # 16S copy number of the internal standard
adj <- (wi/gi)*ci                                                                   # Create adjustment metric
fischeri.1 <- adj/fischeri                                                          # Basically: (wi/gi)*ci / Ri
fischeri.1[is.infinite(fischeri.1)] <- 0                                            # Replace Inf values with 0

asv.table <- as.matrix(can_root.no_mito.no_fischeri@otu_table)                      # Export asv abundance table as matrix
asv.table.norm <- as.data.frame(sweep(asv.table, 2, fischeri.1, `*`))               # Basically Rs * ((wi/gi)*ci / Ri)
tax.table.norm <- as.matrix(can_root.no_mito.no_fischeri@tax_table)                 # This extracts the taxonomic information excluding mitochondria & aliivibrio
can_root.norm.fischeri.std <- phyloseq(otu_table(asv.table.norm, taxa_are_rows = TRUE), 
                                       sample_data(can.asv.sam5), tax_table(tax.table.norm))      # Recreating a phyloseq object (980 taxa, 1319 samples)

# Make dfs for abundances and sample information
asv.root.3 <- as.data.frame(can_root.norm.fischeri.std@otu_table)
sam.root.3 <- as.data.frame(can_root.norm.fischeri.std@sam_data)

# Output data for MAMETetal_SBB_figures.R
step3.hist <- as.data.frame(asv.table)
step3.hist.unlist <- unlist(step3.hist)
par(mar = c(4,4,1,1))
hist(log(step2.hist.unlist+1), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000), main = "", ylab = "", xlab = "", col = "#bc7d39")
par(new = T)
step3.hist.out <- hist(log(step3.hist.unlist+1), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000),  main = "", ylab = "", xlab = "", col = "#6fac5d")
x <- rownames(step3.hist)
x[grepl("^[[:digit:]]", x)] <- paste("X", x[grepl("^[[:digit:]]", x)], sep = "")
rownames(step3.hist) <- x

write.csv(step3.hist, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/step3.taxa.csv")

######################################################################
######################################################################
###### Step 4) MERGE SAMPLES
## Here the samples will combined to make one per site-genotype combination.
## This is to match the coarse resolution given by the yield stability calculation 
## to the finer resolved bacterial community data. 
# Step 1: Merge duplicates. 
# Step 2: Merge replicates.

## Step 1
can_root.merge <- merge_samples(can_root.norm.fischeri.std, "SiteYearLinePlotWeekRootSoil", 
                                fun = sum)                                           # This step will sum everything (factors become NAs, weeks, lines, etc. get summed. Need to fix later)
can_root.no.dup <- can_root.merge
can.asv.root.0 <- as.matrix(can_root.no.dup@otu_table)
can.tax.root.0 <- as.matrix(can_root.no.dup@tax_table)
can.sam.root.0 <- data.frame(can_root.no.dup@sam_data)

# The sample df has been lost as we can't meaningfully sum characters/factors
# Create a new sample df
Sample.ID <- c(rownames(can.sam.root.0))
split.df <- unlist(strsplit(Sample.ID,"[.]"))
Site <- substr(Sample.ID,1,1)
Year <- substr(Sample.ID,3,6)
CanolaLine <- split.df[seq(3,length(split.df),6)]
Plot <- split.df[seq(4,length(split.df),6)]
Week <- split.df[seq(5,length(split.df),6)]
RootSoil <- split.df[seq(6,length(split.df),6)]
SiteYearPlot <- paste(Site, Year, Plot, sep = ".")
SiteYearLineWeekRootSoil <- paste(Site, Year, CanolaLine, Week, RootSoil, sep = ".")
can.sam.root.1 <- cbind.data.frame(Site = Site, Year = Year, Plot = Plot, CanolaLine = CanolaLine, Week = Week, RootSoil = RootSoil, SiteYearPlot = SiteYearPlot, SiteYearLineWeekRootSoil = SiteYearLineWeekRootSoil)
rownames(can.sam.root.1) <- Sample.ID
can_root.no.dup1 <- phyloseq(otu_table(can.asv.root.0, taxa_are_rows = TRUE), 
                             sample_data(can.sam.root.1), tax_table(can.tax.root.0))      # Recreating a phyloseq object (980 x 1057)
can_root.no.dup2 <- prune_taxa(taxa_sums(can_root.no.dup1) != 0, can_root.no.dup1)    # Remove any taxa that have root zero abundances (980 x 1057)

## Step 2
can_root.merge <- merge_samples(can_root.no.dup2, "SiteYearLineWeekRootSoil", 
                                fun = sum)                                           # This step will sum replicates for each line (980 x 360)
can.asv.root.merge <- as.matrix(can_root.merge@otu_table)
can.tax.root.merge <- as.matrix(can_root.merge@tax_table)
can.sam.root.merge <- data.frame(can_root.merge@sam_data)

# Create a new sample df as before
Sample.ID <- c(rownames(can.sam.root.merge))
split.df <- unlist(strsplit(Sample.ID,"[.]"))
Site <- substr(Sample.ID,1,1)
Year <- substr(Sample.ID,3,6)
CanolaLine <- split.df[seq(3,length(split.df),5)]
Week <- split.df[seq(4,length(split.df),5)]
RootSoil <- split.df[seq(5,length(split.df),5)]
SiteYearLine <- paste(Site, Year, CanolaLine, sep = ".")
SiteYearLineWeekRootSoil <- paste(Site, Year, CanolaLine, Week, RootSoil, sep = ".")
can.sam.root.2 <- cbind.data.frame(Site = Site, Year = Year, CanolaLine = CanolaLine, Week = Week, RootSoil = RootSoil, SiteYearLine = SiteYearLine, SiteYearLineWeekRootSoil = SiteYearLineWeekRootSoil)
rownames(can.sam.root.2) <- Sample.ID

# Create a new phyloseq object and make sure there are no zero-sum samples
can_root.no.dup4 <- phyloseq(otu_table(can.asv.root.merge, taxa_are_rows = TRUE), 
                             sample_data(can.sam.root.2), tax_table(can.tax.root.merge))      # Recreating a phyloseq object (980 x 360)
can_root.no.dup5 <- can_root.no.dup4
can_root <- prune_samples(sample_sums(can_root.no.dup5) != 0, can_root.no.dup5)             # Remove zero-sum samples (980 x 360)

# Output data for MAMETetal_SBB_figures.R
step4.hist <- as.data.frame(t(otu_table(can_root@otu_table, taxa_are_rows = FALSE)))
step4.hist.unlist <- unlist(step4.hist)
par(mar = c(4,4,1,1))
hist(log(step3.hist.unlist+1), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000), main = "", ylab = "", xlab = "", col = "#6fac5d")
par(new = T)
step4.hist.out <- hist(log(step4.hist.unlist+1), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000),  main = "", ylab = "", xlab = "", col = "#677ad1")
x <- rownames(step4.hist)
x[grepl("^[[:digit:]]", x)] <- paste("X", x[grepl("^[[:digit:]]", x)], sep = "")
rownames(step4.hist) <- x

write.csv(step4.hist, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/step4.taxa.csv")


######################################################################
######################################################################
# Step 5) In this step, we create a list that is the union of meaningful ASVs present in at least one of our weeks of interest

can.root.0 <- can_root
can.root.w3 <- subset_samples(can.root.0, Week == "03")                     # Subset week 3 (980 x 64)
can.root.w6 <- subset_samples(can.root.0, Week == "06")                     # Subset week 6 (980 x 64)
can.root.w9 <- subset_samples(can.root.0, Week == "09")                     # Subset week 9 (980 x 64)
can.root.w3.1 <- prune_samples(sample_sums(can.root.w3) != 0, can.root.w3)  # Remove zero-sum samples (980 x 64)
can.root.w6.1 <- prune_samples(sample_sums(can.root.w6) != 0, can.root.w6)  # Remove zero-sum samples (980 x 64)
can.root.w9.1 <- prune_samples(sample_sums(can.root.w9) != 0, can.root.w9)  # Remove zero-sum samples (980 x 64)
can.root.w3.1 <- prune_taxa(taxa_sums(can.root.w3) > 2, can.root.w3)        # Remove any taxa with total # reads < 2 (795 x 64)
can.root.w6.1 <- prune_taxa(taxa_sums(can.root.w6) > 2, can.root.w6)        # Remove any taxa with total # reads < 2 (750 x 64)
can.root.w9.1 <- prune_taxa(taxa_sums(can.root.w9) > 2, can.root.w9)        # Remove any taxa with total # reads < 2 (832 x 64)
asv.root.w3 <- as.data.frame(can.root.w3.1@otu_table)                       # Extract the abundance table
asv.root.w6 <- as.data.frame(can.root.w6.1@otu_table)                       # Extract the abundance table
asv.root.w9 <- as.data.frame(can.root.w9.1@otu_table)                       # Extract the abundance table
range(asv.root.w3[asv.root.w3 != 0]) # 1.054407e+01 2.002565e+07            # Range of non-zero absolute abundance in a given cell
range(asv.root.w6[asv.root.w6 != 0]) # 9.550478e+00 4.472160e+07            # Range of non-zero absolute abundance in a given cell
range(asv.root.w9[asv.root.w9 != 0]) # 8.561516e+00 1.421349e+07            # Range of non-zero absolute abundance in a given cell
range(colSums(asv.root.w3 != 0)) # 1 to 60                                  # Range of taxa occurrence in a given sample (N = 64)
range(colSums(asv.root.w6 != 0)) # 1 to 62                                  # Range of taxa occurrence in a given sample (N = 64)
range(colSums(asv.root.w9 != 0)) # 1 to 64                                  # Range of taxa occurrence in a given sample (N = 64)
# hist(unlist(colSums(asv.root.w3 != 0)))
asv.root.w3.names <- names(asv.root.w3)                                     # Extract the list of taxa to union
asv.root.w6.names <- names(asv.root.w6)                                     # Extract the list of taxa to union
asv.root.w9.names <- names(asv.root.w9)                                     # Extract the list of taxa to union
global_union <- Reduce(union, list(asv.root.w3.names, asv.root.w6.names, asv.root.w9.names)) # 977 ASVs across the weeks 3, 6, 9
asv.global <- prune_taxa(global_union, can_root) # 977 x 360                # Subset the phyloseq object to the 977 taxa determined above
range(taxa_sums(asv.global))    # 5.173945e+02 1.482910e+08                 # Range of total abundances for a given taxa
range(sample_sums(asv.global))  # 65.799772e+03 1.211581e+08                # Range of total abundances for a given sample
asv.root <- as.data.frame(asv.global@otu_table)                             # Extract the abundance table
tax.root <- as.matrix(asv.global@tax_table)                                 # Extract the taxonomy table
sam.root <- data.frame(asv.global@sam_data)                                 # Extract the sample table
range(colSums(asv.root != 0)) # 1 to 324: lowest prevalence ASV occurs in 1 samples
range(rowSums(asv.root != 0)) # 4 to 364: lowest richness for a sample is 4 taxa

# Sample quality
# Here we look at the richness of a sample and remove richness <6 across 3 plots - considered a failed sequencing run.
hist(unlist(rowSums(asv.root != 0)), breaks = seq(0,400,1), xlim = c(0,10))
# Remove suspect samples, n=3 samples,
# then remove ASVs that only occur in these samples
asv.root.0 <- asv.root[rowSums(asv.root != 0)>6,] # 357 x 977
range(colSums(asv.root.0 != 0)) # 1 to 324: lowest prevalence ASV occurs in 1 sample
range(rowSums(asv.root.0 != 0)) # 7 to 364: lowest richness for a sample is 7 taxa

# Output data for MAMETetal_SBB_figures.R
step5.hist <- asv.root.0
step5.hist.unlist <- unlist(step5.hist)
par(mar = c(4,4,1,1))
hist(log(step4.hist.unlist+1), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000), main = "", ylab = "", xlab = "", col = "#677ad1")
par(new = T)
step5.hist.out <- hist(log(step5.hist.unlist+1), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000),  main = "", ylab = "", xlab = "", col = "#9750a1")
x <- rownames(step5.hist)
x[grepl("^[[:digit:]]", x)] <- paste("X", x[grepl("^[[:digit:]]", x)], sep = "")
rownames(step5.hist) <- x

write.csv(step5.hist, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/step5.taxa.csv")

# Output the data necessary for the winnowing pipeline, which is the next step
write.csv(asv.root.0, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/asv.root.reps_dupes_merged.prev_and_abun_filt.csv")
write.csv(tax.root,   "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/tax.root.reps_dupes_merged.prev_and_abun_filt.csv")
write.csv(sam.root,   "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/sam.root.reps_dupes_merged.prev_and_abun_filt.csv")

# Read in yield information
can.yield <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/canola.sam.root.2.csv", header = T, row.names = 1)
rownames(can.yield) <- can.yield$SampleID

# ASVs common across devem,opmental stages
global_union1 <- Reduce(union, list(names(asv.root.w3), names(asv.root.w6), names(asv.root.w9))) # 977 ASVs across the weeks 3, 6, 9
sam.root.w3 <- sam.root[rownames(asv.root.w3),]
sam.root.w6 <- sam.root[rownames(asv.root.w6),]
sam.root.w9 <- sam.root[rownames(asv.root.w9),]
sam.root.w3 <- merge(sam.root.w3, can.yield[,c((ncol(can.yield)-2):(ncol(can.yield)))], by = 0)
sam.root.w6 <- merge(sam.root.w6, can.yield[,c((ncol(can.yield)-2):(ncol(can.yield)))], by = 0)
sam.root.w9 <- merge(sam.root.w9, can.yield[,c((ncol(can.yield)-2):(ncol(can.yield)))], by = 0)
range(colSums(asv.root.w3 != 0))
range(colSums(asv.root.w6 != 0))
range(colSums(asv.root.w9 != 0))

write.csv(asv.root.w3, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/asv.root.reps_dupes_merged.prev_and_abun_filt.w3.csv")
write.csv(asv.root.w6, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/asv.root.reps_dupes_merged.prev_and_abun_filt.w6.csv")
write.csv(asv.root.w9, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/asv.root.reps_dupes_merged.prev_and_abun_filt.w9.csv")
write.csv(sam.root.w3, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/sam.root.reps_dupes_merged.prev_and_abun_filt.w3.csv")
write.csv(sam.root.w6, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/sam.root.reps_dupes_merged.prev_and_abun_filt.w6.csv")
write.csv(sam.root.w9, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/sam.root.reps_dupes_merged.prev_and_abun_filt.w9.csv")

# Replace zeroes
asv.root.w3.sq <- cmultRepl(asv.root.w3, method = "SQ", output = "p-counts") # 30 corrected values
asv.root.w6.sq <- cmultRepl(asv.root.w6, method = "SQ", output = "p-counts") # 416 corrected values
asv.root.w9.sq <- cmultRepl(asv.root.w9, method = "SQ", output = "p-counts") # 331 corrected values
write.csv(asv.root.w3.sq, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/asv.root.reps_dupes_merged.prev_and_abun_filt.sq.w3.csv")
write.csv(asv.root.w6.sq, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/asv.root.reps_dupes_merged.prev_and_abun_filt.sq.w6.csv")
write.csv(asv.root.w9.sq, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/asv.root.reps_dupes_merged.prev_and_abun_filt.sq.w9.csv")

