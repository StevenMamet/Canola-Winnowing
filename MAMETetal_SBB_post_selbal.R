## If necessary
# install.packages("remotes")
# remotes::install_github("tidyverse/dplyr")
# devtools::install_github("kylebittinger/usedist")
# devtools::install_github(repo = "malucalle/selbal")
# devtools::install_github("fkeck/phylosignal")
# BiocManager::install("ggtree")

library(vegan)          # Relative abundances
library(tidyverse)      # Tidy data manipulation
library(scales)         # Transparency in plotting
library(selbal)         # Used here for bal.value()
library(zCompositions)  # Zero-replacement
library(car)
library(Metrics)        # RMSE
library(picante)        # Mean nearest taxon distance (mntd, comdistnt)
library(boot)           # Bootstrapped confidence intervals
library(eulerr)         # This is the best package I've found for generating Venn diagrams
library(usedist)        # For working with distance objects (e.g., bMNTD)
library(vioplot)        # Violin plots
library(wesanderson)    # Beautiful color palettes

# Phylogenetic tree packages
library(ggimage)
library(BiocManager)
library(ggtree)
library(picante)
library(phytools)
library(ape)
library(adephylo)
library(phylobase)
library(phylosignal) # Needed for phylo4d plotting

# Path analysis
library(lavaan)
library(semPlot)
library(OpenMx)
library(knitr)
library(kableExtra)
library(GGally)
library(parameters)

rm(list = ls())         # Cleanliness is godliness

##-------------------------------------------------------------------------------------------

#### Workflow
# 1. Read in the data. Calculate relative abundances. Make sure ASV IDs match among the the dfs.
# 2. Read in the selbal output from the filtered full and weekly dataset (sample IDs, yield performance, balance).
# 3. Run univariate linear models to quantify the relationship between the balances and performance.
# 4. Multiple linear regression models to quantify the relationship between weekly communities combined influence on performance.
# 5. Plot linear relationships and visualize variance for each weekly community.
# 6. Read in the selbal bal.tab output for the full filtered and weekly datasets. Extract taxa lists and count the number present in the full-filtered and each weekly community.
# 7. Build dataframe of three weekly communities (balances) over time. Includes SQ zero-replacement on full abundance table, log-transformation, and bal.value() to build balances.
# 8. Linear models to evaluate how well the global balance of the 3 communities predicts performance in other weeks.
# 9. Taxonomic overlap among three communities.
# 10. Absolute abundances, RMSE, and prev-abundances for figure generation.
# 11. Phylogenetic analyses (MNTD and bMNTD), figure.
# 12. SEMs
# 13. Comparison of 63 ASVs to exudates

## 1.-------------------------------------------------------------------------------------------
# Read in the filtered data set. Five filters, 357 samples, 977 taxa.
can.asv <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/asv.root.reps_dupes_merged.prev_and_abun_filt.csv", row.names = 1)
can.tax <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/tax.root.reps_dupes_merged.prev_and_abun_filt.csv", row.names = 1)
can.sam <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/sam.root.reps_dupes_merged.prev_and_abun_filt.csv", row.names = 1)

# Look at means, sums, etc.
mean(rowMeans(can.asv))
sum(can.asv)
mean(rowSums(can.asv != 0))
range(rowSums(can.asv != 0))
sum(can.asv != 0)
dim(can.asv)[1]*dim(can.asv)[2]

# Calculate relative abundance
can.asv.rel0 <- decostand(can.asv, method = "total")*100
rowSums(can.asv.rel0) # Double check that everything sums to 100%

# Do the taxonomy ASV identifiers match the abundance identifiers?
identical(names(can.asv), rownames(can.tax)) # Nope, b/c of the prepended "X" issue.

# Add "X" to any ASV identifier that begins with a number.
x <- rownames(can.tax)
x[grepl("^[[:digit:]]", x)] <- paste("X", x[grepl("^[[:digit:]]", x)], sep = "")
rownames(can.tax) <- x
identical(names(can.asv), rownames(can.tax)) # Yep, now they match.

## 2.-------------------------------------------------------------------------------------------
# Read in the selbal output for the weekly datasets (performance metric 1, weeks 3, 6, and 9)
LMS_GB_full_imp1_w3 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/LMS_2016_17_GlobalBalance_full_imp1.w3.csv", row.names = 1)
LMS_GB_full_imp1_w6 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/LMS_2016_17_GlobalBalance_full_imp1.w6.csv", row.names = 1)
LMS_GB_full_imp1_w9 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/LMS_2016_17_GlobalBalance_full_imp1.w9.csv", row.names = 1)

## 3.-------------------------------------------------------------------------------------------
# Linear models to look at how each community predicts performance separately
LMS_GB_full_imp1_w3.lm <- lm(numy ~ V1, LMS_GB_full_imp1_w3)
LMS_GB_full_imp1_w6.lm <- lm(numy ~ V1, LMS_GB_full_imp1_w6)
LMS_GB_full_imp1_w9.lm <- lm(numy ~ V1, LMS_GB_full_imp1_w9)
summary(LMS_GB_full_imp1_w3.lm)                               # R2adj = 0.8683
summary(LMS_GB_full_imp1_w6.lm)                               # R2adj = 0.8465
summary(LMS_GB_full_imp1_w9.lm)                               # R2adj = 0.8485

## 4.-------------------------------------------------------------------------------------------
# Build a df of the balances of the three weekly communities and see how well they predict yield performance
LMS_GB_full_imp1_w369 <- cbind.data.frame(LMS_GB_full_imp1_w3,
                                          LMS_GB_full_imp1_w6,
                                          LMS_GB_full_imp1_w9)
LMS_GB_full_imp1_w369 <- LMS_GB_full_imp1_w369[,c(1,2,4,6)]   # Remove the repeated performance measures
names(LMS_GB_full_imp1_w369) <- c("imp1","w3","w6","w9")      # Add useful column names
LMS_GB_full_imp1_w369.lm <- lm(imp1 ~ w3 + w6 + w9, 
                               LMS_GB_full_imp1_w369)         # Multiple regression of three weekly communities on performance
summary(LMS_GB_full_imp1_w369.lm)                             # R2adj = 0.9411
vif(LMS_GB_full_imp1_w369.lm)                                 # Surprisingly all VIFs are ~< 5.

## 5.-------------------------------------------------------------------------------------------
# Visualize the relationships between the weekly communities and yield-weighted performance
plot(LMS_GB_full_imp1_w369$imp1, LMS_GB_full_imp1_w369$w3)
plot(LMS_GB_full_imp1_w369$imp1, LMS_GB_full_imp1_w369$w6)
plot(LMS_GB_full_imp1_w369$imp1, LMS_GB_full_imp1_w369$w9)

# Visualize the relationship among the weekly communities
plot(LMS_GB_full_imp1_w369$w3, LMS_GB_full_imp1_w369$w6); cor.test(LMS_GB_full_imp1_w369$w3, LMS_GB_full_imp1_w369$w6)  # r = 0.861
plot(LMS_GB_full_imp1_w369$w3, LMS_GB_full_imp1_w369$w9); cor.test(LMS_GB_full_imp1_w369$w3, LMS_GB_full_imp1_w369$w9)  # r = 0.859
plot(LMS_GB_full_imp1_w369$w6, LMS_GB_full_imp1_w369$w9); cor.test(LMS_GB_full_imp1_w369$w6, LMS_GB_full_imp1_w369$w9)  # r = 0.869

# Calculate and visualize the variance for each weekly community balance
var(LMS_GB_full_imp1_w369$w3)
var(LMS_GB_full_imp1_w369$w6)
var(LMS_GB_full_imp1_w369$w9)
stripchart(LMS_GB_full_imp1_w3$V1, vertical = T, at = 3, xlim = c(1,10), ylim = c(-8,12), pch = 1)
stripchart(LMS_GB_full_imp1_w6$V1, vertical = T, at = 6, xlim = c(1,10), ylim = c(-8,12), pch = 1, add = T)
stripchart(LMS_GB_full_imp1_w9$V1, vertical = T, at = 9, xlim = c(1,10), ylim = c(-8,12), pch = 1, add = T)

## 6.-------------------------------------------------------------------------------------------
# Read in the selbal bal.tab output using the weekly datasets
Bal.Tab.full.w3.imp1 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/Bal.Tab.full.imp1.w3.csv", row.names = 1)
Bal.Tab.full.w6.imp1 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/Bal.Tab.full.imp1.w6.csv", row.names = 1)
Bal.Tab.full.w9.imp1 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/Bal.Tab.full.imp1.w9.csv", row.names = 1)
names(Bal.Tab.full.w3.imp1) <- c("taxa","group")
names(Bal.Tab.full.w6.imp1) <- c("taxa","group")
names(Bal.Tab.full.w9.imp1) <- c("taxa","group")

# Extract the taxa list from each of the above
sub.full.w3.imp1 <- as.character(Bal.Tab.full.w3.imp1$taxa)
sub.full.w6.imp1 <- as.character(Bal.Tab.full.w6.imp1$taxa)
sub.full.w9.imp1 <- as.character(Bal.Tab.full.w9.imp1$taxa)

# Check how many ASVs in each
dim(Bal.Tab.full.w3.imp1)   # 27 ASVs
dim(Bal.Tab.full.w6.imp1)   # 18 ASVs
dim(Bal.Tab.full.w9.imp1)   # 20 ASVs

## 7.-------------------------------------------------------------------------------------------
## Looking at the weeks 3, 6, 9 communities, we want to make sure these are REAL data
## Not like the last round where most of the ASVs were absent at a given site/year
## So plot relative abundance versus prevalence

## The 3 communities are all equally predictive (r2 ~ 0.84). Let's look at the predictive abilities
## in the alternate weeks. I.e., how predictive is the week 3 community in weeks 6 and 9?

# Read in the zero-replaced data that were used in the server selbal pipeline
asv.imp.full.w3.sq <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/asv.root.reps_dupes_merged.prev_and_abun_filt.sq.w3.csv", header = T, row.names = 1)
asv.imp.full.w6.sq <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/asv.root.reps_dupes_merged.prev_and_abun_filt.sq.w6.csv", header = T, row.names = 1)
asv.imp.full.w9.sq <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/asv.root.reps_dupes_merged.prev_and_abun_filt.sq.w9.csv", header = T, row.names = 1)

# Replace zeroes in the full dataset to use to fill in the missing taxa in other weeks.
# E.g., There are 27 ASVs in the week 3 community. 2 of those taxa are not present in the week 6 data, so let's
# fill those in from the full dataset.
can.asv.sq <- cmultRepl(can.asv, method = "SQ", output = "p-counts") # 6972 corrected values, range: 0.0000387 to 44,721,600

#### Ridiculous fact: the taxa list output from selbal.cv MUST be strings, not factors
Bal.Tab.full.w3.imp1$taxa <- as.character(Bal.Tab.full.w3.imp1$taxa)
Bal.Tab.full.w6.imp1$taxa <- as.character(Bal.Tab.full.w6.imp1$taxa)
Bal.Tab.full.w9.imp1$taxa <- as.character(Bal.Tab.full.w9.imp1$taxa)

# Subset the zero-replaced data by weekly relevant taxa
w3.comm.sq <- can.asv.sq[,Bal.Tab.full.w3.imp1$Taxa]
w6.comm.sq <- can.asv.sq[,Bal.Tab.full.w6.imp1$Taxa]
w9.comm.sq <- can.asv.sq[,Bal.Tab.full.w9.imp1$Taxa]

# Output the taxonomy for the 63 relevant ASVs and the abundances for each weekly community
write.csv(w3.comm.sq, "~/Desktop/w3.comm.27.csv")
write.csv(w6.comm.sq, "~/Desktop/w6.comm.18.csv")
write.csv(w9.comm.sq, "~/Desktop/w9.comm.20.csv")

## First check if we can recreate the week 3 balance using bal.value()
w3comm.w3.imp1.bal <- as.data.frame(bal.value(Bal.Tab.full.w3.imp1, log(asv.imp.full.w3.sq))) # Extract the balance
identical(rownames(LMS_GB_full_imp1_w3), rownames(w3comm.w3.imp1.bal))                        # Check to see if the row names match between the selbal and bal.value outputs

# There appears to be precision errors between the two balances (at the 15th decimal place), so the identical() function won't yield "TRUE"
identical(w3comm.w3.imp1.bal$`bal.value(Bal.Tab.full.w3.imp1, log(asv.imp.full.w3.sq))`, LMS_GB_full_imp1_w3$V1)  # FALSE

# But we can round the data to 14 decimal places and they will be identical
identical(round(w3comm.w3.imp1.bal$`bal.value(Bal.Tab.full.w3.imp1, log(asv.imp.full.w3.sq))`,14), round(LMS_GB_full_imp1_w3$V1,14))  # TRUE

## We can successfully recreate the global balance, so are comfortable using that to create
## balances in other weeks. c|:-)-{--<

## Get the taxa names from each community that are missing in other weeks
# Week 3 community missing in the week 6 sq-abundance df
w3.notin.w6 <- setdiff(Bal.Tab.full.w3.imp1$taxa, names(asv.imp.full.w6.sq))
w3.notin.w6.abu <- can.asv.sq[,w3.notin.w6]
asv.imp.full.w6.sq.w3 <- merge(asv.imp.full.w6.sq, w3.notin.w6.abu, by = 0, all.x = T)
rownames(asv.imp.full.w6.sq.w3) <- asv.imp.full.w6.sq.w3$Row.names
asv.imp.full.w6.sq.w3$Row.names <- NULL

# Week 3 community missing in the week 9 sq-abundance df
w3.notin.w9 <- setdiff(Bal.Tab.full.w3.imp1$taxa, names(asv.imp.full.w9.sq))
w3.notin.w9.abu <- can.asv.sq[,w3.notin.w9]
asv.imp.full.w9.sq.w3 <- merge(asv.imp.full.w9.sq, w3.notin.w9.abu, by = 0, all.x = T)
rownames(asv.imp.full.w9.sq.w3) <- asv.imp.full.w9.sq.w3$Row.names
asv.imp.full.w9.sq.w3$Row.names <- NULL

# Week 6 community missing in the week 3 sq-abundance df
w6.notin.w3 <- setdiff(Bal.Tab.full.w6.imp1$taxa, names(asv.imp.full.w3.sq))
w6.notin.w3.abu <- can.asv.sq[,w6.notin.w3]
asv.imp.full.w3.sq.w6 <- merge(asv.imp.full.w3.sq, w6.notin.w3.abu, by = 0, all.x = T)
rownames(asv.imp.full.w3.sq.w6) <- asv.imp.full.w3.sq.w6$Row.names
asv.imp.full.w3.sq.w6$Row.names <- NULL

# Week 6 community missing in the week 9 sq-abundance df - no missing taxa

## Now extract the balances
w3comm.w3.imp1.bal <- as.data.frame(bal.value(Bal.Tab.full.w3.imp1, log(asv.imp.full.w3.sq)))     # Week 3 community in week 3
w3comm.w6.imp1.bal <- as.data.frame(bal.value(Bal.Tab.full.w3.imp1, log(asv.imp.full.w6.sq.w3)))  # Week 3 community in week 6
w3comm.w9.imp1.bal <- as.data.frame(bal.value(Bal.Tab.full.w3.imp1, log(asv.imp.full.w9.sq.w3)))  # Week 3 community in week 9
w6comm.w3.imp1.bal <- as.data.frame(bal.value(Bal.Tab.full.w6.imp1, log(asv.imp.full.w3.sq.w6)))  # Week 6 community in week 3
w6comm.w6.imp1.bal <- as.data.frame(bal.value(Bal.Tab.full.w6.imp1, log(asv.imp.full.w6.sq)))     # Week 6 community in week 6
w6comm.w9.imp1.bal <- as.data.frame(bal.value(Bal.Tab.full.w6.imp1, log(asv.imp.full.w9.sq)))     # Week 6 community in week 9
w9comm.w3.imp1.bal <- as.data.frame(bal.value(Bal.Tab.full.w9.imp1, log(asv.imp.full.w3.sq)))     # Week 9 community in week 3
w9comm.w6.imp1.bal <- as.data.frame(bal.value(Bal.Tab.full.w9.imp1, log(asv.imp.full.w6.sq)))     # Week 6 community in week 6
w9comm.w9.imp1.bal <- as.data.frame(bal.value(Bal.Tab.full.w9.imp1, log(asv.imp.full.w9.sq)))     # Week 6 community in week 9

# Quick look at the row names to make sure they are in order
rownames(w3comm.w3.imp1.bal)
rownames(w3comm.w6.imp1.bal)
rownames(w3comm.w9.imp1.bal)

## 8.-------------------------------------------------------------------------------------------
# Build a sample df
Sample.ID <- c(rownames(LMS_GB_full_imp1_w3))
split.df <- unlist(strsplit(Sample.ID,"[.]"))
Site <- substr(Sample.ID,1,1)
Year <- substr(Sample.ID,3,6)
CanolaLine <- split.df[seq(3,length(split.df),5)]
SiteYearLine <- paste(Site, Year, CanolaLine, sep = ".")
can.sam.root.0 <- cbind.data.frame(Site = Site, Year = Year, CanolaLine = CanolaLine)
rownames(can.sam.root.0) <- SiteYearLine
can.sam.root.0$SiteYear <- paste(can.sam.root.0$Site, can.sam.root.0$Year, sep = ".")
can.sam.root.0 <- can.sam.root.0 %>%
  relocate(SiteYear, .before = Site)    # This just moves the "SiteYear" column to precede "Site"

# Build the df to use in the lms
comm.df <- cbind.data.frame(can.sam.root.0, 
                            imp1 = LMS_GB_full_imp1_w3$numy,
                            w3comm.w3 = w3comm.w3.imp1.bal[,1],
                            w3comm.w6 = w3comm.w6.imp1.bal[,1],
                            w3comm.w9 = w3comm.w9.imp1.bal[,1],
                            w6comm.w3 = w6comm.w3.imp1.bal[,1],
                            w6comm.w6 = w6comm.w6.imp1.bal[,1],
                            w6comm.w9 = w6comm.w9.imp1.bal[,1],
                            w9comm.w3 = w9comm.w3.imp1.bal[,1],
                            w9comm.w6 = w9comm.w6.imp1.bal[,1],
                            w9comm.w9 = w9comm.w9.imp1.bal[,1])

# Linear models to compare the global balances from the week 3, 6, and 9 communities in all time periods to performance
w3comm.w3.mod <- lm(imp1 ~ w3comm.w3, comm.df)
w3comm.w6.mod <- lm(imp1 ~ w3comm.w6, comm.df)
w3comm.w9.mod <- lm(imp1 ~ w3comm.w9, comm.df)
w6comm.w3.mod <- lm(imp1 ~ w6comm.w3, comm.df)
w6comm.w6.mod <- lm(imp1 ~ w6comm.w6, comm.df)
w6comm.w9.mod <- lm(imp1 ~ w6comm.w9, comm.df)
w9comm.w3.mod <- lm(imp1 ~ w9comm.w3, comm.df)
w9comm.w6.mod <- lm(imp1 ~ w9comm.w6, comm.df)
w9comm.w9.mod <- lm(imp1 ~ w9comm.w9, comm.df)

# Need predicted values to calculate RMSE
w3comm.w3.pred <- predict(w3comm.w3.mod, comm.df)
w3comm.w6.pred <- predict(w3comm.w6.mod, comm.df)
w3comm.w9.pred <- predict(w3comm.w9.mod, comm.df)
w6comm.w3.pred <- predict(w6comm.w3.mod, comm.df)
w6comm.w6.pred <- predict(w6comm.w6.mod, comm.df)
w6comm.w9.pred <- predict(w6comm.w9.mod, comm.df)
w9comm.w3.pred <- predict(w9comm.w3.mod, comm.df)
w9comm.w6.pred <- predict(w9comm.w6.mod, comm.df)
w9comm.w9.pred <- predict(w9comm.w9.mod, comm.df)

# Build a df of the results of the RMSE analyses
rmse.df <- cbind.data.frame(community = rep(c("W3","W6","W9"), each = 3),
                            week = rep(c(3,6,9), 3),
                            rmse = c(rmse(comm.df$imp1, w3comm.w3.pred),
                                     rmse(comm.df$imp1, w3comm.w6.pred),
                                     rmse(comm.df$imp1, w3comm.w9.pred),
                                     rmse(comm.df$imp1, w6comm.w3.pred),
                                     rmse(comm.df$imp1, w6comm.w6.pred), 
                                     rmse(comm.df$imp1, w6comm.w9.pred),
                                     rmse(comm.df$imp1, w9comm.w3.pred),
                                     rmse(comm.df$imp1, w9comm.w6.pred),
                                     rmse(comm.df$imp1, w9comm.w9.pred)),
                            r2adj = c(RsquareAdj(w3comm.w3.mod)[[2]],
                                      RsquareAdj(w3comm.w6.mod)[[2]],
                                      RsquareAdj(w3comm.w9.mod)[[2]],
                                      RsquareAdj(w6comm.w3.mod)[[2]],
                                      RsquareAdj(w6comm.w6.mod)[[2]],
                                      RsquareAdj(w6comm.w9.mod)[[2]],
                                      RsquareAdj(w9comm.w3.mod)[[2]],
                                      RsquareAdj(w9comm.w6.mod)[[2]],
                                      RsquareAdj(w9comm.w9.mod)[[2]]))

# Add a column that will be used to size points in the plot later on
rmse.df$vcex <- scales::rescale(abs(rmse.df$r2adj), to = c(0.5, 2.6))

# Subset to the three communities
rmse.df.w3 <- subset(rmse.df, community == "W3")
rmse.df.w6 <- subset(rmse.df, community == "W6")
rmse.df.w9 <- subset(rmse.df, community == "W9")

write.csv(rmse.df.w3, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/rmse.df.w3.csv")
write.csv(rmse.df.w6, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/rmse.df.w6.csv")
write.csv(rmse.df.w9, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/rmse.df.w9.csv")

## 9.-------------------------------------------------------------------------------------------
# Look at overlap in weeks 3, 6, and 9 winnowed communities
bal.w369.list <- Reduce(union, list(sub.full.w3.imp1, 
                                    sub.full.w6.imp1, 
                                    sub.full.w9.imp1))      # No ASVs are shared among all weeks
length(bal.w369.list)                                       # 63 unique ASVs altogether in the 3 weeks
sum(length(sub.full.w3.imp1), 
    length(sub.full.w6.imp1), 
    length(sub.full.w9.imp1))                               # 65 ASVs when adding all three communities together
                                                            # This means there are 65 - 63 = 2 non-unique ASVs
Reduce(intersect, list(sub.full.w3.imp1, sub.full.w6.imp1)) # 0 ASVs shared between 3 and 6
Reduce(intersect, list(sub.full.w3.imp1, sub.full.w9.imp1)) # 1 ASV shared between 3 and 9
Reduce(intersect, list(sub.full.w6.imp1, sub.full.w9.imp1)) # 1 ASV shared between 6 and 9

#### Generate Venn diagrams to look at overlap among three communities
w3.comm <- cbind.data.frame(sub.full.w3.imp1, rep(1,length(sub.full.w3.imp1)))
w6.comm <- cbind.data.frame(sub.full.w6.imp1, rep(1,length(sub.full.w6.imp1)))
w9.comm <- cbind.data.frame(sub.full.w9.imp1, rep(1,length(sub.full.w9.imp1)))
rownames(w3.comm) <- w3.comm$sub.full.w3.imp1
rownames(w6.comm) <- w6.comm$sub.full.w6.imp1
rownames(w9.comm) <- w9.comm$sub.full.w9.imp1
w3.comm$sub.full.w3.imp1 <- NULL
w6.comm$sub.full.w6.imp1 <- NULL
w9.comm$sub.full.w9.imp1 <- NULL
w369.comm <- merge(w3.comm, w6.comm, by = 0, all = T)
rownames(w369.comm) <- w369.comm$Row.names
w369.comm$Row.names <- NULL
w369.comm <- merge(w369.comm, w9.comm, by = 0, all = T)
names(w369.comm) <- c("asv", paste("w", c(3,6,9), sep = ""))
w369.comm[is.na(w369.comm)] <- 0
rownames(w369.comm) <- w369.comm$asv
w369.comm$asv <- NULL

## 10.------------------------------------------------------------------------

# Subset relative abundances for each weekly community
asv.relabu.union.w3 <- can.asv.rel0[,sub.full.w3.imp1]
asv.relabu.union.w6 <- can.asv.rel0[,sub.full.w6.imp1]
asv.relabu.union.w9 <- can.asv.rel0[,sub.full.w9.imp1]

# Subset absolute abundances for each weekly community
asv.absabu.union.w3 <- can.asv[,sub.full.w3.imp1]
asv.absabu.union.w6 <- can.asv[,sub.full.w6.imp1]
asv.absabu.union.w9 <- can.asv[,sub.full.w9.imp1]

# Make a sample df to subset with
Sample.ID <- c(rownames(asv.absabu.union.w3))
split.df <- unlist(strsplit(Sample.ID,"[.]"))
Site <- substr(Sample.ID,1,1)
Year <- substr(Sample.ID,3,6)
CanolaLine <- split.df[seq(3,length(split.df),5)]
Week <- split.df[seq(4,length(split.df),5)]
RootSoil <- split.df[seq(5,length(split.df),5)]
SiteYear <- paste(Site, Year, sep = ".")
SiteYearLineWeekRootSoil <- paste(Site, Year, CanolaLine, Week, RootSoil, sep = ".")
can.sam.root.1 <- cbind.data.frame(Site = Site, Year = Year, CanolaLine = CanolaLine, Week = Week, RootSoil = RootSoil, SiteYearLineWeekRootSoil = SiteYearLineWeekRootSoil, SiteYear = SiteYear)
rownames(can.sam.root.1) <- can.sam.root.1$SiteYearLineWeekRootSoil

## Relative abundance/sample merge
# Week 3 community
asv.relabu.union.w3.1 <- merge(can.sam.root.1, asv.relabu.union.w3, by = 0, all.y = T)  # Merge the sample df (above) to the relative abundance table
rownames(asv.relabu.union.w3.1) <- asv.relabu.union.w3.1$Row.names                      # Assign row names
asv.relabu.union.w3.1$Row.names <- NULL                                                 # Remove row names column
# Week 6 community
asv.relabu.union.w6.1 <- merge(can.sam.root.1, asv.relabu.union.w6, by = 0, all.y = T)  # Merge the sample df (above) to the relative abundance table
rownames(asv.relabu.union.w6.1) <- asv.relabu.union.w6.1$Row.names                      # Assign row names
asv.relabu.union.w6.1$Row.names <- NULL                                                 # Remove row names column
# Week 6 community
asv.relabu.union.w9.1 <- merge(can.sam.root.1, asv.relabu.union.w9, by = 0, all.y = T)  # Merge the sample df (above) to the relative abundance table
rownames(asv.relabu.union.w9.1) <- asv.relabu.union.w9.1$Row.names                      # Assign row names
asv.relabu.union.w9.1$Row.names <- NULL                                                 # Remove row names column

## Calculate weekly mean relative abundances per-taxa for each line within each site
# Week 3 community
relabu.union.w3.2 <- asv.relabu.union.w3.1[,-c(1,2,5,6)] %>%
  group_by(SiteYear, Week, CanolaLine) %>% 
  summarise_all(list(mean = mean))
relabu.union.w3.3 <- rowMeans(relabu.union.w3.2[-c(1:3)])                               # Mean abundance of winnowed community
# Week 6 community
relabu.union.w6.2 <- asv.relabu.union.w6.1[,-c(1,2,5,6)] %>%
  group_by(SiteYear, Week, CanolaLine) %>% 
  summarise_all(list(mean = mean))
relabu.union.w6.3 <- rowMeans(relabu.union.w6.2[-c(1:3)])                               # Mean abundance of winnowed community
# Week 9 community
relabu.union.w9.2 <- asv.relabu.union.w9.1[,-c(1,2,5,6)] %>%
  group_by(SiteYear, Week, CanolaLine) %>% 
  summarise_all(list(mean = mean))
relabu.union.w9.3 <- rowMeans(relabu.union.w9.2[-c(1:3)])                               # Mean abundance of winnowed community

## Create a df of mean relative abundances for the three communities
relabu.df <- cbind.data.frame(relabu.union.w3.2[,c(1:3)], w3 = relabu.union.w3.3, w6 = relabu.union.w6.3, w9 = relabu.union.w9.3)
relabu.df$Week <- as.integer(relabu.df$Week)

## Mean abundance and variance of each community in each week for each site
relabu.df <- relabu.df[,-3] %>%
  group_by(Week, SiteYear) %>% 
  summarise_all(list(mean = mean, var = var))

## Absolute abundance/sample merge
# Week 3 community
asv.absabu.union.w3.1 <- merge(can.sam.root.1, asv.absabu.union.w3, by = 0, all.y = T)  # Merge the sample df (above) to the relative abundance table
rownames(asv.absabu.union.w3.1) <- asv.absabu.union.w3.1$Row.names                      # Assign row names
asv.absabu.union.w3.1$Row.names <- NULL                                                 # Remove row names column
# Week 6 community
asv.absabu.union.w6.1 <- merge(can.sam.root.1, asv.absabu.union.w6, by = 0, all.y = T)  # Merge the sample df (above) to the relative abundance table
rownames(asv.absabu.union.w6.1) <- asv.absabu.union.w6.1$Row.names                      # Assign row names
asv.absabu.union.w6.1$Row.names <- NULL                                                 # Remove row names column
# Week 9 community
asv.absabu.union.w9.1 <- merge(can.sam.root.1, asv.absabu.union.w9, by = 0, all.y = T)  # Merge the sample df (above) to the relative abundance table
rownames(asv.absabu.union.w9.1) <- asv.absabu.union.w9.1$Row.names                      # Assign row names
asv.absabu.union.w9.1$Row.names <- NULL                                                 # Remove row names column

## Calculate weekly mean absolute abundances per-taxa for each line within each site
# Week 3 community
absabu.union.w3.2 <- asv.absabu.union.w3.1[,-c(1,2,5,6)] %>%
  group_by(SiteYear, Week, CanolaLine) %>% 
  summarise_all(list(mean = mean))
absabu.union.w3.3 <- rowMeans(absabu.union.w3.2[-c(1:3)])
# Week 6 community
absabu.union.w6.2 <- asv.absabu.union.w6.1[,-c(1,2,5,6)] %>%
  group_by(SiteYear, Week, CanolaLine) %>% 
  summarise_all(list(mean = mean))
absabu.union.w6.3 <- rowMeans(absabu.union.w6.2[-c(1:3)])
# Week 9 community
absabu.union.w9.2 <- asv.absabu.union.w9.1[,-c(1,2,5,6)] %>%
  group_by(SiteYear, Week, CanolaLine) %>% 
  summarise_all(list(mean = mean))
absabu.union.w9.3 <- rowMeans(absabu.union.w9.2[-c(1:3)])

## Create a df of mean absolute abundances for the three communities
absabu.df <- cbind.data.frame(absabu.union.w3.2[,c(1:3)], w3 = absabu.union.w3.3, w6 = absabu.union.w6.3, w9 = absabu.union.w9.3)
absabu.df$Week <- as.integer(absabu.df$Week)

## Mean abundance and variance of each community in each week for each site
absabu.df2 <- absabu.df[,-c(1,3)] %>%
  group_by(Week) %>% 
  summarise_all(list(mean = mean, var = var))

# Function for bootstrapped confidence intervals
my_boot = function(x, times=1000) {
  # Get column name from input object
  var = deparse(substitute(x))
  var = gsub("^\\.\\$","", var)
  # Bootstrap 95% CI
  cis = quantile(replicate(times, mean(sample(x, replace=TRUE))), probs=c(0.025,0.975))
  # Return data frame of results
  data.frame(var, n=length(x), mean=mean(x), lower.ci=cis[1], upper.ci=cis[2])
}

## Calculate means and CIs for absolute abundances in each week
set.seed(42)                      # Have reproducible results
# Week 3 community
absabu.df.boot.w3 <- absabu.df[,-c(1,3)] %>%
  group_by(Week) %>%
  do(my_boot(.$w3))
# Week 6 community
absabu.df.boot.w6 <- absabu.df[,-c(1,3)] %>%
  group_by(Week) %>%
  do(my_boot(.$w6))
# Week 9 community
absabu.df.boot.w9 <- absabu.df[,-c(1,3)] %>%
  group_by(Week) %>%
  do(my_boot(.$w9))

# Make a df of the means and CIs for each weekly community
absabu.df3 <- cbind.data.frame(absabu.df2, absabu.df.boot.w3[,c(3:6)], absabu.df.boot.w6[,c(3:6)], absabu.df.boot.w9[,c(3:6)])

# Add meaningful column names
names(absabu.df3) <- c("week", paste("w",c(3,6,9), rep(c("_mean","_var"), each = 3), sep = ""), paste("w",rep(c(3,6,9), each = 4), "_", names(absabu.df3)[8:19], sep = ""))

# Keep only the columns needed for plotting
absabu.df4 <- cbind.data.frame(week = absabu.df3$week, log(absabu.df3[,c(9:11,13:15,17:19)]))

# Subset to weeks 3, 6, and 9
absabu.df5 <- subset(absabu.df4, week == 3 | week == 6 | week == 9)

write.csv(absabu.df5, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/std.read.counts.csv")

# Make prevalence-abundance df
mean.abu <- as.data.frame(colMeans(can.asv.rel0))   # Mean relative abundance of each taxa across samples
prev <- as.data.frame(colSums(can.asv.rel0 != 0))   # Prevalence of each taxa across samples
prev.abu.df <- merge(mean.abu, prev, by = 0)        # Merge abundance and prevalence dfs from above
rownames(prev.abu.df) <- prev.abu.df$Row.names      # Add row names
prev.abu.df$Row.names <- NULL                       # Remove row names column from merge
names(prev.abu.df) <- c("abundance","prevalence")   # Meaningful column names
write.csv(prev.abu.df, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/pre.abu.df.csv")

# Subset to members of the weekly communities
prev.abu.df.w3 <- prev.abu.df[names(asv.relabu.union.w3),]
prev.abu.df.w6 <- prev.abu.df[names(asv.relabu.union.w6),]
prev.abu.df.w9 <- prev.abu.df[names(asv.relabu.union.w9),]

## 11.------------------------------------------------------------------------
#### Phylogenetic analyses

# Another way of thinking about the phylogenetic relatedness of species in a community is to ask ’how closely related are the average pair of species or individuals in a
# community?’, and relate the patterns we observe to what we’d expect under various
# null models of evolution and community assembly.
# The function mpd will calculate the mean pairwise distance (MPD) between all
# species in each community. Similarly, the mntd function calculates the mean nearest
# taxon distance (MNTD), the mean distance separating each species in the community
# from its closest relative. The mpd and mntd functions differs slightly from the pd
# function in that they take a distance matrix as input rather than a phylogeny object.
# A phylo object can be converted to a interspecific phylogenetic distance matrix using
# the cophenetic function. Since the mpd and mntd functions can use any distance
# matrix as input, we could easily calculate trait diversity measures by substituting a
# trait distance matrix for the phylogenetic distance matrix.

# Read in the tree and make sure ASV identifiers starting with a numeric have an "X" prepended
# root.tree <- phytools::read.newick("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/insertion-tree.nwk")
root.tree <- phyloseq::read_tree_greengenes("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/insertion-tree.nwk")
x <- root.tree$tip.label
x[grepl("^[[:digit:]]", x)] <- paste("X", x[grepl("^[[:digit:]]", x)], sep = "")
root.tree$tip.label <- x

# Subset the Aliivibio-standardized data set to the ASVs from each balance
w3.w3.comm <- can.asv[,Bal.Tab.full.w3.imp1$taxa]
w6.w6.comm <- can.asv[,Bal.Tab.full.w6.imp1$taxa]
w9.w9.comm <- can.asv[,Bal.Tab.full.w9.imp1$taxa]
w3.w3.comm2 <- can.asv[,unique(c(Bal.Tab.full.w3.imp1$taxa, Bal.Tab.full.w6.imp1$taxa, Bal.Tab.full.w9.imp1$taxa))]

# Now subset to include only the samples from each week (3,6,9)
w3.rows <- rownames(asv.imp.full.w3.sq)
w6.rows <- rownames(asv.imp.full.w6.sq)
w9.rows <- rownames(asv.imp.full.w9.sq)

can.asv.w3 <- can.asv[w3.rows,]
can.asv.w6 <- can.asv[w6.rows,]
can.asv.w9 <- can.asv[w9.rows,]
can.asv.63.w3 <- w3.w3.comm2[w3.rows,]
can.asv.63.w6 <- w3.w3.comm2[w6.rows,]
can.asv.63.w9 <- w3.w3.comm2[w9.rows,]

w3.w3.comm <- w3.w3.comm[w3.rows,]
w6.w6.comm <- w6.w6.comm[w6.rows,]
w9.w9.comm <- w9.w9.comm[w9.rows,]
w3.w3.comm.2 <- w3.w3.comm2[w3.rows,]
w6.w6.comm.2 <- w3.w3.comm2[w6.rows,]
w9.w9.comm.2 <- w3.w3.comm2[w9.rows,]
w369.comm.2 <- w3.w3.comm2[c(w3.rows, w6.rows, w9.rows),]

# Intra-comparisons

# Subset the tree to include only the ASVs of interest (this is to avoid errors due to the large tree size and computer memory)
# This is the error you would receive if the tree is too large:
# Error in double(nm * nm) : vector size cannot be NA 
# In addition: Warning message: In nm * nm : NAs produced by integer overflow

# Weeks 3, 6, and 9 tree
w369.tree <- drop.tip(root.tree, root.tree$tip.label[-match(unique(c(Bal.Tab.full.w3.imp1$taxa, Bal.Tab.full.w6.imp1$taxa, Bal.Tab.full.w9.imp1$taxa)), root.tree$tip.label)])

## Calculate empirical betaMNTD (https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/bNTI_Local_Machine.r)

match.phylo.otu.w369 <- match.phylo.data(as.phylo(w369.tree), t(w369.comm.2))

# Calculate empirical betaMNTD
beta.mntd.weighted.63 <- as.matrix(comdistnt(t(match.phylo.otu.w369$data),cophenetic(match.phylo.otu.w369$phy), abundance.weighted=T))

# Just a check, should = TRUE
identical(colnames(match.phylo.otu.w369$data),colnames(beta.mntd.weighted.63))
identical(colnames(match.phylo.otu.w369$data),rownames(beta.mntd.weighted.63))

write.csv(beta.mntd.weighted.63, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/beta.mntd.weighted.63.csv")

#### 12.------------------------------------------------------------------------
#### SEMs using flowering, genetic distance (everything w.r.t. NAM-79 - the lowest yielding line), and the three balances versus yield performance

## Read in the data
flowering <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/flowering_info.csv", header = T)
flowering$SiteYearLine <- paste(flowering$Site, flowering$Year, flowering$CanolaLine, sep = ".")  # Add column to group by

## Calculate means for each site-year-line and add rownames
flowering.df <- flowering %>%
  group_by(SiteYearLine) %>%
  summarise(w3.start.diff = mean(w3.start.diff),
            w6.start.diff = mean(w6.start.diff),
            w9.start.diff = mean(w9.start.diff), 
            w3.end.diff = mean(w3.end.diff),
            w6.end.diff = mean(w6.end.diff),
            w9.end.diff = mean(w9.end.diff),
            w369.tmp = mean(w369.tmp),
            w369.pre = mean(w369.pre),
            w3.tmp = mean(w3.tmp),
            w3.pre = mean(w3.pre),
            w6.tmp = mean(w6.tmp),
            w6.pre = mean(w6.pre),
            w9.tmp = mean(w9.tmp),
            w9.pre = mean(w9.pre))
flowering.df <- as.data.frame(flowering.df)
rownames(flowering.df) <- flowering.df$SiteYearLine

## Read in the file to extract canola genetic distance from (distances are relative to NAM-79 - the lowest yield-performing line)
can.bal.df <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/Modelling_GB.GD.FT.int.imp1.2.csv", header = T)[,-1]

## Merge flowering, distance, and balance dfs, and tidy up
# Extract site-year-line information for sample IDs and add to the distance df
split.df <- unlist(strsplit(can.bal.df$SampleID,"[.]"))
site.1 <- split.df[seq(1,length(split.df),5)]
year.1 <- split.df[seq(2,length(split.df),5)]
line.1 <- split.df[seq(3,length(split.df),5)]
distance.df <- cbind.data.frame(SampleID = paste(site.1, year.1, line.1, sep = "."),
                                Distance = can.bal.df$Distance)
rownames(distance.df) <- distance.df$SampleID

flowdist <- merge(flowering.df, distance.df, by = 0)
rownames(flowdist) <- flowdist$Row.names
flowdist <- flowdist[,-c(1,2)]
flowdist$SampleID <- NULL

# Extract site-year-line information for sample IDs and add to the balance df
balance.df <- cbind.data.frame(LMS_GB_full_imp1_w3,
                               LMS_GB_full_imp1_w6,
                               LMS_GB_full_imp1_w9)
split.df <- unlist(strsplit(rownames(balance.df),"[.]"))
site.1 <- split.df[seq(1,length(split.df),5)]
year.1 <- split.df[seq(2,length(split.df),5)]
line.1 <- split.df[seq(3,length(split.df),5)]
rownames(balance.df) <- paste(site.1, year.1, line.1, sep = ".")
balance.df$numy <- NULL
balance.df$numy <- NULL
names(balance.df) <- c("w3.bal","w6.bal","performance","w9.bal")

# Merge and tidy
flowdistbal.df <- merge(flowdist, balance.df, by = 0)
flowdistbal.df <- flowdistbal.df %>%
  relocate(performance, .before = w3.start.diff)
flowdistbal.df <- cbind.data.frame(flowdistbal.df[,c(1,2)], decostand(flowdistbal.df[,c(3:5,9:20)], "standardize"))

write.csv(flowdistbal.df, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/flowdistbal.df.csv")

# Check for collinearity
require(psych)
pairs.panels(flowdistbal.df[c("Distance",
                              "w3.start.diff","w6.start.diff","w9.start.diff",
                              "w3.bal","w6.bal","w9.bal",
                              "w369.tmp","w369.pre","w3.tmp","w3.pre","w6.tmp","w6.pre","w9.tmp","w9.pre")], lm = TRUE, stars = TRUE) # Diff and start diff correlated
flowdistbal.sub <- flowdistbal.df[,c(1:5,8:17)]

check_factorstructure(flowdistbal.sub[,-c(1,2)])
n_factors(flowdistbal.sub[,-c(1,2)])

plot(flowdistbal.df$w3.tmp, ylim = c(-1.5,1.5))
points(flowdistbal.df$w369.tmp, col = "blue")
plot(flowdistbal.df$w3.pre, ylim = c(-1.5,1.5))
points(flowdistbal.df$w369.pre, col = "blue")

# Use matched community metrics (e.g., w3 tmp and pre, w3 flowering diff, etc.)

# Week 3
mod.w3.full <- '
performance ~ Distance + w3.start.diff + w3.tmp + w3.pre + w3.bal
w3.bal ~ w3.tmp + w3.pre + Distance
w3.start.diff ~ w3.tmp + w3.pre + Distance
w3.tmp ~~ w3.pre'

fit.w3.full <- cfa(mod.w3.full, data = flowdistbal.df, fixed.x = FALSE)
summary(fit.w3.full, fit.measures = TRUE, standardized = T, rsquare = T)
MI <- modificationIndices(fit.w3.full)
subset(MI, mi>5)

# Week 6
mod.w6.full <- '
performance ~ Distance + w6.start.diff + w6.tmp + w6.pre + w6.bal
w6.bal ~ w6.tmp + w6.pre + Distance
w6.start.diff ~ w6.tmp + w6.pre + Distance
w6.tmp ~~ w6.pre'

fit.w6.full <- cfa(mod.w6.full, data = flowdistbal.df, fixed.x = FALSE)
summary(fit.w6.full, fit.measures = TRUE, standardized = T, rsquare = T)
MI <- modificationIndices(fit.w6.full)
subset(MI, mi>5)

# Week 9
mod.w9.full <- '
performance ~ Distance + w9.start.diff + w9.tmp + w9.pre + w9.bal
w9.bal ~ w9.tmp + w9.pre + Distance
w9.start.diff ~ w9.tmp + w9.pre + Distance
w9.tmp ~~ w9.pre'

fit.w9.full <- cfa(mod.w9.full, data = flowdistbal.df, fixed.x = FALSE)
summary(fit.w9.full, fit.measures = TRUE, standardized = T, rsquare = T)
MI <- modificationIndices(fit.w9.full)
subset(MI, mi>5)

## Check for non-linear relationships
summary(glm(performance ~ Distance + I(Distance^2), data = flowdistbal.df, family = gaussian))
summary(glm(performance ~ w3.start.diff + I(w3.start.diff^2), data = flowdistbal.df, family = gaussian))
summary(glm(performance ~ w6.start.diff + I(w6.start.diff^2), data = flowdistbal.df, family = gaussian))
summary(glm(performance ~ w9.start.diff + I(w9.start.diff^2), data = flowdistbal.df, family = gaussian))
summary(glm(performance ~ w3.bal + I(w3.bal^2), data = flowdistbal.df, family = gaussian))
summary(glm(performance ~ w6.bal + I(w6.bal^2), data = flowdistbal.df, family = gaussian))
summary(glm(performance ~ w9.bal + I(w9.bal^2), data = flowdistbal.df, family = gaussian))

# Check for collinearity
require(psych)
pairs.panels(flowdistbal.df[c("Distance","w3.start.diff","w6.start.diff","w9.start.diff","w3.bal","w6.bal","w9.bal")], lm = TRUE, stars = TRUE) # Diff and start diff correlated

# There is significant collinearity among the flowering times and balances. Let's just use week 3 for the coefficient comparison

## Model yield metrics as a function of genetic distance, flowering date differential, and global balance
can.yield.stab.m2 <- glm(performance ~ Distance + 
                           w3.start.diff +
                           w3.bal, data = flowdistbal.df, family = gaussian)       # Interaction b/w GB and diff
summary(can.yield.stab.m2)            # Distance: P = 0.8508; GlobalBalance: P < 2e-16; Differential: P = 0.2168; GB:diff: P = 0.0117
RsquareAdj(can.yield.stab.m2)         # R2adj = 0.924
can.yield.stab.m3 <- lm(performance ~ Distance + 
                           w3.start.diff +
                           w3.bal, data = flowdistbal.df)       # Interaction b/w GB and diff
summary(can.yield.stab.m3)
coef.df <- cbind.data.frame(coefs = coef(can.yield.stab.m3),
                            confint(can.yield.stab.m3))
coef.df <- coef.df[order(coef.df$coefs),]
coef.df <- coef.df[-which.max(coef.df$coefs),] # Remove the intercept
names(coef.df) <- c("coef","ci2p5","ci975")

af1 <- anova(can.yield.stab.m3)
afss1 <- af1$"Sum Sq"
print(cbind(af1,PctExp=afss1/sum(afss1)*100))
# Df      Sum Sq     Mean Sq     F value       Pr(>F)      PctExp
# Distance                  1  2.83736422  2.83736422  37.6248925 7.743629e-08  4.50375274
# GlobalBalance             1 55.15230021 55.15230021 731.3475476 6.101669e-35 87.54333366
# JD.w6.diff                1  0.05016368  0.05016368   0.6651958 4.180121e-01  0.07962488
# GlobalBalance:JD.w6.diff  1  0.51087043  0.51087043   6.7744017 1.167688e-02  0.81090545
# Residuals                59  4.44930146  0.07541189          NA           NA  7.06238327

cor.test(can.bal.df$Importance1, can.bal.df$GlobalBalance) # r = 0.9594

#### 13.------------------------------------------------------------------------
#### ASVs vs exudates

## Read in the data
can.env <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/canola_environmental.csv", row.names = 1, header = T)
can.tax.63 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/taxonomy63.csv", row.names = 1)

# Merge sample and abundance information and match to the exudates
can.sam.sub <- can.sam[can.sam$Week %in% c(3,6,9) & can.sam$Year == 2016,]
can.asv.sub <- can.asv[rownames(can.sam.sub), rownames(can.tax.63)]
can.asv.mat <- merge(can.sam.sub, can.asv.sub, by = 0)
rownames(can.asv.mat) <- can.asv.mat$Row.names
can.asv.mat$Row.names <- NULL
rownames(can.env) <- rownames(can.asv.mat)

# Hellinger transformation
can.asv.mat.hel <- decostand(can.asv.mat[,c(8:70)], method = "hellinger")

# RDA comparing exudates to community composition
can.rda <- rda(can.asv.mat.hel ~ ., data = can.env[,-c(1,2)])
set.seed(42); anova(can.rda, by = "terms")
set.seed(42); anova(can.rda, test = "F")

# Look at R2s, VIFs, and summary information
can.rda.R2 <- RsquareAdj(can.rda)$r.squared
can.rda.R2adj <- RsquareAdj(can.rda)$adj.r.squared
summary(can.rda)

