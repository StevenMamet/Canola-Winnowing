library(vegan)
library(cluster)
library(scales)
library(ape)
library(eulerr)
library(tidyverse)

rm(list = ls())

# Read in the sample df
can.sam <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/canola.sam.root.reps.merged.csv", header = T, row.names = 1)

# Add yield and line information to the sample df (was lost in sample_merge in step 1 R code)
can.metadata <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/canola_metadata.csv", header = T)
names(can.metadata)[2] <- "SiteYearPlot"
can.sample <- can.metadata[can.metadata$CanolaLine %in% c("NAM-0","NAM-5","NAM-13","NAM-14","NAM-17","NAM-23","NAM-30","NAM-32",
                                           "NAM-37","NAM-43","NAM-46","NAM-48","NAM-72","NAM-76","NAM-79","YN04-C1213"),]

## ******************
## CV calcs

# Create site-year metric for grouping calculations
can.sample$SiteYear <- as.factor(paste(can.sample$site.short, can.sample$year, sep = "."))

# Fill the missing value with the mean for that line
can.sample[is.na(can.sample$Yield.kg.ha),] # 5 missing values
can.sample[can.sample$SiteYearPlot == "L.2016.1228","Yield.kg.ha"] <- aggregate(Yield.kg.ha ~ SiteYear + CanolaLine, data = can.sample, mean)[130,3]  # L.2016, NAM-46
can.sample[can.sample$SiteYearPlot == "S.2017.1106","Yield.kg.ha"] <- aggregate(Yield.kg.ha ~ SiteYear + CanolaLine, data = can.sample, mean)[40,3]   # S.2017, NAM-14
can.sample[can.sample$SiteYearPlot == "S.2017.1113","Yield.kg.ha"] <- aggregate(Yield.kg.ha ~ SiteYear + CanolaLine, data = can.sample, mean)[96,3]   # S.2017, NAM-32
can.sample[can.sample$SiteYearPlot == "S.2017.1306","Yield.kg.ha"] <- aggregate(Yield.kg.ha ~ SiteYear + CanolaLine, data = can.sample, mean)[192,3]  # S.2017, NAM-76
can.sample[can.sample$SiteYearPlot == "S.2017.1313","Yield.kg.ha"] <- aggregate(Yield.kg.ha ~ SiteYear + CanolaLine, data = can.sample, mean)[164,3]  # S.2017, NAM-5

## Calculate the site-genotype mean and standard deviation to use for CV calcs
RCBD<-function(x,r,t)
{
  d=data.frame(x)
  r=as.numeric(r)
  t=as.numeric(t)
  CanolaLine=d[,1]
  replication=d[,2]
  response=d[,3]
  CanolaLine=as.factor(CanolaLine)
  replication=as.factor(replication)
  fit1=aov(response~CanolaLine+replication)
  result=summary(fit1)
  meantables=model.tables(fit1,type="means",se=T)
  sed=meantables$se$CanolaLine # standard error of difference of mean
  edf=(r-1)*(t-1)
  cd=qt(.975,edf)*sed # critical difference 
  list(result=result,meantables=meantables,sed=sed,cd=cd)
}

## Subset into 14 trials
sam.B.2015 <- subset(can.sample, SiteYear == "B.2015")
sam.B.2017 <- subset(can.sample, SiteYear == "B.2017")
sam.L.2014 <- subset(can.sample, SiteYear == "L.2014")
sam.L.2016 <- subset(can.sample, SiteYear == "L.2016")
sam.L.2017 <- subset(can.sample, SiteYear == "L.2017")
sam.M.2015 <- subset(can.sample, SiteYear == "M.2015")
sam.M.2017 <- subset(can.sample, SiteYear == "M.2017")
sam.N.2017 <- subset(can.sample, SiteYear == "N.2017")
sam.O1.2015 <- subset(can.sample, SiteYear == "O1.2015")
sam.O2.2015 <- subset(can.sample, SiteYear == "O2.2015")
sam.S.2015 <- subset(can.sample, SiteYear == "S.2015")
sam.S.2017 <- subset(can.sample, SiteYear == "S.2017")
sam.S1.2014 <- subset(can.sample, SiteYear == "S1.2014")
sam.S1.2015 <- subset(can.sample, SiteYear == "S1.2015")

## Calculate number of yield reps at each trial
sam.B.2015 <- sam.B.2015 %>% 
  arrange(SiteYear, CanolaLine, Yield.kg.ha) %>% 
  group_by(CanolaLine) %>% 
  mutate(rep = row_number())

sam.B.2017 <- sam.B.2017 %>% 
  arrange(SiteYear, CanolaLine, Yield.kg.ha) %>% 
  group_by(CanolaLine) %>% 
  mutate(rep = row_number())

sam.L.2014 <- sam.L.2014 %>% 
  arrange(SiteYear, CanolaLine, Yield.kg.ha) %>% 
  group_by(CanolaLine) %>% 
  mutate(rep = row_number())

sam.L.2016 <- sam.L.2016 %>% 
  arrange(SiteYear, CanolaLine, Yield.kg.ha) %>% 
  group_by(CanolaLine) %>% 
  mutate(rep = row_number())

sam.L.2017 <- sam.L.2017 %>% 
  arrange(SiteYear, CanolaLine, Yield.kg.ha) %>% 
  group_by(CanolaLine) %>% 
  mutate(rep = row_number())

sam.M.2015 <- sam.M.2015 %>% 
  arrange(SiteYear, CanolaLine, Yield.kg.ha) %>% 
  group_by(CanolaLine) %>% 
  mutate(rep = row_number())

sam.M.2017 <- sam.M.2017 %>% 
  arrange(SiteYear, CanolaLine, Yield.kg.ha) %>% 
  group_by(CanolaLine) %>% 
  mutate(rep = row_number())

sam.N.2017 <- sam.N.2017 %>% 
  arrange(SiteYear, CanolaLine, Yield.kg.ha) %>% 
  group_by(CanolaLine) %>% 
  mutate(rep = row_number())

sam.O1.2015 <- sam.O1.2015 %>% 
  arrange(SiteYear, CanolaLine, Yield.kg.ha) %>% 
  group_by(CanolaLine) %>% 
  mutate(rep = row_number())

sam.O2.2015 <- sam.O2.2015 %>% 
  arrange(SiteYear, CanolaLine, Yield.kg.ha) %>% 
  group_by(CanolaLine) %>% 
  mutate(rep = row_number())

sam.S.2015 <- sam.S.2015 %>% 
  arrange(SiteYear, CanolaLine, Yield.kg.ha) %>% 
  group_by(CanolaLine) %>% 
  mutate(rep = row_number())

sam.S.2017 <- sam.S.2017 %>% 
  arrange(SiteYear, CanolaLine, Yield.kg.ha) %>% 
  group_by(CanolaLine) %>% 
  mutate(rep = row_number())

sam.S1.2014 <- sam.S1.2014 %>% 
  arrange(SiteYear, CanolaLine, Yield.kg.ha) %>% 
  group_by(CanolaLine) %>% 
  mutate(rep = row_number())

sam.S1.2015 <- sam.S1.2015 %>% 
  arrange(SiteYear, CanolaLine, Yield.kg.ha) %>% 
  group_by(CanolaLine) %>% 
  mutate(rep = row_number())

# Create new data frames
sam.B.2015  <- cbind.data.frame(CanolaLine = sam.B.2015$CanolaLine, replication = sam.B.2015$rep, response = sam.B.2015$Yield.kg.ha)
sam.B.2017  <- cbind.data.frame(CanolaLine = sam.B.2017$CanolaLine, replication = sam.B.2017$rep, response = sam.B.2017$Yield.kg.ha)
sam.L.2014  <- cbind.data.frame(CanolaLine = sam.L.2014$CanolaLine, replication = sam.L.2014$rep, response = sam.L.2014$Yield.kg.ha)
sam.L.2016  <- cbind.data.frame(CanolaLine = sam.L.2016$CanolaLine, replication = sam.L.2016$rep, response = sam.L.2016$Yield.kg.ha)
sam.L.2017  <- cbind.data.frame(CanolaLine = sam.L.2017$CanolaLine, replication = sam.L.2017$rep, response = sam.L.2017$Yield.kg.ha)
sam.M.2015  <- cbind.data.frame(CanolaLine = sam.M.2015$CanolaLine, replication = sam.M.2015$rep, response = sam.M.2015$Yield.kg.ha)
sam.M.2017  <- cbind.data.frame(CanolaLine = sam.M.2017$CanolaLine, replication = sam.M.2017$rep, response = sam.M.2017$Yield.kg.ha)
sam.N.2017  <- cbind.data.frame(CanolaLine = sam.N.2017$CanolaLine, replication = sam.N.2017$rep, response = sam.N.2017$Yield.kg.ha)
sam.O1.2015 <- cbind.data.frame(CanolaLine = sam.O1.2015$CanolaLine, replication = sam.O1.2015$rep, response = sam.O1.2015$Yield.kg.ha)
sam.O2.2015 <- cbind.data.frame(CanolaLine = sam.O2.2015$CanolaLine, replication = sam.O2.2015$rep, response = sam.O2.2015$Yield.kg.ha)
sam.S.2015  <- cbind.data.frame(CanolaLine = sam.S.2015$CanolaLine, replication = sam.S.2015$rep, response = sam.S.2015$Yield.kg.ha)
sam.S.2017  <- cbind.data.frame(CanolaLine = sam.S.2017$CanolaLine, replication = sam.S.2017$rep, response = sam.S.2017$Yield.kg.ha)
sam.S1.2014 <- cbind.data.frame(CanolaLine = sam.S1.2014$CanolaLine, replication = sam.S1.2014$rep, response = sam.S1.2014$Yield.kg.ha)
sam.S1.2015 <- cbind.data.frame(CanolaLine = sam.S1.2015$CanolaLine, replication = sam.S1.2015$rep, response = sam.S1.2015$Yield.kg.ha)

# Use RCBD code to calculate site-genotype mean and standard deviation to use for CV calcs.
# The first number is the number of replicates, the second is number of genotypes.
sam.B.2015.rcbd <- RCBD(sam.B.2015, max(sam.B.2015$replication), nlevels(droplevels(as.factor(sam.B.2015$CanolaLine))))
sam.B.2017.rcbd <- RCBD(sam.B.2017, max(sam.B.2017$replication), nlevels(droplevels(as.factor(sam.B.2017$CanolaLine))))
sam.L.2014.rcbd <- RCBD(sam.L.2014, max(sam.L.2014$replication), nlevels(droplevels(as.factor(sam.L.2014$CanolaLine))))
sam.L.2016.rcbd <- RCBD(sam.L.2016, max(sam.L.2016$replication), nlevels(droplevels(as.factor(sam.L.2016$CanolaLine))))
sam.L.2017.rcbd <- RCBD(sam.L.2017, max(sam.L.2017$replication), nlevels(droplevels(as.factor(sam.L.2017$CanolaLine))))
sam.M.2015.rcbd <- RCBD(sam.M.2015, max(sam.M.2015$replication), nlevels(droplevels(as.factor(sam.M.2015$CanolaLine))))
sam.M.2017.rcbd <- RCBD(sam.M.2017, max(sam.M.2017$replication), nlevels(droplevels(as.factor(sam.M.2017$CanolaLine))))
sam.N.2017.rcbd <- RCBD(sam.N.2017, max(sam.N.2017$replication), nlevels(droplevels(as.factor(sam.N.2017$CanolaLine))))
sam.O1.2015.rcbd <- RCBD(sam.O1.2015, max(sam.O1.2015$replication), nlevels(droplevels(as.factor(sam.O1.2015$CanolaLine))))
sam.O2.2015.rcbd <- RCBD(sam.O2.2015, max(sam.O2.2015$replication), nlevels(droplevels(as.factor(sam.O2.2015$CanolaLine))))
sam.S.2015.rcbd <- RCBD(sam.S.2015, max(sam.S.2015$replication), nlevels(droplevels(as.factor(sam.S.2015$CanolaLine))))
sam.S.2017.rcbd <- RCBD(sam.S.2017, max(sam.S.2017$replication), nlevels(droplevels(as.factor(sam.S.2017$CanolaLine))))
sam.S1.2014.rcbd <- RCBD(sam.S1.2014, max(sam.S1.2014$replication), nlevels(droplevels(as.factor(sam.S1.2014$CanolaLine))))
sam.S1.2015.rcbd <- RCBD(sam.S1.2015, max(sam.S1.2015$replication), nlevels(droplevels(as.factor(sam.S1.2015$CanolaLine))))

# Make a yield data frame. Not all trials had all 16 genotypes of interest. Here I still include
# all 16 for each site, but codify an NA in the missing genotypes place.
can.agg <- rbind.data.frame(data.frame(Yield.kg.ha = c(sam.B.2015.rcbd$meantables$tables$CanolaLine[1:15], NA)),
                            data.frame(Yield.kg.ha = c(sam.B.2017.rcbd$meantables$tables$CanolaLine[1:15], NA)),
                            data.frame(Yield.kg.ha = c(sam.L.2014.rcbd$meantables$tables$CanolaLine[1:10], NA, 
                                                       sam.L.2014.rcbd$meantables$tables$CanolaLine[11:14], NA)),
                            data.frame(Yield.kg.ha = sam.L.2016.rcbd$meantables$tables$CanolaLine[1:16]),
                            data.frame(Yield.kg.ha = sam.L.2017.rcbd$meantables$tables$CanolaLine[1:16]),
                            data.frame(Yield.kg.ha = c(sam.M.2015.rcbd$meantables$tables$CanolaLine[1:15], NA)),
                            data.frame(Yield.kg.ha = sam.M.2017.rcbd$meantables$tables$CanolaLine[1:16]),
                            data.frame(Yield.kg.ha = c(sam.N.2017.rcbd$meantables$tables$CanolaLine[1:15], NA)),
                            data.frame(Yield.kg.ha = c(sam.O1.2015.rcbd$meantables$tables$CanolaLine[1:15], NA)),
                            data.frame(Yield.kg.ha = c(sam.O2.2015.rcbd$meantables$tables$CanolaLine[1:15], NA)),
                            data.frame(Yield.kg.ha = c(sam.S.2015.rcbd$meantables$tables$CanolaLine[1:15], NA)),
                            data.frame(Yield.kg.ha = sam.S.2017.rcbd$meantables$tables$CanolaLine[1:16]),
                            data.frame(Yield.kg.ha = c(sam.S1.2014.rcbd$meantables$tables$CanolaLine[1:10], NA, 
                                                       sam.S1.2014.rcbd$meantables$tables$CanolaLine[11:14], NA)),
                            data.frame(Yield.kg.ha = c(sam.S1.2015.rcbd$meantables$tables$CanolaLine[1:15], NA)))
can.agg$SiteYear <- rep(c("B.2015","B.2017","L.2014","L.2016","L.2017","M.2015","M.2017",
                          "N.2017","O1.2015","O2.2015","S.2015","S.2017","S1.2014","S1.2015"), each = 16)
can.agg$CanolaLine <- rep(c(paste("NAM-",c(0,13,14,17,23,30,32,37,43,46,48,5,72,76,79),sep = ""), "YN04-C1213"), 14)

## ******************
## Yield stability (S2 & S) calcs
## Following Lin et al. 1986. Contribution no. 1-770 from the Engineering and Statistical Res. Centre, Res. Branch, Agric. Canada, Ottawa, Canada

## Calculate the mean yield for each genotype (mi)
# Add column of rep numbers
can.sam2 <- can.sample %>% 
  arrange(CanolaLine, Yield.kg.ha) %>% 
  group_by(CanolaLine) %>% 
  mutate(rep = row_number())

# Make new df to use in function to calculate genotype means in an RCBD framework (Jeelani et al., 2018. IJISS)
can.sam3 <- droplevels(cbind.data.frame(CanolaLine = can.sam2$CanolaLine, replication = can.sam2$rep, response = can.sam2$Yield.kg.ha))
summary(can.sam3)

# Will have three yield measures in the end:
# 1. Per genotype (mi)
# 2. Per genotype within a trial (Rij)
# 3. Per plot (Yield.kg.ha)

can.sam.rcbd <- RCBD(can.sam3, max(can.sam3$replication),nlevels(as.factor(can.sam3$CanolaLine)))
can.agg2 <- data.frame(CanolaLine = names(can.sam.rcbd$meantables$tables$CanolaLine),
                       Yield.kg.ha = can.sam.rcbd$meantables$tables$CanolaLine[1:16])

# Create merging indices and merge
can.agg$SiteYearLine <- paste(can.agg$SiteYear, can.agg$CanolaLine, sep = ".")
can.sample$SiteYearLine <- paste(can.sample$site.short, can.sample$year, can.sample$CanolaLine, sep = ".")
can.sample2 <- merge(can.sample, can.agg[-c(2,3)], by = "SiteYearLine", all.x = TRUE)
can.sample3 <- merge(can.sample2, can.agg2, by = "CanolaLine", all.x = TRUE)

# Rename to something meaningful
names(can.sample3)[9] <- "Yield.kg.ha"
names(can.sample3)[11] <- "Rij"
names(can.sample3)[12] <- "mi"

# Set the number of environments - 1 (i.e., 14 - 1)
can.sample3$q_1 <- nlevels(can.sample3$SiteYear) - 1

# Intermediate calcs for s2
can.sample3$s2pre <- (can.sample3$Rij - can.sample3$mi)^2
can.agg3 <- aggregate(s2pre ~ SiteYear + CanolaLine, data = can.sample3, sum)
can.agg3$SiteYearLine <- paste(can.agg3$SiteYear, can.agg3$CanolaLine, sep = ".")
can.sample4 <- merge(can.sample3, can.agg3[,-1], by = "SiteYearLine", all.x = TRUE)
can.sample4[15] <- NULL
names(can.sample4)[14] <- "s2pre"
names(can.sample4)[15] <- "sum.s2pre"

# Calculate s2 and sqrt(s2)
can.sample4$s2 <- can.sample4$sum.s2pre / can.sample4$q_1
can.sample4$s <- sqrt(can.sample4$s2)
nlevels(can.sample4$SiteYear)

# Rescale stability be reverse values (i.e., high values indicate high stability rather than high variance)
can.sample5 <- can.sample4
can.sample5$s.rev <- scales::rescale(can.sample5$s, to = c(max(can.sample5$s),min(can.sample5$s)))

# Subset to only the 4 sites that we have microbiome data for
can.sample6 <- subset(can.sample5, SiteYear == "L.2016" | SiteYear == "L.2017" | SiteYear == "M.2017" | SiteYear == "S.2017")

# Calculate mean yield, stability
yield.mean <- mean(can.sample5$Yield.kg.ha)
stab.mean <- mean(can.sample5$s.rev)
can.sample5 <- droplevels(can.sample5)

write.csv(can.sample5, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/canola.2016_2017.plot_level_lines.csv")
write.csv(can.sample6, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/canola.2016_2017.plot_level_lines_microbiome.csv")

# Calculate means for each genotype within each trial
can.sample7 <- can.sample6 %>% 
  group_by(SiteYearLine) %>% 
  summarise(CanolaLine = CanolaLine.x[1], Year = year[1], Site = site.short[1], SiteYear = SiteYear[1], s2 = s2[1], s = s[1], s.rev = s.rev[1], Yield=mean(Yield.kg.ha))

## Here we will calculate our yield performance metrics

# Create columns to populate
can.sample7$s.std <- NA
can.sample7$yield.std <- NA

# Scale and centre the yield metrics
can.sample7$SiteYear <- droplevels(can.sample7$SiteYear)
can.sample7$s.std[can.sample7$SiteYear == "L.2016"] <- scale(can.sample7$s.rev[can.sample7$SiteYear == "L.2016"], scale = T)
can.sample7$s.std[can.sample7$SiteYear == "L.2017"] <- scale(can.sample7$s.rev[can.sample7$SiteYear == "L.2017"], scale = T)
can.sample7$s.std[can.sample7$SiteYear == "M.2017"] <- scale(can.sample7$s.rev[can.sample7$SiteYear == "M.2017"], scale = T)
can.sample7$s.std[can.sample7$SiteYear == "S.2017"] <- scale(can.sample7$s.rev[can.sample7$SiteYear == "S.2017"], scale = T)

can.sample7$yield.std[can.sample7$SiteYear == "L.2016"] <- scale(can.sample7$Yield[can.sample7$SiteYear == "L.2016"], scale = T)
can.sample7$yield.std[can.sample7$SiteYear == "L.2017"] <- scale(can.sample7$Yield[can.sample7$SiteYear == "L.2017"], scale = T)
can.sample7$yield.std[can.sample7$SiteYear == "M.2017"] <- scale(can.sample7$Yield[can.sample7$SiteYear == "M.2017"], scale = T)
can.sample7$yield.std[can.sample7$SiteYear == "S.2017"] <- scale(can.sample7$Yield[can.sample7$SiteYear == "S.2017"], scale = T)

can.sample7 <- as.data.frame(can.sample7)
can.sample7$s.std <- scales::rescale(can.sample7$s.std, to = c(0, 1))
can.sample7$yield.std <- scales::rescale(can.sample7$yield.std, to = c(0, 1))

plot(can.sample7$yield.std, can.sample7$s.std, pch = 16, col = scales::alpha("blue",0.4))

# Calculate a quadrant coefficient
can.sample7$gamma <- NA
can.sample7$gamma[can.sample7$yield.std > 0.5 & can.sample7$s.std > 0.5] <- (1/2)
can.sample7$gamma[can.sample7$yield.std < 0.5 & can.sample7$s.std < 0.5] <- (1/6)
can.sample7$gamma[is.na(can.sample7$gamma)] <- (1/3)
sum(is.na(can.sample7$gamma))   # Check for any NAs

# Calculate a yield/s weighting
can.sample7$alpha1 <- 0.75
can.sample7$alpha2 <- 0.25
can.sample7$alpha3 <- 0.5

# Calculate metrics of importance: 1. Yield-weighted. 2. Stability-weighted. 3. Unweighted.
can.sample7$importance1 <- sqrt((can.sample7$gamma * can.sample7$alpha1 * can.sample7$yield.std)^2 + (can.sample7$gamma * (1-can.sample7$alpha1) * can.sample7$s.std)^2)
can.sample7$importance2 <- sqrt((can.sample7$gamma * can.sample7$alpha2 * can.sample7$yield.std)^2 + (can.sample7$gamma * (1-can.sample7$alpha2) * can.sample7$s.std)^2)
can.sample7$importance3 <- sqrt((can.sample7$gamma * can.sample7$alpha3 * can.sample7$yield.std)^2 + (can.sample7$gamma * (1-can.sample7$alpha3) * can.sample7$s.std)^2)
can.sample7$SiteYearLine <- as.factor(can.sample7$SiteYearLine)

write.csv(can.sample7, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/canola_yield_performance.csv")

# Make the yield performance df
can.sample8 <- can.sample7 %>% 
  group_by(SiteYearLine) %>% 
  summarise(CanolaLine = CanolaLine[1], Year = Year[1], Site = Site[1], SiteYear = SiteYear[1], s2 = s2[1], s = s[1], s.std = s.std[1], s.rev = s.rev[1], 
            Yield=mean(Yield), yield.std = yield.std[1], importance1 = importance1[1], importance2 = importance2[1], importance3 = importance3[1])

# Merge the sample information with the performance df, curate, and name the columns appropriately
can.sam4 <- merge(can.sam, can.sample8, by = "SiteYearLine", all.y = TRUE)
range(can.sam4$s.std)                   # Should be [0,1]
can.sam4 <- can.sam4[,c(1:5,7,11:20)]   # Only keep the relevant columns
names(can.sam4) <- c("SiteYearLine","Site","Year","CanolaLine","Week","SampleID","SiteYear",
                     "s2","s","s.std","s.rev","Yield","yield.std",
                     "importance1","importance2","importance3")

write.csv(can.sam4, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/canola.sam.root.2.csv")
