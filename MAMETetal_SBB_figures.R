library(scales)         # Transparency in plotting and rescaling data
library(phyloseq)       # Working with phylogenetic trees
library(ape)            # Working with phylogenetic trees
library(ggtree)         # Plotting phylogenetic trees
library(wesanderson)    # Beautiful color palettes
library(vegan)          # Used for decostand
library(tidyverse)      # Tidyverse is lovely
library(dplyr)          

rm(list = ls())

####-------------------------------------------------------------------------------------------
#### Figure 2

## Histograms
## Read in the data
step0.hist <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/step0.taxa.csv", row.names = 1)
step1.hist <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/step1.taxa.csv", row.names = 1)
step2.hist <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/step2.taxa.csv", row.names = 1)
step3.hist <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/step3.taxa.csv", row.names = 1)
step4.hist <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/step4.taxa.csv", row.names = 1)
step5.hist <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/step5.taxa.csv", row.names = 1)

# 1. Noise
jpeg("hist01.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(4.75,4.75,1,1))
hist(log(as.matrix(step0.hist+1)), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000), main = "", ylab = "", xlab = "", col = "grey90")
par(new = T)
step1.hist.out <- hist(log(as.matrix(step1.hist+1)), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000), main = "", ylab = "", xlab = "", col = "#b94663")
dev.off()

# 2. Trace
jpeg("hist02.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(4.75,4.75,1,1))
hist(log(as.matrix(step1.hist+1)), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000), main = "", ylab = "", xlab = "", col = "#b94663")
par(new = T)
step2.hist.out <- hist(log(as.matrix(step2.hist+1)), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000), main = "", ylab = "", xlab = "", col = "#bc7d39")
dev.off()

# 3. Majority
jpeg("hist03.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(4.75,4.75,1,1))
hist(log(as.matrix(step2.hist+1)), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000), main = "", ylab = "", xlab = "", col = "#bc7d39")
par(new = T)
step3.hist.out <- hist(log(as.matrix(step3.hist+1)), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000), main = "", ylab = "", xlab = "", col = "#6fac5d")
dev.off()

# 4. Merge
jpeg("hist04.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(4.75,4.75,1,1))
hist(log(as.matrix(step3.hist+1)), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000), main = "", ylab = "", xlab = "", col = "#6fac5d")
par(new = T)
step4.hist.out <- hist(log(as.matrix(step4.hist+1)), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000), main = "", ylab = "", xlab = "", col = "#677ad1")
dev.off()

# 5. Predictive
jpeg("hist05.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(4.75,4.75,1,1))
hist(log(as.matrix(step4.hist+1)), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000), main = "", ylab = "", xlab = "", col = "#677ad1")
par(new = T)
step5.hist.out <- hist(log(as.matrix(step5.hist+1)), freq = TRUE, breaks = seq(0,20,0.5), xlim = c(0,20), ylim = c(0,50000), main = "", ylab = "", xlab = "", col = "#9750a1")
dev.off()

## Scatterplots
## Build relative abundance-prevalence df
can.asv.rel <- decostand(t(step0.hist), method = "total")*100
abu.prev.df <- cbind.data.frame(rel.abundance = colMeans(can.asv.rel), prevalence = colSums(can.asv.rel != 0), read.abundance = colSums(t(step0.hist)))

## Build the dfs to use for plotting
step0.df <- abu.prev.df
step1.df <- step0.df[rownames(step1.hist),]
step2.df <- step0.df[rownames(step2.hist),]
step3.df <- step0.df[rownames(step3.hist),]
step4.df <- step0.df[rownames(step4.hist),]
step5.df <- step0.df[colnames(step5.hist),]

## Create dfs for the composite samples (steps 4 and 5 in SBB figure 2)
# Step 4
step4.prev <- rowSums(step4.hist != 0)
step4.absabu <- rowSums(step4.hist)
step4.df1 <- cbind.data.frame(rel.abundance = NA, prevalence = rowSums(step4.hist != 0), read.abundance = rowSums(step4.hist))
step4.df1 <- merge(step4.df, step4.df1, by = 0)

# Step 5
step5.prev <- colSums(step5.hist != 0)
step5.absabu <- colSums(step5.hist)
step5.df1 <- cbind.data.frame(rel.abundance = NA, prevalence = colSums(step5.hist != 0), read.abundance = colSums(step5.hist))
step5.df1 <- merge(step5.df, step5.df1, by = 0)

## Step 1
jpeg("scatter01.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(4.75,4.75,1,1))
plot(rel.abundance ~ prevalence, log = "xy", yaxt = "n", pch = 16, cex = 0.7, axes = F, ylab = "", xlab = "", 
     col = scales::alpha("grey70", 0.5), data = step0.df[rowSums(step0.df) != 0,])
box("plot", lwd = 2)
points(step1.df$prevalence, step1.df$rel.abundance, cex = 0.75, pch = 21, bg = scales::alpha("#b94663",0.75), 
       col = "black", lwd = 0.75)
axis(1, at = 1*2^(0:10), labels = rep("",11),lwd = 2, cex.axis = 2)
axis(1, at = 1*2^(0:10), lwd = 0, cex.axis = 2.1, line = 0.6)
axis(2, 1*2^(0:10), labels = rep("",11), lwd = 2, cex.axis = 2)
axis(2, at = 1*2^(0:10), lwd = 0, cex.axis = 2.1, line = 0)
dev.off()

## Step 2
jpeg("scatter02.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(4.75,4.75,1,1))
plot(rel.abundance ~ prevalence, log = "xy", yaxt = "n", pch = 16, cex = 0.7, axes = F, ylab = "", xlab = "", 
     col = scales::alpha("white", 0.5), data = step0.df[rowSums(step0.df) != 0,])
points(rel.abundance ~ prevalence, yaxt = "n", pch = 16, cex = 0.7, col = scales::alpha("#b94663", 0.5), data = step1.df)
box("plot", lwd = 2)
points(step2.df$prevalence, step2.df$rel.abundance, cex = 0.75, pch = 21, bg = scales::alpha("#bc7d39",0.75), col = "black", lwd = 0.75)
axis(1, at = 1*2^(0:10), labels = rep("",11),lwd = 2, cex.axis = 2)
axis(1, at = 1*2^(0:10), lwd = 0, cex.axis = 2.1, line = 0.6)
axis(2, 1*2^(0:10), labels = rep("",11), lwd = 2, cex.axis = 2)
axis(2, at = 1*2^(0:10), lwd = 0, cex.axis = 2.1, line = 0)
dev.off()

## Step 3
jpeg("scatter03.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(4.75,4.75,1,1))
plot(rel.abundance ~ prevalence, log = "xy", yaxt = "n", pch = 16, cex = 0.7, axes = F, ylab = "", xlab = "", 
     col = scales::alpha("white", 0.5), data = step0.df[rowSums(step0.df) != 0,])
points(rel.abundance ~ prevalence, pch = 16, cex = 0.7, col = scales::alpha("#bc7d39", 0.5), data = step2.df)
box("plot", lwd = 2)
points(step3.df$prevalence, step3.df$rel.abundance, cex = 0.75, pch = 21, bg = scales::alpha("#6fac5d",0.75), col = "black", lwd = 0.75)
axis(1, at = 1*2^(0:10), labels = rep("",11),lwd = 2, cex.axis = 2)
axis(1, at = 1*2^(0:10), lwd = 0, cex.axis = 2.1, line = 0.6)
axis(2, 1*2^(0:10), labels = rep("",11), lwd = 2, cex.axis = 2)
axis(2, at = 1*2^(0:10), lwd = 0, cex.axis = 2.1, line = 0)
dev.off()

## Step 4
jpeg("scatter04.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(4.75,4.75,1,1))
plot(rel.abundance ~ prevalence, log = "xy", yaxt = "n", pch = 16, cex = 0.7, axes = F, ylab = "", xlab = "", 
     col = scales::alpha("white", 0.5), data = step0.df[rowSums(step0.df) != 0,])
points(rel.abundance ~ prevalence, pch = 16, cex = 0.7, col = scales::alpha("#6fac5d", 0.5), data = step3.df)
box("plot", lwd = 2)
points(step4.df1$rel.abundance.x ~ step4.df1$prevalence.y, cex = 0.75, pch = 21, bg = scales::alpha("#677ad1",0.75), col = "black", lwd = 0.75)
axis(1, at = 1*2^(0:10), labels = rep("",11),lwd = 2, cex.axis = 2)
axis(1, at = 1*2^(0:10), lwd = 0, cex.axis = 2.1, line = 0.6)
axis(2, 1*2^(0:10), labels = rep("",11), lwd = 2, cex.axis = 2)
axis(2, at = 1*2^(0:10), lwd = 0, cex.axis = 2.1, line = 0)
dev.off()

## Step 5
jpeg("scatter05.jpg", width = 4.25, height = 4, units = "in", res = 300)
par(mar = c(4.75,4.75,1,1))
plot(rel.abundance ~ prevalence, log = "xy", yaxt = "n", pch = 16, cex = 0.7, axes = F, ylab = "", xlab = "", 
     col = scales::alpha("white", 0.5), data = step0.df[rowSums(step0.df) != 0,])
points(rel.abundance.x ~ prevalence.y, pch = 16, cex = 0.7, col = scales::alpha("#677ad1", 0.5), data = step4.df1)
box("plot", lwd = 2)
points(step5.df1$rel.abundance.x ~ step5.df1$prevalence.y, cex = 0.75, pch = 21, bg = scales::alpha("#9750a1",0.75), col = "black", lwd = 0.75)
axis(1, at = 1*2^(0:10), labels = rep("",11),lwd = 2, cex.axis = 2)
axis(1, at = 1*2^(0:10), lwd = 0, cex.axis = 2.1, line = 0.6)
axis(2, 1*2^(0:10), labels = rep("",11), lwd = 2, cex.axis = 2)
axis(2, at = 1*2^(0:10), lwd = 0, cex.axis = 2.1, line = 0)
dev.off()

####-------------------------------------------------------------------------------------------
#### Figure 3a: Full yield versus yield variance for 14 trials

fig03a1.df <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/canola.2016_2017.plot_level_lines.csv")
fig03a2.df <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/canola.2016_2017.plot_level_lines_microbiome.csv")

jpeg("Fig03a.jpg", width = 4.6, height = 4, units = "in", res = 300)
par(mar = c(3,3,0.1,0.1))
plot.new(); plot.window(xlim = c(0,5000), ylim = c(0,800))
points(fig03a1.df$Yield.kg.ha, fig03a1.df$s, pch = 16, col = scales::alpha("gray", 0.6))
points(fig03a2.df$Yield.kg.ha, fig03a2.df$s, pch = 16, col = scales::alpha("blue", 0.6))
box()
axis(side = 1, at = seq(0,5000,1000), labels = F, cex.axis = 1.5)
axis(side = 2, at = seq(0,800,200), labels = F, cex.axis = 1.5)
dev.off()


####-------------------------------------------------------------------------------------------
#### Figure 3b: Unweighted yield performance

# Read in the data
can.sample7 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/canola_yield_performance.csv")

jpeg("Fig03b.jpg", width = 4.6, height = 4, units = "in", res = 300)
par(mar = c(3,3,0.1,0.1))
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
abline(h = 0.5, v = 0.5, lty = 3, col = "gray25")
points(can.sample7$yield.std, can.sample7$s.std, pch = 16, cex = can.sample7$importance1*10, col = scales::alpha("blue", 0.6), )
box()
axis(side = 1, at = seq(0,1,0.25), labels = F, cex.axis = 1.5)
axis(side = 2, at = seq(0,1,0.25), labels = F, cex.axis = 1.5)
dev.off()


####-------------------------------------------------------------------------------------------
#### Figure 4b: RMSE and r2adj

## Read in the data
rmse.df.w3 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/rmse.df.w3.csv")[-1]
rmse.df.w6 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/rmse.df.w6.csv")[-1]
rmse.df.w9 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/rmse.df.w9.csv")[-1]

jpeg("Fig04b.jpg", width = 5.5, height = 5, units = "in", res = 300)
par(mar = c(3.2,3.2,0.5,0.5))
plot(rmse ~ week, rmse.df.w3, pch = 16, type = "b", col = "#3bb5ff", xlim = c(2.5,9.5), ylim = c(0.02,0.09), axes = F, xlab = "", ylab = ""); par(new = T)
plot(rmse ~ week, rmse.df.w6, pch = 16, type = "b", col = "#be4e08", xlim = c(2.5,9.5), ylim = c(0.02,0.09), axes = F, xlab = "", ylab = ""); par(new = T)
plot(rmse ~ week, rmse.df.w9, pch = 16, type = "b", col = "#80a512", xlim = c(2.5,9.5), ylim = c(0.02,0.09), axes = F, xlab = "", ylab = "")
text(3, rmse.df.w3$rmse[1] - 0.005, label = bquote({italic("r")^2}[italic(adj)] == .(format(rmse.df.w3$r2adj[1], digits = 3))), col = "#3bb5ff",)
text(6, rmse.df.w6$rmse[2] - 0.005, label = bquote({italic("r")^2}[italic(adj)] == .(format(rmse.df.w6$r2adj[2], digits = 3))), col = "#be4e08",)
text(9, rmse.df.w9$rmse[3] - 0.005, label = bquote({italic("r")^2}[italic(adj)] == .(format(rmse.df.w9$r2adj[3], digits = 3))), col = "#80a512")
box()
legend("topright", legend = c("Leaf","Anthesis","Seed"), 
       horiz = T, pch = 16, col = c("#3bb5ff","#be4e08","#80a512"), 
       bty = "n", x.intersp = 0.75, y.intersp = .75, text.width = c(1.09,1,0.4), inset = c(0.2,0))
axis(side = 1, at = seq(3,9,3), labels = c("Leaf","Anthesis","Seed"))
axis(side = 2, seq(0.02,0.1,0.02))
mtext("Developmental stage", side = 1, line = 2.2, cex = 1)
mtext("RMSE", side = 2, line = 2.2, cex = 1)
legend("topleft", legend = "b", text.font = 2, bty = "n", inset = c(-0.0,0))
dev.off()

####-------------------------------------------------------------------------------------------
#### Figure 4c: Absolute abundances

## Read in the data
std.read.counts <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/std.read.counts.csv")

jpeg("Fig04c.jpg", width = 5.5, height = 5, units = "in", res = 300)
par(mar = c(3.2,3.2,0.5,0.5))
plot(std.read.counts$week, std.read.counts$w3_mean, type = "n", lwd = 1.5, xlim = c(2.5,9), ylim = c(4,12), axes = F, ylab = "", xlab = "")
points(std.read.counts$week-0.1, std.read.counts$w3_mean, col = "#3bb5ff", pch = 16)
points(std.read.counts$week, std.read.counts$w6_mean, col = "#be4e08", pch = 16)
points(std.read.counts$week+0.1, std.read.counts$w9_mean, col = "#80a512", pch = 16)
arrows(std.read.counts$week-0.1, std.read.counts$w3_mean, std.read.counts$week-0.1, std.read.counts$w3_upper.ci, col = "#3bb5ff", length = 0.02, angle = 90)
arrows(std.read.counts$week-0.1, std.read.counts$w3_mean, std.read.counts$week-0.1, std.read.counts$w3_lower.ci, col = "#3bb5ff", length = 0.02, angle = 90)
arrows(std.read.counts$week, std.read.counts$w6_mean, std.read.counts$week, std.read.counts$w6_upper.ci, col = "#be4e08", length = 0.02, angle = 90)
arrows(std.read.counts$week, std.read.counts$w6_mean, std.read.counts$week, std.read.counts$w6_lower.ci, col = "#be4e08", length = 0.02, angle = 90)
arrows(std.read.counts$week+0.1, std.read.counts$w9_mean, std.read.counts$week+0.1, std.read.counts$w9_upper.ci, col = "#80a512", length = 0.02, angle = 90)
arrows(std.read.counts$week+0.1, std.read.counts$w9_mean, std.read.counts$week+0.1, std.read.counts$w9_lower.ci, col = "#80a512", length = 0.02, angle = 90)
lines(std.read.counts$week-0.1, std.read.counts$w3_mean, type = "l", lwd = 1.5, col = scales::alpha("#3bb5ff",0.9))
lines(std.read.counts$week, std.read.counts$w6_mean, type = "l", lwd = 1.5, col = scales::alpha("#be4e08",0.9))
lines(std.read.counts$week+0.1, std.read.counts$w9_mean, type = "l", lwd = 1.5, col = scales::alpha("#80a512",0.9))
box()
axis(side = 1, at = seq(3,9,3), labels = c("Leaf","Anthesis","Seed"))
axis(side = 2, at = seq(0,12,2))
mtext("Developmental stage", side = 1, line = 2.2, cex = 1)
mtext("ln(standardized read counts)", side = 2, line = 2.2, cex = 1)
legend("topleft", legend = "c", text.font = 2, bty = "n", inset = c(-0.0,0))
dev.off()

####-------------------------------------------------------------------------------------------
#### Figure 4d: Ecological niches

## Read in the data
can.env <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/canola_environmental.csv", row.names = 1, header = T)

# Use Gower distance - works well for datasets with mixed variables (i.e., plant cover and environmental variables; Legendre and Legendre 2012).
can.env.dist <- vegdist(can.env[,-c(1:2)], "gower")

# Function to evaluate the number of dimensions to use (k)
NMDS.scree <- function(x){
  plot(rep(1,10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress),
       xlim = c(1,10), ylim = c(0,0.30), xlab = "X of Dimensions",
       ylab = "Stress", meain = "NMDS Stress Plot")
  for(i in 1:10){
    points(rep(i+1,10), replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

# Run the scree function
NMDS.scree(can.env.dist)

# 3 looks like a happy balance between low stress and too many dimensions
can.k3 <- metaMDS(can.env.dist, k = 3, autotransform = F)

# Extract the scores for plotting; ranges for plot limits
can.scrs <- vegan::scores(can.k3)
can.xlim <- range(can.scrs[,1])
can.ylim <- range(can.scrs[,2])

jpeg("Fig04d.jpg", width = 5.5, height = 5, units = "in", res = 300)
par(mar = c(3.2,3.2,0.5,0.5))
plot.new()
plot.window(xlim = can.xlim, ylim = can.ylim, asp = 1)
abline(h = 0, lty = "dotted")
abline(v = 0, lty = "dotted")
colvec <- c("#3bb5ff","#be4e08","#80a512") # Set a color for each of my four classes below
vpch <- c(0:2,5,6,15:25) # Plot a different symbol for each line (n = 16)
can.env$line <- as.factor(can.env$line)
can.env$week <- as.factor(can.env$week)
with(can.env, points(can.scrs[,1], can.scrs[,2], pch = 16, 
                     col = colvec[can.env$week]))
ordiellipse(can.k3, groups = can.env$week, label = FALSE, col = colvec, lty = 1, lwd = 1)
axis(side = 1)
axis(side = 2)
mtext("NMDS1", side = 1, line = 2.2, cex = 0.9)
mtext("NMDS2", side = 2, line = 2.2, cex = 0.9)
legend("topleft", legend = "e", text.font = 2, bty = "n", inset = c(0,-0.11))
legend("bottomright", legend = "Stress = 0.057", bty = "n", inset = c(0.06,-0.0))
box()
dev.off()

####-------------------------------------------------------------------------------------------
#### Figure 4e: Prevalence versus abundance

## Read in the data
prev.abu.df <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/pre.abu.df.csv")
prev.abu.df.w3 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/pre.abu.df.w3.csv", row.names = 1)
prev.abu.df.w6 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/pre.abu.df.w6.csv", row.names = 1)
prev.abu.df.w9 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/pre.abu.df.w9.csv", row.names = 1)

# Calculate mean relative abundance
mean(prev.abu.df$abundance)

# Codify generalists (present in >50% of samples) and specialists (present in <50% of samples, >0.1% mean relative abundance)
prev.abu.df$Gen_Spe <- NA
prev.abu.df.w3$Gen_Spe <- NA
prev.abu.df.w6$Gen_Spe <- NA
prev.abu.df.w9$Gen_Spe <- NA
prev.abu.df$Gen_Spe[prev.abu.df$prevalence >= 96] <- "Generalist"
prev.abu.df.w3$Gen_Spe[prev.abu.df.w3$prevalence >= 96] <- "Generalist"
prev.abu.df.w6$Gen_Spe[prev.abu.df.w6$prevalence >= 96] <- "Generalist"
prev.abu.df.w9$Gen_Spe[prev.abu.df.w9$prevalence >= 96] <- "Generalist"
prev.abu.df$Gen_Spe[prev.abu.df$abundance >= 0.1 & prev.abu.df$prevalence < 19] <- "Specialist"
prev.abu.df.w3$Gen_Spe[prev.abu.df.w3$abundance >= 0.1 & prev.abu.df.w3$prevalence < 19] <- "Specialist"
prev.abu.df.w6$Gen_Spe[prev.abu.df.w6$abundance >= 0.1 & prev.abu.df.w6$prevalence < 19] <- "Specialist"
prev.abu.df.w9$Gen_Spe[prev.abu.df.w9$abundance >= 0.1 & prev.abu.df.w9$prevalence < 19] <- "Specialist"

# Make the figure
jpeg("Fig04e.jpg", width = 5.5, height = 5, units = "in", res = 300)
par(mar = c(3.2,3.2,0.5,0.5))
plot(abundance ~ prevalence, prev.abu.df, log = "xy", pch = 16, col = scales::alpha("white",0.15), xlab = "", ylab = "", axes = F)
points(abundance ~ prevalence, prev.abu.df.w3, pch = 16, cex = 1.5, col = "#3bb5ff")
points(abundance ~ prevalence, prev.abu.df.w3[prev.abu.df.w3$Gen_Spe == "Generalist",], pch = 15, cex = 2, col = scales::alpha("#3bb5ff",1))
points(abundance ~ prevalence, prev.abu.df.w6, pch = 16, cex = 1.5, col = "#be4e08")
points(abundance ~ prevalence, prev.abu.df.w6[prev.abu.df.w6$Gen_Spe == "Generalist",], pch = 15, cex = 2, col = scales::alpha("#be4e08",1))
points(abundance ~ prevalence, prev.abu.df.w9, pch = 16, cex = 1.5, col = "#80a512")
points(abundance ~ prevalence, prev.abu.df.w9[prev.abu.df.w9$Gen_Spe == "Generalist",], pch = 15, cex = 2, col = scales::alpha("#80a512",1))
axis(1, labels = F)
axis(2, labels = F)
axis(2)
box()
dev.off()

####-------------------------------------------------------------------------------------------
#### Figure 5: SEM scatterplots

## Read in the data
flowdistbal.df <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/flowdistbal.df.csv", row.names = 2)[-1]

# Run linear models for the significant relationships found in the post-selbal pipeline (P < 0.05)
w3.bal.lm <- lm(performance ~ w3.bal, flowdistbal.df)
w6.bal.lm <- lm(performance ~ w6.bal, flowdistbal.df)
w9.bal.lm <- lm(performance ~ w9.bal, flowdistbal.df)
dist.lm <- lm(performance ~ Distance, flowdistbal.df)

jpeg("Fig05.jpg", width = 7.75, height = 7.5, units = "in", res = 300)
par(mfrow = c(3,3), mar = c(1,1,0,0), oma = c(3,3,1,1))
# Temperature
plot(performance ~ w369.tmp, data = flowdistbal.df, pch = 16, col = scales::alpha("cadetblue3",0.7), ylab = "", xlab = "", axes = F, xlim = c(-3.7,3.2))
legend("topleft", legend = substitute(paste(bold("a"), " Mean temperature")), bty = "n", inset = c(0,0))
box()
axis(side = 1, labels = F)
axis(side = 2)

# Precipitation
plot(performance ~ w369.pre, data = flowdistbal.df, pch = 16, col = scales::alpha("cadetblue3",0.7), ylab = "", xlab = "", axes = F, xlim = c(-3.7,3.2))
legend("topleft", legend = substitute(paste(bold("b"), " Total precipitation")), bty = "n", inset = c(0,0))
box()
axis(side = 1, labels = F)
axis(side = 2, labels = F)

# Genetic distance
plot(performance ~ Distance, data = flowdistbal.df, pch = 16, col = scales::alpha("#84a955",0.7), ylab = "", xlab = "", axes = F, xlim = c(-3.7,3.2))
legend("topleft", legend = substitute(paste(bold("c"), " Genetic distance")), bty = "n", inset = c(0,0))
box()
axis(side = 1, labels = F)
axis(side = 2, labels = F)

# Flowering differential w3
plot(performance ~ w3.start.diff, data = flowdistbal.df, pch = 16, col = scales::alpha("#84a955",0.7), ylab = "", xlab = "", axes = F, xlim = c(-3.7,3.2))
legend("topleft", legend = substitute(paste(bold("d"), " Leaf")), bty = "n", inset = c(0,0))
legend("topleft", legend = "development flowering", bty = "n", inset = c(0.038,0.045))
legend("topleft", legend = "time differential", bty = "n", inset = c(0.038,0.09))
box()
axis(side = 1, labels = F)
axis(side = 2)

# Flowering differential w6
plot(performance ~ w6.start.diff, data = flowdistbal.df, pch = 16, col = scales::alpha("#84a955",0.7), ylab = "", xlab = "", axes = F, xlim = c(-3.7,3.2))
legend("topleft", legend = substitute(paste(bold("e"), " Anthesis")), bty = "n", inset = c(0,0))
legend("topleft", legend = "flowering time", bty = "n", inset = c(0.038,0.045))
legend("topleft", legend = "differential", bty = "n", inset = c(0.038,0.09))
box()
axis(side = 1, labels = F)
axis(side = 2, labels = F)

# Flowering differential w9
plot(performance ~ w9.start.diff, data = flowdistbal.df, pch = 16, col = scales::alpha("#84a955",0.7), ylab = "", xlab = "", axes = F, xlim = c(-3.7,3.2))
legend("topleft", legend = substitute(paste(bold("f"), " Seed")), bty = "n", inset = c(0,0))
legend("topleft", legend = "development flowering", bty = "n", inset = c(0.028,0.045))
legend("topleft", legend = "time differential", bty = "n", inset = c(0.028,0.09))
box()
axis(side = 1, labels = F)
axis(side = 2, labels = F)

# Week 3 balance
plot(performance ~ w3.bal, data = flowdistbal.df, pch = 16, col = scales::alpha("#965da7",0.7), ylab = "", xlab = "", axes = F, xlim = c(-3.7,3.2))
abline(coef(w3.bal.lm), lty = 3)
legend("topleft", legend = substitute(paste(bold("g"), " Leaf development")), bty = "n", inset = c(0,0))
legend("topleft", legend = "global balance", bty = "n", inset = c(0.038,0.045))
box()
axis(side = 1, labels = T)
axis(side = 2)

# Week 6 balance
plot(performance ~ w6.bal, data = flowdistbal.df, pch = 16, col = scales::alpha("#965da7",0.7), ylab = "", xlab = "", axes = F, xlim = c(-3.7,3.2))
abline(coef(w6.bal.lm), lty = 3)
legend("topleft", legend = substitute(paste(bold("h"), " Anthesis global balance")), bty = "n", inset = c(0,0))
box()
axis(side = 1, labels = T)
axis(side = 2, labels = F)

# Week 9 balance
plot(performance ~ w9.bal, data = flowdistbal.df, pch = 16, col = scales::alpha("#965da7",0.7), ylab = "", xlab = "", axes = F, xlim = c(-3.7,3.2))
abline(coef(w9.bal.lm), lty = 3)
legend("topleft", legend = substitute(paste(bold("i"), " Seed development")), bty = "n", inset = c(0,0))
legend("topleft", legend = "global balance", bty = "n", inset = c(0.026,0.045))
box()
axis(side = 1, labels = T)
axis(side = 2, labels = F)

mtext("Standardized explanatory variable", side = 1, outer = T, line = 1.25, cex = 0.75)
mtext("Yield-weighted performance", side = 2, outer = T, line = 1.25, cex = 0.75)
dev.off()

####-------------------------------------------------------------------------------------------
#### Figure 6a: Phylogenetic tree for 63 ASVs

# Read in the selbal bal.tab output using the weekly datasets, and taxonomy
Bal.Tab.full.w3.imp1 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/Bal.Tab.full.imp1.w3.csv", row.names = 1)
Bal.Tab.full.w6.imp1 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/Bal.Tab.full.imp1.w6.csv", row.names = 1)
Bal.Tab.full.w9.imp1 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/Bal.Tab.full.imp1.w9.csv", row.names = 1)
names(Bal.Tab.full.w3.imp1) <- c("taxa","group")
names(Bal.Tab.full.w6.imp1) <- c("taxa","group")
names(Bal.Tab.full.w9.imp1) <- c("taxa","group")
LMS_GB_full_imp1_w3 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/LMS_2016_17_GlobalBalance_full_imp1.w3.csv", row.names = 1)
LMS_GB_full_imp1_w6 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/LMS_2016_17_GlobalBalance_full_imp1.w6.csv", row.names = 1)
LMS_GB_full_imp1_w9 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/LMS_2016_17_GlobalBalance_full_imp1.w9.csv", row.names = 1)
can.tax <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/tax.root.reps_dupes_merged.prev_and_abun_filt.csv", row.names = 1)
can.asv <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/asv.root.reps_dupes_merged.prev_and_abun_filt.csv", row.names = 1)

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

# Read in the tree and make sure ASV identifiers starting with a numeric have an "X" prepended
# root.tree <- phytools::read.newick("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/insertion-tree.nwk")
root.tree <- phyloseq::read_tree_greengenes("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/insertion-tree.nwk")
x <- root.tree$tip.label
x[grepl("^[[:digit:]]", x)] <- paste("X", x[grepl("^[[:digit:]]", x)], sep = "")
root.tree$tip.label <- x

# Subset the tree to include only the ASVs of interest (this is to avoid errors due to the large tree size and computer memory)
# This is the error you would receive if the tree is too large:
# Error in double(nm * nm) : vector size cannot be NA 
# In addition: Warning message: In nm * nm : NAs produced by integer overflow

# Weeks 3, 6, and 9 tree
w369.tree <- ape::drop.tip(root.tree, root.tree$tip.label[-match(unique(c(Bal.Tab.full.w3.imp1$taxa, Bal.Tab.full.w6.imp1$taxa, Bal.Tab.full.w9.imp1$taxa)), root.tree$tip.label)])

# Make a df that has each of the 63 ASVs and their selbal group (numerator, denominator)
rownames(Bal.Tab.full.w3.imp1) <- Bal.Tab.full.w3.imp1$taxa
rownames(Bal.Tab.full.w6.imp1) <- Bal.Tab.full.w6.imp1$taxa
rownames(Bal.Tab.full.w9.imp1) <- Bal.Tab.full.w9.imp1$taxa
w369.bal.tab <- merge(Bal.Tab.full.w3.imp1, Bal.Tab.full.w6.imp1, by = 0, all = T)
rownames(w369.bal.tab) <- w369.bal.tab$Row.names
w369.bal.tab$Row.names <- NULL
w369.bal.tab <- merge(w369.bal.tab, Bal.Tab.full.w9.imp1, by = 0, all = T)
rownames(w369.bal.tab) <- w369.bal.tab$Row.names
w369.bal.tab$Row.names <- w369.bal.tab$Taxa.x <- w369.bal.tab$Taxa.y <- w369.bal.tab$Taxa <- NULL
w369.bal.tab <- w369.bal.tab[,c(2,4,6)]
names(w369.bal.tab) <- c("w3.group","w6.group","w9.group")

# Extract the taxonomy to a file for reference purposes
w369.tax <- can.tax[w369.tree$tip.label,]
write.csv(w369.tax, "~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/taxonomy63.csv")

# Subset relative abundance to 63 ASVs from the three communities,
# and create a df that contains all the relevant plotting information for the tree (relative abundances, color, etc)
w369.rel0 <- can.asv.rel0[,unique(c(Bal.Tab.full.w3.imp1$taxa, Bal.Tab.full.w6.imp1$taxa, Bal.Tab.full.w9.imp1$taxa))]
w3.relabu <- w369.rel0[rownames(LMS_GB_full_imp1_w3),]
w6.relabu <- w369.rel0[rownames(LMS_GB_full_imp1_w6),]
w9.relabu <- w369.rel0[rownames(LMS_GB_full_imp1_w9),]
w3.relabu.means <- as.data.frame(colMeans(w3.relabu))
w6.relabu.means <- as.data.frame(colMeans(w6.relabu))
w9.relabu.means <- as.data.frame(colMeans(w9.relabu))
identical(rownames(w3.relabu.means), rownames(w6.relabu.means))
identical(rownames(w3.relabu.means), rownames(w9.relabu.means))
w369.relabu.means <- cbind.data.frame(w3.relabu.means, w6.relabu.means, w9.relabu.means)
names(w369.relabu.means) <- c("w3.comm","w6.comm","w9.comm")
w369.tax.rel.merge <- merge(w369.tax, w369.relabu.means, by = 0)
rownames(w369.tax.rel.merge) <- w369.tax.rel.merge$Row.names
w369.tax.rel.bal.merge <- merge(w369.tax.rel.merge, w369.bal.tab, by = 0)
rownames(w369.tax.rel.bal.merge) <- w369.tax.rel.bal.merge$Row.names
w369.tax.rel.bal.merge$Row.names <- w369.tax.rel.bal.merge$Row.names <- NULL
w369.tax.rel.bal.merge <- w369.tax.rel.bal.merge[w369.tree$tip.label,]

# Make sure the df matches the tree. Add in plotting information to be used in ggtree.
identical(rownames(w369.tax.rel.bal.merge), w369.tree$tip.label)
w369.tax.rel.bal.merge$w3.group <- as.character(w369.tax.rel.bal.merge$w3.group)
w369.tax.rel.bal.merge$w6.group <- as.character(w369.tax.rel.bal.merge$w6.group)
w369.tax.rel.bal.merge$w9.group <- as.character(w369.tax.rel.bal.merge$w9.group)
w369.tax.rel.bal.merge$w3.abu.den <- 0
w369.tax.rel.bal.merge$w6.abu.den <- 0
w369.tax.rel.bal.merge$w9.abu.den <- 0
w369.tax.rel.bal.merge$w3.abu.num <- 0
w369.tax.rel.bal.merge$w6.abu.num <- 0
w369.tax.rel.bal.merge$w9.abu.num <- 0
w369.tax.rel.bal.merge$w3.abu.den[which(w369.tax.rel.bal.merge$w3.group == "DEN")] <- w369.tax.rel.bal.merge$w3.comm[which(w369.tax.rel.bal.merge$w3.group == "DEN")]
w369.tax.rel.bal.merge$w3.abu.num[which(w369.tax.rel.bal.merge$w3.group == "NUM")] <- w369.tax.rel.bal.merge$w3.comm[which(w369.tax.rel.bal.merge$w3.group == "NUM")]
w369.tax.rel.bal.merge$w6.abu.den[which(w369.tax.rel.bal.merge$w6.group == "DEN")] <- w369.tax.rel.bal.merge$w6.comm[which(w369.tax.rel.bal.merge$w6.group == "DEN")]
w369.tax.rel.bal.merge$w6.abu.num[which(w369.tax.rel.bal.merge$w6.group == "NUM")] <- w369.tax.rel.bal.merge$w6.comm[which(w369.tax.rel.bal.merge$w6.group == "NUM")]
w369.tax.rel.bal.merge$w9.abu.den[which(w369.tax.rel.bal.merge$w9.group == "DEN")] <- w369.tax.rel.bal.merge$w9.comm[which(w369.tax.rel.bal.merge$w9.group == "DEN")]
w369.tax.rel.bal.merge$w9.abu.num[which(w369.tax.rel.bal.merge$w9.group == "NUM")] <- w369.tax.rel.bal.merge$w9.comm[which(w369.tax.rel.bal.merge$w9.group == "NUM")]
w369.tax.rel.bal.merge <- droplevels(w369.tax.rel.bal.merge)
w369.tax.rel.bal.merge$w3.col <- as.character(NA)
w369.tax.rel.bal.merge$w6.col <- as.character(NA)
w369.tax.rel.bal.merge$w9.col <- as.character(NA)
w369.tax.rel.bal.merge$w3.col <- "#b94663"
w369.tax.rel.bal.merge$w6.col <- "#bc7d39"
w369.tax.rel.bal.merge$w9.col <- "#6fab5a"
w369.tax.rel.bal.merge$w3.pch <- as.numeric(NA)
w369.tax.rel.bal.merge$w6.pch <- as.numeric(NA)
w369.tax.rel.bal.merge$w9.pch <- as.numeric(NA)
w369.tax.rel.bal.merge$w3.pch[w369.tax.rel.bal.merge$w3.group == "DEN"] <- 16
w369.tax.rel.bal.merge$w3.pch[w369.tax.rel.bal.merge$w3.group == "NUM"] <- 17
w369.tax.rel.bal.merge$w6.pch[w369.tax.rel.bal.merge$w6.group == "DEN"] <- 16
w369.tax.rel.bal.merge$w6.pch[w369.tax.rel.bal.merge$w6.group == "NUM"] <- 17
w369.tax.rel.bal.merge$w9.pch[w369.tax.rel.bal.merge$w9.group == "DEN"] <- 16
w369.tax.rel.bal.merge$w9.pch[w369.tax.rel.bal.merge$w9.group == "NUM"] <- 17
w369.tax.rel.bal.merge$w3.col[w369.tax.rel.bal.merge$w3.group == "DEN"] <- "#b94663"
w369.tax.rel.bal.merge$w3.col[w369.tax.rel.bal.merge$w3.group == "NUM"] <- "#5794d7"
w369.tax.rel.bal.merge$w6.col[w369.tax.rel.bal.merge$w6.group == "DEN"] <- "#b94663"
w369.tax.rel.bal.merge$w6.col[w369.tax.rel.bal.merge$w6.group == "NUM"] <- "#5794d7"
w369.tax.rel.bal.merge$w9.col[w369.tax.rel.bal.merge$w9.group == "DEN"] <- "#b94663"
w369.tax.rel.bal.merge$w9.col[w369.tax.rel.bal.merge$w9.group == "NUM"] <- "#5794d7"
w369.tax.rel.bal.merge$w3.col[is.na(w369.tax.rel.bal.merge$w3.col)] <- "white"
w369.tax.rel.bal.merge$w6.col[is.na(w369.tax.rel.bal.merge$w6.col)] <- "white"
w369.tax.rel.bal.merge$w9.col[is.na(w369.tax.rel.bal.merge$w9.col)] <- "white"
w369.tax.rel.bal.merge$genus <- as.character(w369.tax.rel.bal.merge$genus)
w369.tax.rel.bal.merge$genus2 <- make.unique(w369.tax.rel.bal.merge$genus)  # Make genera unique (appends numbers to matches)
identical(rownames(w369.tax.rel.bal.merge), w369.tree$tip.label)

# Numerator and denominator abundance
w369.tax.rel.bal.merge$w3.abundance <- c(w369.tax.rel.bal.merge$w3.abu.num + w369.tax.rel.bal.merge$w3.abu.den)
w369.tax.rel.bal.merge$w6.abundance <- c(w369.tax.rel.bal.merge$w6.abu.num + w369.tax.rel.bal.merge$w6.abu.den)
w369.tax.rel.bal.merge$w9.abundance <- c(w369.tax.rel.bal.merge$w9.abu.num + w369.tax.rel.bal.merge$w9.abu.den)

# Now can finally make the tree
merge2a <- w369.tax.rel.bal.merge[w369.tree$tip.label,]     # Order the merged file by the tip labels (i.e., species)
identical(rownames(merge2a), w369.tree$tip.label)           # Always double check they match
w369.tree$tip.label <- merge2a$genus2                       # Add unique genera as the tip labels

View(merge2a[c(rownames(prev.abu.df.w3[!is.na(prev.abu.df.w3$Gen_Spe),]),
          rownames(prev.abu.df.w6[!is.na(prev.abu.df.w6$Gen_Spe),]),
          rownames(prev.abu.df.w9[!is.na(prev.abu.df.w9$Gen_Spe),])),])

## Make the tree and label the nodes sequentially to determine coloring, etc. below
p <- ggtree(w369.tree, layout = "circular", size = 0.2)
p + geom_label(aes(x=branch, label=node)) + 
  geom_tiplab(aes(angle = angle, cex = 1, label = gsub(" ","_",label)), offset = .15, parse = T, show.legend=F)

nlevels(as.factor(merge2a$phylum)) # n = 5
nlevels(as.factor(merge2a$genus)) # n = 38

# I inspected the tree from above and manually assigned the colors based on phylum
acidobacteria <- c(52,53,114)
actinobacteria <- c(54:62,115:123)
bacteroidetes <- c(43:51,105:112)
chloroflexi <- 63
proteobacteria <- c(1:42,64:103)
black.col <- c(104,113)

# This is just a check that I've got everything (the row number should match with the assigned number)
sort(c(acidobacteria,
       actinobacteria,
       bacteroidetes,
       chloroflexi,
       proteobacteria,
       black.col))

# Prepare the phylum color for the tree
d <- data.frame(node = 1:123, color = "black")
d$color[acidobacteria] <- "#f806a0"
d$color[actinobacteria] <- "#3eb73c"
d$color[bacteroidetes] <- "#a53ec7"
d$color[chloroflexi] <- "#3f62fb"
d$color[proteobacteria] <- "#ca7b0a"

# Now to make a new tree that includes the color information
p <- ggtree(w369.tree,layout = "circular", size = 0.2) %<+% d + aes(color=I(color)) + xlim(NA, 2.05)

# Add relative abundances
p$data$leaf <- NA
p$data$leaf[c(1:63)] <- merge2a$w3.comm
p$data$leaf[which(p$data$leaf == 0)] <- NA
p$data$anthesis <- NA
p$data$anthesis[c(1:63)] <- merge2a$w6.comm
p$data$anthesis[which(p$data$anthesis == 0)] <- NA
p$data$seed <- NA
p$data$seed[c(1:63)] <- merge2a$w9.comm
p$data$seed[which(p$data$seed == 0)] <- NA
p$data$phylum <- NA
p$data$phylum[c(1:63)] <- as.character(gsub(" ", "_", merge2a$phylum))
p$data$genus <- NA
p$data$genus[c(1:63)] <- merge2a$genus
p$data$genus[c(1:63)] <- as.character(gsub("-", "_", p$data$genus[c(1:63)]))
p$data$genus[c(1:63)] <- as.character(gsub(" ", "_", p$data$genus[c(1:63)]))
p$data$genus[63] <- "X1959_1"
p$data$leaf.col <- "gray"
p$data$anthesis.col <- "gray"
p$data$seed.col <- "gray"
w3.labels.num <- merge2a$genus2[merge2a$w3.group == "NUM" & is.na(merge2a$w3.group) == FALSE]
w3.labels.den <- merge2a$genus2[merge2a$w3.group == "DEN" & is.na(merge2a$w3.group) == FALSE]
p$data$leaf.col[p$data$label %in% w3.labels.num] <- "#5794d7"
p$data$leaf.col[p$data$label %in% w3.labels.den] <- "#b94663"
w6.labels.num <- merge2a$genus2[merge2a$w6.group == "NUM" & is.na(merge2a$w6.group) == FALSE]
w6.labels.den <- merge2a$genus2[merge2a$w6.group == "DEN" & is.na(merge2a$w6.group) == FALSE]
p$data$anthesis.col[p$data$label %in% w6.labels.num] <- "#5794d7"
p$data$anthesis.col[p$data$label %in% w6.labels.den] <- "#b94663"
w9.labels.num <- merge2a$genus2[merge2a$w9.group == "NUM" & is.na(merge2a$w9.group) == FALSE]
w9.labels.den <- merge2a$genus2[merge2a$w9.group == "DEN" & is.na(merge2a$w9.group) == FALSE]
p$data$seed.col[p$data$label %in% w9.labels.num] <- "#5794d7"
p$data$seed.col[p$data$label %in% w9.labels.den] <- "#b94663"

## I import this tree to Photoshop to beautify it
jpeg("tree1.jpg", width = 8, height = 8, units = "in", res = 600)
p2 <- p + geom_tiplab(aes(x = rep(0.5,123), angle = angle, cex = 1, label = paste0('italic(', genus,')')), size = 2, offset = 0.7, parse = T, show.legend = F) +
  theme_tree(plot.margin = margin(0, 120, 0, 60)) +
  geom_tippoint(aes(x = rep(0.76,123), size = leaf, color = leaf.col), alpha = 0.65, stroke = 0, show.legend = T) +
  geom_tippoint(aes(x = rep(0.91,123), size = anthesis, color = anthesis.col), alpha = 0.65, stroke = 0, show.legend = F) +
  geom_tippoint(aes(x = rep(1.06,123), size = seed, color = seed.col), alpha = 0.65, stroke = 0, show.legend = F) +
  scale_size_continuous(range = c(0.5, 5.5)) +
  geom_text(aes(x = 0, y = 0, label = ""), show.legend=F) +
  theme(legend.position = c(1,0.5)) +
  labs(size = "Relative abundance (%)", colour = "Phylum")
p2
dev.off()

##-------------------------------------------------------------------------------------------
## Figure 6b: Phylogenetic dissimilarity (abundance-weighted beta mean nearest taxon distance; βMNTD)

beta.mntd.weighted.63 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/beta.mntd.weighted.63.csv", row.names = 1)

# NMDS
monoNMDS.scree <- function(x){
  plot(rep(1,10), replicate(10, monoMDS(x, autotransform = F, k = 1)$stress),
       xlim = c(1,10), ylim = c(0,0.30), xlab = "X of Dimensions",
       ylab = "Stress", main = "NMDS Stress Plot")
  for(i in 1:10){
    points(rep(i+1,10), replicate(10, monoMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

monoNMDS.scree(beta.mntd.weighted.63)

set.seed(42)
w369.nmds <- monoMDS(beta.mntd.weighted.63, k = 4, autotransform = F) # stress = 0.195

# Extract the scores for plotting; ranges for plot limits
w369.scrs <- vegan::scores(w369.nmds)
w369.xlim <- range(w369.scrs[,1])
w369.ylim <- range(w369.scrs[,2])
w369.meta <- data.frame(comm = rep(c("w3","w6","w9"), each = 64))

jpeg("Fig06B.jpg", width = 5.5, height = 5, units = "in", res = 300)
par(mar = c(3.2,3.2,0.5,0.5))
plot.new()
plot.window(xlim = w369.xlim, ylim = w369.ylim, asp = 1)
par(mar = c(3.2, 3.2, 0.5, 0.5))
points(w369.scrs[c(1:64),c(1,2)], pch = 16, col = scales::alpha("#3bb5ff",0.6), cex = 1)
points(w369.scrs[c(65:128),c(1,2)], pch = 16, col = scales::alpha("#be4e08",0.6), cex = 1)
points(w369.scrs[c(129:192),c(1,2)], pch = 16, col = scales::alpha("#80a512",0.6), cex = 1)
ordiellipse(w369.nmds, groups = w369.meta$comm, label = FALSE, col = c("#3bb5ff","#be4e08","#80a512"), lty = 1, lwd = 1.5)
legend("topleft", legend = c("L","F","S"), col = c("#3bb5ff","#be4e08","#80a512"), pch = 16, bty = "n")
axis(side = 1); axis(side = 2)
box()
mtext("NMDS1", side = 1, line = 2.2, cex = 1)
mtext("NMDS2", side = 2, line = 2.2, cex = 1)
dev.off()

####-------------------------------------------------------------------------------------------
#### Figure S2: Rarefaction

load("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/FigureS2.RData")

## Make the figure locally
jpeg("FigS2.jpg", width = 5, height = 4, units = "in", res = 300)
par(mar = c(3.3,3.2,0.1,0.9))
plot.new()
plot.window(xlim = c(15,10000), ylim = c(25,1100))
for (i in seq_along(rare.df0)) {
  # for (i in c(1:100)) {
  N <- attr(rare.df0[[i]], "Subsample")
  lines(N, rare.df0[[i]], col = scales::alpha("gray10",0.5))
}
box()
axis(1); axis(2)
mtext("Sample size", side = 1, line = 2.2)
mtext("Richness", side = 2, line = 2.2)
abline(v = seq(0,10000,2000)-20, lty = 3, lwd = 0.75, col = "gray85")
abline(h = seq(0,1000,200)-20, lty = 3, lwd = 0.75, col = "gray85")
dev.off()

####-------------------------------------------------------------------------------------------
#### Figure S3: Yield performance metrics

jpeg("FigS3.jpg", width = 12, height = 4.0, units = "in", res = 300)
par(mfrow = c(1,3), mar = c(1,1,0.1,0.1), oma = c(2.75,2.75,0.1,0.1))
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
abline(h = 0.5, v = 0.5, lty = 3, col = "gray25")
points(can.sample7$yield.std, can.sample7$s.std, pch = 16, cex = can.sample7$importance1*10, col = scales::alpha("blue", 0.6), )
box()
axis(side = 1, at = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1"), cex.axis = 1.5)
axis(side = 2, at = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1"), cex.axis = 1.5)
legend("bottomleft", legend = substitute(paste(bold("a"), " Yield-weighted")), bty = "n", cex = 1.5, inset = c(0,0))

plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
abline(h = 0.5, v = 0.5, lty = 3, col = "gray25")
points(can.sample7$yield.std, can.sample7$s.std, pch = 16, cex = can.sample7$importance2*10, col = scales::alpha("blue", 0.6), )
box()
axis(side = 1, at = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1"), cex.axis = 1.5)
axis(side = 2, at = seq(0,1,0.25), labels = F)
legend("bottomleft", legend = substitute(paste(bold("b"), " Stability-weighted")), bty = "n", cex = 1.5, inset = c(0,0))

plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
abline(h = 0.5, v = 0.5, lty = 3, col = "gray25")
points(can.sample7$yield.std, can.sample7$s.std, pch = 16, cex = can.sample7$importance3*10, col = scales::alpha("blue", 0.6), )
box()
axis(side = 1, at = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1"), cex.axis = 1.5)
axis(side = 2, at = seq(0,1,0.25), labels = F)
legend("bottomleft", legend = substitute(paste(bold("c"), " Unweighted")), bty = "n", cex = 1.5, inset = c(0,0))

mtext("Yield (scaled)", side = 1, line = 1.45, outer = T)
mtext("Yield stability (scaled)", side = 2, line = 1.45, outer = T)
dev.off()

####-------------------------------------------------------------------------------------------
#### Figure S4

# Read in the results of the parallelized selbal pipeline for each performance metric
bal.imp1 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/LMS_2016_17_GlobalBalance_full_imp1.2.csv", header = T)
bal.imp2 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/LMS_2016_17_GlobalBalance_full_imp2.2.csv", header = T)
bal.imp3 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/LMS_2016_17_GlobalBalance_full_imp3.2.csv", header = T)

# Output r2.adj for each yield performance metric
imp1.r2 <- round(RsquareAdj(lm(numy ~ V1, bal.imp1))$adj.r.squared,3)
imp2.r2 <- round(RsquareAdj(lm(numy ~ V1, bal.imp2))$adj.r.squared,3)
imp3.r2 <- round(RsquareAdj(lm(numy ~ V1, bal.imp3))$adj.r.squared,3)

jpeg("FigS4.jpg", width = 3, height = 7.5, units = "in", res = 300)
par(mfrow = c(3,1), mar = c(1.5,2,1,1), oma = c(2,1.5,0.1,0.1))

plot(numy ~ V1, bal.imp1, pch = 16, col = scales::alpha("black",0.6), axes = F, xlab = "", ylab = "", ylim = c(0.02878908, 0.38530741))
abline(coef(lm(numy ~ V1, bal.imp1)), lty = 2, lwd = 1.5)
axis(1); axis(2); box()
legend("topleft", legend = substitute(paste(bold("a"), " Yield-weighted")), bty = "n", inset = c(-0.01,0.05))
legend("bottomright", legend = bquote({italic("r")^2}[italic(adj)] == .(format(imp1.r2, digits = 3))), bty = "n", inset = c(0,0))

plot(numy ~ V1, bal.imp2, pch = 16, col = scales::alpha("black",0.6), axes = F, xlab = "", ylab = "", ylim = c(0.02878908, 0.38530741))
abline(coef(lm(numy ~ V1, bal.imp2)), lty = 2, lwd = 1.5)
axis(1); axis(2); box()
mtext("Yield performance", side = 2, line = 2.2, cex = 0.7)
legend("topleft", legend = substitute(paste(bold("b"), " Yield stability-weighted")), bty = "n", inset = c(-0.01,0.05))
legend("bottomright", legend = bquote({italic("r")^2}[italic(adj)] == .(format(imp2.r2, digits = 3))), bty = "n", inset = c(0,0))

plot(numy ~ V1, bal.imp3, pch = 16, col = scales::alpha("black",0.6), axes = F, xlab = "", ylab = "", ylim = c(0.02878908, 0.38530741))
abline(coef(lm(numy ~ V1, bal.imp3)), lty = 2, lwd = 1.5)
axis(1); axis(2); box()
mtext("Global balance", side = 1, line = 2.2, cex = 0.7)
legend("topleft", legend = substitute(paste(bold("c"), " Unweighted")), bty = "n", inset = c(-0.01,0.05))
legend("bottomright", legend = bquote({italic("r")^2}[italic(adj)] == .(format(imp3.r2, digits = 3))), bty = "n", inset = c(0,0))
dev.off()

####-------------------------------------------------------------------------------------------
#### Figure S5: B. napus genetic distances and global balances for each line

## Create the balance df
split.df <- unlist(strsplit(rownames(LMS_GB_full_imp1_w3),"[.]"))
line.2 <- split.df[seq(3,length(split.df),5)]
balance.df2 <- cbind.data.frame(CanolaLine = rep(line.2,3),
                                rbind.data.frame(LMS_GB_full_imp1_w3,
                                                 LMS_GB_full_imp1_w6,
                                                 LMS_GB_full_imp1_w9))
balance.df.w3 <- balance.df2[c(1:64),]
balance.df.w6 <- balance.df2[c(65:128),]
balance.df.w9 <- balance.df2[c(129:192),]
balance.df.w3$CanolaLine <- as.factor(balance.df.w3$CanolaLine)
balance.df.w6$CanolaLine <- as.factor(balance.df.w6$CanolaLine)
balance.df.w9$CanolaLine <- as.factor(balance.df.w9$CanolaLine)

# Re-arrange the lines
balance.df.w3$CanolaLine = factor(balance.df.w3$CanolaLine,
                                  levels(balance.df.w3$CanolaLine)[c(15,12,11,3,8,6,13,9,1,2,5,16,10,14,7,4)])   # Reorder the factors so that NAM-5 is after NAM-0
balance.df.w6$CanolaLine = factor(balance.df.w6$CanolaLine,
                                  levels(balance.df.w6$CanolaLine)[c(15,12,11,3,8,6,13,9,1,2,5,16,10,14,7,4)])   # Reorder the factors so that NAM-5 is after NAM-0
balance.df.w9$CanolaLine = factor(balance.df.w9$CanolaLine,
                                  levels(balance.df.w9$CanolaLine)[c(15,12,11,3,8,6,13,9,1,2,5,16,10,14,7,4)])   # Reorder the factors so that NAM-5 is after NAM-0

# Make a color vector and add it to the df
can_line_cols <- c("#8d2759","#543f96","#6e76dd","#c8533b",
                   "#55a959","#c3ae3b","#cf7bce","#45c097",
                   "#b9475f","#ba464d","#ca8c43","#d86395",
                   "#6d8dd7","#82317c","#86a042","#8e5628")
balance.df.w3$colour <- NA                                          # Create a colour column
balance.df.w3$colour <- can_line_cols[balance.df.w3$CanolaLine]     # Assign colours
balance.df.w6$colour <- NA                                          # Create a colour column
balance.df.w6$colour <- can_line_cols[balance.df.w6$CanolaLine]     # Assign colours
balance.df.w9$colour <- NA                                          # Create a colour column
balance.df.w9$colour <- can_line_cols[balance.df.w9$CanolaLine]     # Assign colours
bymedian.w3 <- with(balance.df.w3, reorder(CanolaLine, V1, median, na.rm = T))
bymedian.w6 <- with(balance.df.w6, reorder(CanolaLine, V1, median, na.rm = T))
bymedian.w9 <- with(balance.df.w9, reorder(CanolaLine, V1, median, na.rm = T))
xlabs <- levels(bymedian.w9)
xlabs[xlabs == "NAM-48"] <- "DH27298"

## Read in the genetic distance data
nam.dist <- read.table("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/Select_NAM_genetic_distances.txt", 
                         header = T, stringsAsFactors = F)

# Remove 25576-2 and 25576-3, use 25576-1 as YN04-C1213
na.dist.sub <- subset(nam.dist, row != "25576-2" & col != "25576-3")
na.dist.sub <- subset(na.dist.sub, row != "25576-3" & col != "25576-2")

# Rename "25576-1" to "YN04-C1213", and replace underscores with dashes
na.dist.sub[na.dist.sub$row == "25576-1",1] <- "YN04-C1213"
na.dist.sub[na.dist.sub$col == "25576-1",2] <- "YN04-C1213"
na.dist.sub$row <- gsub("_", "-", na.dist.sub$row)
na.dist.sub$col <- gsub("_", "-", na.dist.sub$col)

# Recast the subset data into a wide format matrix for plotting
nam.dist.wide <- spread(na.dist.sub, col, value)
nam.dist.wide.mat <- as.matrix(nam.dist.wide[-1])
rownames(nam.dist.wide.mat) <- nam.dist.wide$row

# Hierarchical clustering based on distance
nam.dist.euc <- vegdist(nam.dist.wide.mat, "euclidean")
nam.hc <- hclust(nam.dist.euc, "average")
chwo.spe <- reorder(nam.hc, as.dist(nam.dist.wide.mat))

# Assign the clusters to the df
tree.labels <- c("NAM-46","NAM-79","NAM-72","NAM-17","NAM-32","NAM-0","NAM-30","YN04-C1213",
                 "DH27298","NAM-23","NAM-13","NAM-76","NAM-5","NAM-37","NAM-14","NAM-43")

# Make the dendrogram
hcd <- as.dendrogram(chwo.spe)

# Make the figure. Panel a will need some post-hoc tidying of tip labels
jpeg("FigS5.jpg", width = 3, height = 7.5, units = "in", res = 300)
par(mfrow = c(4,1), mar = c(3.5,1,0.2,1), oma = c(5.5,1,0,1))

# Genetic distance
plot(hcd,  xlab = "", horiz = TRUE)
mtext(side = 1, "Genetic distance (cM)", line = 2.2, cex = 0.65)
legend("topleft", legend = substitute(paste(bold("a"), "")), bty = "n", inset = c(0.02,-0.05))

# Week 3 balance
par(mar = c(1,3,0.5,1))
bp <- boxplot(V1 ~ CanolaLine, data = balance.df.w3, col = can_line_cols, xlab = "", ylab = "", axes = F)
stripchart(V1 ~ bymedian.w3, data = balance.df.w3, vertical = TRUE, method = "jitter", 
           add = TRUE, bg = "gray50", pch = 21, cex = 0.75)
axis(side = 1, at = c(1:16), labels = F)
axis(side = 2, seq(-4,12,4))
box()
legend("topleft", legend = substitute(paste(bold("b"), " Leafing")), bty = "n", inset = c(0.02,0.01))

# Week 6 balance
par(mar = c(1,3,0.5,1))
bp <- boxplot(V1 ~ bymedian.w6, data = balance.df.w6, col = can_line_cols, xlab = "", ylab = "", axes = F)
stripchart(V1 ~ bymedian.w6, data = balance.df.w6, vertical = TRUE, method = "jitter", 
           add = TRUE, bg = "gray50", pch = 21, cex = 0.75)
axis(side = 1, at = c(1:16), labels = F)
axis(side = 2, seq(-4,12,4))
box()
legend("topleft", legend = substitute(paste(bold("c"), " Flowering")), bty = "n", inset = c(0.02,0.01))
mtext("Global balance", side = 2, line = 2.2, cex = 0.65)

# Week 6 balance
par(mar = c(1,3,0.5,1))
bp <- boxplot(V1 ~ bymedian.w9, data = balance.df.w9, col = can_line_cols, xlab = "", ylab = "", axes = F)
stripchart(V1 ~ bymedian.w9, data = balance.df.w9, vertical = TRUE, method = "jitter", 
           add = TRUE, bg = "gray50", pch = 21, cex = 0.75)
axis(side = 1, at = c(1:16), labels = xlabs, las = 2)
axis(side = 2, seq(-4,12,4))
box()
legend("topleft", legend = substitute(paste(bold("c"), " Seeding")), bty = "n", inset = c(0.02,0.01))
dev.off()

####-------------------------------------------------------------------------------------------
#### Figure S7: Stacked barplots

## Read in the data
can.asvs <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/asv.root.reps_dupes_merged.prev_and_abun_filt.csv", row.names = 1)
tax.63 <- read.csv("~/Dropbox/CFREF Work/SteveM/Canola winnowing/Rcode/Curated_Code/Data/taxonomy63.csv", row.names = 1)

# Make the sample df
Sample.Code <- rownames(can.asvs)
split.df <- unlist(strsplit(Sample.Code,"[.]"))
Site <- split.df[seq(1,length(split.df),5)]
Year <- split.df[seq(2,length(split.df),5)]
CanolaLine <- split.df[seq(3,length(split.df),5)]
Week <- split.df[seq(4,length(split.df),5)]
can.sam <- cbind.data.frame(SampleID = paste(Site, Year, CanolaLine, Week, sep = "."),
                            Site = Site,
                            Year = Year,
                            CanolaLine = CanolaLine,
                            Week = Week,
                            stringsAsFactors = F)
rownames(can.sam) <- Sample.Code

# Merge sample information with taxonomy
identical(rownames(can.sam), rownames(can.asvs))
can.tax.abu <- merge(can.sam, can.asvs, by = 0)
can.tax.abu$Row.names <- NULL

# Subset to weeks of interest
can.tax.abu1 <- can.tax.abu[can.tax.abu$Week %in% c("03","06","09"),]

# Calculate relative abundances
can.root.prop <- decostand(can.tax.abu1[,-c(1:5)], "total", MARGIN = 1)*100 
range(rowSums(can.root.prop)) # Check to see if the math was right (all should sum to 100)

# Bind the taxonomy and relative abundances
can.root.prop1 <- cbind.data.frame(can.tax.abu1[,c(1:5)], can.root.prop[rownames(tax.63)])
identical(names(can.root.prop1)[6:68], rownames(tax.63))
can.root.prop2 <- t(can.root.prop1[,c(6:68)])
colnames(can.root.prop2) <- can.root.prop1$SampleID
can.root.prop2 <- cbind.data.frame(tax.63, can.root.prop2)

# Sum relative abundances by genus
can.root.prop3 <- can.root.prop2[,c(6,8:ncol(can.root.prop2))] %>%
  group_by(genus) %>%
  summarise_all(list(~ sum(.)))

# Make a df with sample IDs as rows, weeks and taxonomy in columns
can.root.prop3.t <- t(can.root.prop3[,-1])
colnames(can.root.prop3.t) <- can.root.prop3$genus
can.root.prop3.t <- data.frame(Week = str_sub(rownames(can.root.prop3.t),start = -2), can.root.prop3.t)

# Calculate mean abundances for each taxon in each week
can.root.prop4 <- can.root.prop3.t %>%
  group_by(Week) %>%
  summarise_all(list(~ mean(.)))
can.root.prop4.t <- t(can.root.prop4[-1])
colnames(can.root.prop4.t) <- c(3,6,9)

# Add canola line information and calculate mean abundances for each line in each week
split.df <- unlist(strsplit(rownames(can.root.prop3.t),"[.]"))
CanolaLine <- split.df[seq(3,length(split.df),4)]
can.root.prop3.t1 <- can.root.prop3.t
can.root.prop3.t1 <- data.frame(CanolaLine, can.root.prop3.t)
can.root.prop3.t1$Week <- NULL
can.root.prop5 <- can.root.prop3.t1 %>%
  group_by(CanolaLine) %>%
  summarise_all(list(~ mean(.)))
can.root.prop5.t <- t(can.root.prop5[-1])
colnames(can.root.prop5.t) <- can.root.prop5$CanolaLine     # Add line as column names
can.root.prop5.t <- can.root.prop5.t[,c(1,12,2:11,13:16)]   # Reorder columns
colnames(can.root.prop5.t)[12] <- "DH27298"                 # Rename NAM-48

# Make relative abundance df sorted by mean abundance of each taxon
genus.mean <- data.frame(sort(rowMeans(can.root.prop4.t), decreasing = TRUE))
names(genus.mean) <- "mean"
can.root.prop4.t <- as.matrix(can.root.prop4.t[rownames(genus.mean),])
rownames(can.root.prop4.t) <- gsub("\\.","-", rownames(can.root.prop4.t))         # Replace periods with dashes in rownames
rownames(can.root.prop4.t)[rownames(can.root.prop4.t) == "X1959-1"] <- "1959-1"   # Tidy row names
rownames(can.root.prop4.t)[rownames(can.root.prop4.t) == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Allorhizobium-N-P-R"
rownames(can.root.prop4.t)[rownames(can.root.prop4.t) == "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia-C-P"
rownames(can.root.prop4.t)[rownames(can.root.prop4.t) == "uncultured-forest-soil-bacterium"] <- "uncultured forest soil bacterium"

# Make a color palette
col.pals <- rev(c("#b0005a","#ff90a9","#920032","#600015","#fa585c",
                  "#de4f37","#651000","#b13509","#ffae74","#a65a00",
                  "#ffc458","#edb030","#9c843f","#766900","#bacf41",
                  "#c7df68","#457700","#6ecd58","#004905","#96e78c",
                  "#006717","#01ca86","#1eefbf","#0067b5","#577cc0",
                  "#006cd1","#93a7ff","#452190","#ce8dff","#9d4bbe",
                  "#3d114e","#f4b5ff","#a72292","#8b0068","#f19ccd",
                  "#bf167f","#740042","#f2468d"))

# Sort relative abundances by means and tidy up genera names
can.root.prop5.t <- as.matrix(can.root.prop5.t[rownames(genus.mean),])            # Sort
rownames(can.root.prop5.t) <- gsub("\\.","-", rownames(can.root.prop5.t))         # Replace periods with dashes
rownames(can.root.prop5.t)[rownames(can.root.prop5.t) == "X1959-1"] <- "1959-1"   # Tidy names
rownames(can.root.prop5.t)[rownames(can.root.prop5.t) == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Allorhizobium-N-P-R"
rownames(can.root.prop5.t)[rownames(can.root.prop5.t) == "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia-C-P"
rownames(can.root.prop5.t)[rownames(can.root.prop5.t) == "uncultured-forest-soil-bacterium"] <- "uncultured forest soil bacterium"

## Create stacked barplots for developmental stage and line

jpeg("FigS7.jpg", width = 7.75, height = 6.5, units = "in", res = 300)
par(mfrow = c(1,2), mar = c(0.7,0.6,0.5,1.5), oma = c(5, 2.5, 0.5, 12), xpd = T)
# Developmental stage
bp1 <- barplot(can.root.prop4.t, col = col.pals, border = NA, width = 2, ylim = c(0,30), axes = F, xaxt = "n")
axis(side = 1, at = bp1, labels = NA, lwd = 1, tck = -0.02)
axis(side = 1, at = bp1, labels = c("Leaf","Anthesis","Seed"), cex.axis = 0.9)
axis(2, cex.axis = 0.9)
mtext("Mean relative abundance (%)", side = 2, line = 2.2, cex = 0.9)
legend("topleft", legend = substitute(paste(bold("a"), " Developmental stage")), bty = "n", inset = c(-0.05,0), cex = 0.9)
box()

# Line
par(mar=c(0.7,0.6,0.5,0.25), xpd=NA)
bp2 <- barplot(can.root.prop5.t, col = col.pals, border = NA, width = 2, ylim = c(0,30), axes = F, xaxt = "n")
axis(side = 1, at = bp2, labels = NA, lwd = 1, tck = -0.02)
axis(side = 1, at = bp2, labels = colnames(can.root.prop5.t), cex.axis = 0.85, las = 2)
axis(2, cex.axis = 0.9)
legend("topleft", legend = substitute(paste(bold("b"), " Cultivar")), bty = "n", inset = c(-0.05,0), cex = 0.9)
legend("right",inset = c(-1.1,2.5), fill = rev(col.pals), legend=rev(rownames(can.root.prop5.t)), border = NA, bty = "n", text.font = 3, cex = 0.87, y.intersp = 0.82)
box()
dev.off()
