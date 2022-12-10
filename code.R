library(readxl)
library(dplyr) 
library(data.table)
library(tibble)
library(limma)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)
library(plotly)
library(pracma)
library(scales)
library(purrr)

#set working directory
setwd("~/genetic_project")

#read expression data from excel 
exp_df <- read_excel("data.xlsx", sheet = 1, skip = 3)
View(exp_df)

# column names to vector
col_names <- colnames(exp_df)

#set dataframe rowname
exp_df <- column_to_rownames(exp_df, var ="...1")

#grouping by protocol
gr1 <- map_chr(col_names[2:length(col_names)],function(x) substr(x, 0, 3))

#group by naive or prime
vec <- grepl("\\.N",col_names[2:length(col_names)])
fac <- factor(vec, labels = c('P', 'N') )

#check data range
max(exp_df)
min(exp_df)

#### draw plotbox , showing that data is normalized
pdf("boxplot.pdf", width = 200, height = 50)
par(cex.axis = 2)
boxplot(exp_df)
dev.off()


### Correlation Heatmap
pdf("CorHeatmap.pdf", width = 200, height = 200)
pheatmap(cor(exp_df), labels_row = colnames(exp_df), labels_col = colnames(exp_df), color = bluered(256), border_color = NA)
dev.off()


## principal component analysis for all cells
exp_df.scale <- t(scale(t(exp_df), scale = F))
pc1 <- prcomp(exp_df.scale)
pdf("PC_scaled.pdf")
plot(pc1)
plot(pc1$x[,1:2])
dev.off()

#plot PCA1 in file
pcr1 <- data.frame(pc1$rotation[,1:2], Batch = gr1, Sample = fac)
pdf("PCA.pdf", width = 10, height = 10)
cols <- hue_pal()(2)
ggplot(pcr1,aes(PC1, PC2,shape = Batch , color = Sample))+  scale_shape_manual(values=seq(1,10)) + geom_point(size = 5) +  scale_color_manual(values = rev(cols)) + theme_bw()
    dev.off()

#selecting naive cells from dataframe
naive_df <- exp_df[vec]

#PCA for naive cells 
naive_df.scale <- t(scale(t(naive_df), scale = F))
pc2 <- prcomp(naive_df.scale)
pdf("PC_scaled_naive.pdf")
plot(pc2)
plot(pc2$x[,1:2])
dev.off()

cluster <- kmeans(t(naive_df.scale), 2)$cluster

pcr2 <- data.frame(pc2$rotation[,1:3], Cluster = cluster)
pdf("PCA2.pdf", width = 10, height = 10)
ggplot(pcr2,aes(PC1, PC2, color = Cluster)) + geom_point(size = 5)  + theme_bw()
dev.off()

