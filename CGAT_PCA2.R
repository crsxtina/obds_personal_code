library(reshape2)
library(ggplot2)
library(pheatmap)
## Dataset 1
## Look at the data
system("head /Users/cristina/Cristina/DPhil/CGAT/OBDS_R_week/joined_final.txt")

## Load the data into R
a = read.csv("/Users/cristina/Cristina/DPhil/CGAT/OBDS_R_week/joined_final.txt", header=T, row.names = 1, sep=' ')
sample_info = read.csv("/Users/cristina/Cristina/DPhil/CGAT/OBDS_R_week/sample_info.txt", sep='\t')

## check the size of the data
head(a)
head(sample_info)
dim(a) # 12 samples
colnames(a)

## Split Sample_accession into different parts and create columns for each of them in the data frame
m = matrix(unlist(strsplit(as.character(a$sample_title),"_",fixed=TRUE)),byrow=T,ncol=4)
##double knockout DKO mice, knockin Kin mice, TRF transcriptional regulation factor
a$TRF = m[,1]
a$Mice = m[,2]
a$Marker = m[,3]
a$Rep = m[,4]

## create a wide data format
#data_wide <- dcast(a, TRF+Mice+Rep ~ Marker, value.var="read_count")
#data_wide # 2 types of mice, 3 replicates each, 2 markers CD4 and CD8

## run PCA
#project.pca <- prcomp(data_wide[,4:5], center=TRUE, scale=TRUE)
#summary(project.pca)

## Determine the proportion of variance of each component
#project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100
#barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))
## plot PCA components
#ggplot(as.data.frame(project.pca$x), aes(PC1,PC2,colour=data_wide$Mice))+geom_point()+ggtitle("Principal components analysis")
#par(cex=1.0, cex.axis=0.8, cex.main=0.8)
#pairs(project.pca$x[,1:2], col=as.numeric(as.factor(data_wide$Mice)), pch=19, main="Principal components analysis bi-plot\nPCs 1-2") 

##################
# Second dataset

# Look at the data
system("head /Users/vi/Dropbox/WIMM/CGAT/data_david/symonds_counts_UCSC_refGene_mm10.txt")
# Load the data into R
a = read.table("/Users/vi/Dropbox/WIMM/CGAT/data_david/symonds_counts_UCSC_refGene_mm10.txt", head=T, row.names=1)
# check the size of the data
dim(a) # 12 samples, 24453 genes
# Split Sample_accession into different parts and create columns for each of them in the data frame
m = as.data.frame(matrix(unlist(strsplit(names(a),"_",fixed=TRUE)),byrow=T,ncol=4))
#double knockout DKO mice, knockin Kin mice, TRF transcriptional regulation factor
names(m) = c("TRF", "Mice", "Marker", "Rep")

# total expression per sample
colSums(a)
# total expression per gene
rowSums(a)
hist(as.numeric(rowSums(a)),
     breaks = 100, main = "Expression sum per gene",
     xlab = "Sum expression")
abline(v=median(as.numeric(rowSums(a))),col="red")

# remove genes not expressed in any sample
keep_feature <- rowSums(a) > 0
a <- a[keep_feature, ]

# density plot of raw read counts (log10)
logcounts <- log(a[,1]+1,10) 
d <- density(logcounts)
plot(d,xlim=c(1,8),main="",ylim=c(0,.45),xlab="Raw read counts per gene (log10)", ylab="Density")
for (s in 2:ncol(a)){
  logcounts <- log(a[,s],10) 
  d <- density(logcounts)
  lines(d)
}

a_log = log(a+1)
hist(as.numeric(rowSums(a_log)),
     breaks = 100, main = "log Expression sum per gene",
     xlab = "Sum expression")
abline(v=median(as.numeric(rowSums(a_log))),col="red")



# select data for the 11000 most highly expressed genes
select = order(rowMeans(a), decreasing=TRUE)[1:11000]
select = order(rowMeans(a), decreasing=TRUE)[1:100]
highexprgenes_counts <- a[select,]
highexprgenes_logcounts <- a_log[select,]

hist(as.numeric(rowSums(highexprgenes_logcounts)),
     breaks = 100, main = "log Expression sum per gene",
     xlab = "Sum expression")
abline(v=median(as.numeric(rowSums(highexprgenes_logcounts))),col="red")


par(mfrow=c(2,2)) ## if you want to have multiple plots on the same window.


# heatmap with sample name on X-axis
heatmap(as.matrix(highexprgenes_logcounts), col=topo.colors(50), margin=c(10,6))
# heatmap with condition group as labels
#colnames(highexprgenes_logcounts)<- m$Mice
# plot
#heatmap(as.matrix(highexprgenes_logcounts), col = topo.colors(50), margin=c(10,6))


data_for_PCA <- t(highexprgenes_logcounts)
dim(data_for_PCA)

############### Back to slides #################
data_for_PCA_scaled = scale(data_for_PCA,center = TRUE, scale = FALSE)

## by hand
ctmp = t(data_for_PCA_scaled)%*%(data_for_PCA_scaled) ## X^T%*%X matrix
ceig = eigen(ctmp) ## eigenvalue decomposition
ceigenvectors = ceig$vectors ## eigenvectors
ceigenvalues = ceig$values ## eigenvalues - diagonal of L matrix from slides
# projection of data onto PC space
cPC= (data_for_PCA_scaled)%*%ceigenvectors%*%diag(1/sqrt(ceigenvalues))
# percentage of variance explained by each PC
print(round(ceigenvalues/sum(ceigenvalues) * 100, digits = 2))
par(mfrow=c(2,1))
plot(ceigenvalues/sum(ceigenvalues) * 100, xlab="PC",ylab="% Variance explained")
# cumulative % of variance explained by each PC
round(cumsum(ceigenvalues)/sum(ceigenvalues) * 100, digits = 2) 
plot(cumsum(ceigenvalues/sum(ceigenvalues) * 100), xlab="PC",ylab="Cumulative % Variance explained")
#
par(mfrow=c(1,1))
plot(cPC[,1],cPC[,2], xlab="PC1",ylab="PC2", main="PC manual")


# using prcomp
bpc = prcomp(data_for_PCA, center=TRUE, scale=FALSE) ## stats package 
beigenvalues = bpc$sdev^2 # eigenvalues
beigenvectors = bpc$rotation # eigenvectors
par(mfrow=c(2,1))
plot(beigenvalues/sum(beigenvalues) * 100,xlab="PC",ylab="% Variance explained") 
plot(cumsum(beigenvalues)/sum(beigenvalues) * 100, xlab="PC",ylab="Cumulative % Variance explained")

plot(bpc$x[,1]/(bpc$sdev[1]*sqrt(12)),bpc$x[,2]/(bpc$sdev[2]*sqrt(12)), xlab="PC1",ylab="PC2", main="PC prcomp")
library(ggfortify)
autoplot(bpc, main = "PC prcomp, autoplot")


#source("https://bioconductor.org/biocLite.R")
#biocLite("pcaExplorer")
library(pcaExplorer)
# extract genes with the highest/lowest loadings per PC
hi_loadings(bpc, whichpc = 1, topN = 10, exprTable = NULL,
            annotation = NULL, title = "Top/bottom loadings - ")

hi_loadings(bpc, whichpc = 2, topN = 10, exprTable = NULL,
            annotation = NULL, title = "Top/bottom loadings - ")

s = data_for_PCA[,which(colnames(data_for_PCA)=="213742")]
s2 = cbind(s,m)
ggplot(s2, aes(Mice,s, fill=Marker))+ geom_boxplot(position="dodge") +geom_point(alpha=0.6, aes(group=Marker), data=s2, position = position_dodge(width=0.75))+ylab("log Expression")+xlab("")+ggtitle("Gene 213742")

s = data_for_PCA[,which(colnames(data_for_PCA)=="21414")]
s2 = cbind(s,m)
ggplot(s2, aes(Mice,s, fill=Marker))+ geom_boxplot(position="dodge") +geom_point(alpha=0.6, aes(group=Marker), data=s2, position = position_dodge(width=0.75))+ylab("log Expression")+xlab("")+ggtitle("Gene 21414")

# doing the above with custom code
whichpc = 1
topN = 10
geneloadings_sorted <- sort(bpc$rotation[, whichpc])
geneloadings_extreme <- c(tail(geneloadings_sorted, topN), head(geneloadings_sorted, topN))
barplot(geneloadings_extreme, las = 2, col = c(rep("steelblue", topN), rep("coral", topN)), main = paste0("PC", whichpc," loadings"))

# plotting a heatmap of expression values for the top topN genes
s = data_for_PCA[,which(colnames(data_for_PCA)%in%names(geneloadings_extreme))]
s2 = cbind(s,m)
data_narrow <- melt(s2, id =names(m))
ggplot(data_narrow, aes(Mice,variable))+geom_tile(aes(fill=value))+facet_grid(.~Marker)+ylab("Gene")+scale_fill_gradientn(colours = rainbow(20))  #+ scale_fill_gradient2( low = "blue", high = "red", na.value="black", name = "" )
pheatmap(s)

# doing the same for PC2
whichpc = 2
topN = 10
geneloadings_sorted <- sort(bpc$rotation[, whichpc])
geneloadings_extreme <- c(tail(geneloadings_sorted, topN), head(geneloadings_sorted, topN))
barplot(geneloadings_extreme, las = 2, col = c(rep("steelblue", topN), rep("coral", topN)), main = paste0("PC", whichpc," loadings"))
s = data_for_PCA[,which(colnames(data_for_PCA)%in%names(geneloadings_extreme))]
s2 = cbind(s,m)
data_narrow <- melt(s2, id =names(m))
ggplot(data_narrow, aes(Mice,variable))+geom_tile(aes(fill=value))+facet_grid(.~Marker)+ylab("Gene")+scale_fill_gradientn(colours = rainbow(20))
pheatmap(s)



###
# plot gene expression on top of PCA plot
s = data_for_PCA[,which(colnames(data_for_PCA)=="21414")]
logExpression = as.numeric(s)
shape = paste(m$Mice,m$Marker,sep="-")
ggplot(bpc$x, aes(PC1,PC2,colour=logExpression)) +geom_point(aes(shape=shape)) +ggtitle("Gene 21414 expression")

s = data_for_PCA[,which(colnames(data_for_PCA)=="12010")]
logExpression = as.numeric(s)
shape = paste(m$Mice,m$Marker,sep="-")
ggplot(bpc$x, aes(PC1,PC2,colour=logExpression)) +geom_point(aes(shape=shape)) +ggtitle("Gene 12010 expression")






####### Back to slides ###########
# other ways to perform tsne - slow so not run here for the scRNA-seq data
library("Rtsne")
iris_matrix <- as.matrix(data_for_PCA)
set.seed(42) # Set a seed if you want reproducible results
tsne_out <- Rtsne(iris_matrix, dims = 2, perplexity = 3) # Run TSNE
# Show the objects in the 2D tsne representation
#plot(tsne_out$Y,col=m$Mice, pch=as.numeric(m$Marker))

dat = tsne_out$Y
colnames(dat) = c("tSNE1","tSNE2")
ggplot(dat, aes(tSNE1,tSNE2))+geom_point() +ggtitle("tSNE")

s = data_for_PCA[,which(colnames(data_for_PCA)=="21414")]
logExpression = as.numeric(s)
shape = paste(m$Mice,m$Marker,sep="-")
ggplot(dat, aes(tSNE1,tSNE2,colour=logExpression)) +geom_point(aes(shape=shape)) +ggtitle("Gene 21414 expression")



s = data_for_PCA[,which(colnames(data_for_PCA)=="12010")]
logExpression = as.numeric(s)
s2 = cbind(s,m)
s2$cluster = paste(s2$Mice,":",s2$Marker, sep="")
ggplot(s2,aes(cluster, s,fill=cluster))+geom_boxplot()+geom_point()

s = data_for_PCA[,which(colnames(data_for_PCA)=="213742")]
logExpression = as.numeric(s)
s2 = cbind(s,m)
s2$cluster = paste(s2$Mice,":",s2$Marker, sep="")
ggplot(s2,aes(cluster, s,fill=cluster))+geom_boxplot()+geom_point()


#library(edgeR)
#library(limma)
#design <- model.matrix(~0+m$Mice+m$Marker)
## substitute "group" from the design column names
#colnames(design)<- c("DKO","Kin","CD8")
# check your design matrix
#fit <- lmFit(t(data_for_PCA),design)
#
#cont.matrix <- makeContrasts(DKO-Kin,levels=design)
#fit <- contrasts.fit(fit, cont.matrix)
#fit <- eBayes(fit)
#options(digits=3)

# check the output fit
#dim(fit)
#mypval=0.01
#myfc=3

# get the coefficient name for the comparison  of interest
#colnames(fit$coefficients)
#mycoef="DKO - Kin"
# get the output table for the 10 most significant DE genes for this comparison
#topTable(fit,coef=mycoef)
#limma.res <- topTable(fit,coef=mycoef,n=dim(fit)[1])

# get significant DE genes only (adjusted p-value < mypval)
#limma.res.pval <- topTable(fit,coef=mycoef,n=dim(fit)[1],p.val=mypval)
#dim(limma.res.pval)

# most differentially expressed genes
#s = data_for_PCA[,which(colnames(data_for_PCA)%in%rownames(limma.res.pval))]
#s2 = cbind(s,m)
#pheatmap(s)

#  http://monashbioinformaticsplatform.github.io/RNAseq-DE-analysis-with-R/RNAseq_DE_analysis_with_R.html

############# Back to slides #####################

### Clustering
# kmeans clustering https://www.r-bloggers.com/k-means-clustering-in-r/
Cluster <- kmeans(data_for_PCA[,1:2], 3, nstart = 20)
Cluster$cluster <- as.factor(Cluster$cluster)
# get cluster means 
aggregate(data_for_PCA[,1:2],by=list(Cluster$cluster),FUN=mean)

ggplot(data_for_PCA[,1:2], aes(`13627`,`13629`, color = Cluster$cluster, shape=shape)) + geom_point()

Cluster <- kmeans(data_for_PCA[,1:2], 4, nstart = 20)
Cluster$cluster <- as.factor(Cluster$cluster)
ggplot(data_for_PCA[,1:2], aes(`13627`,`13629`, color = Cluster$cluster, shape=shape)) + geom_point()

par(mfrow=c(1,1))
measure= NULL
for (k in 1:11){
  Clust <- kmeans(data_for_PCA[,1:2],k,nstart=12)
  measure = c(measure,Clust$betweenss/Clust$totss)
}
plot(measure, xlab="k", ylab="between_SS / total_SS", type="b", pch=19)


# Important
#Within cluster sum of squares by cluster:
#  [1] 13.05769 16.29167  2.02200
#(between_SS / total_SS =  94.3 %)

#totss, tot.withinss, betweenss info on full clustering

#clusters 
Cluster$cluster
# their centers
Cluster$centers
# numbers of points per cluster
Cluster$size




# https://www.r-bloggers.com/finding-optimal-number-of-clusters/

### NOT run on full data!
##################

#hierarchical clustering

d_pca <- dist(data_for_PCA) # euclidean by default
hc_pca <- hclust(d_pca, method = "complete") ## hierarchical cluster analysis
plot(hc_pca)
groups <- cutree(hc_pca,k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(hc_pca, k=5, border="red")

library(pheatmap)  
pheatmap(data_for_PCA[1:12,1:10])

pheatmap(data_for_PCA)

# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
fit <- pvclust(t(data_for_PCA), method.hclust="ward",
               method.dist="euclidean")
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)



# Model Based Clustering
#library(mclust)
#fit <- Mclust(t(data_for_PCA))
#plot(fit) # plot results 
#summary(fit) # display the best model



############## Exercise ################


#Download the two files from the links below
#https://www.dropbox.com/s/au4a4snlhxa39fq/GTEx_subset_1000_subset.txt?dl=0
#https://www.dropbox.com/s/m6hhef9nv0vwka5/GTEx_subset_1000_subset_info.txt?dl=0

#File 1 contains a read count matrix (rows=genes, cols=samples)
#File 2 includes information on the tissues of each sample

#1. Perform PCA. How many PCs do you think you should keep for follow up analysis?
#2. Use the tissue information for plotting to see which tissues cluster together
#3. Investigate which tissues (if any) the first few principal components are separating
#4. Find the top genes correlated to PC1/PC2 and plot the Expression values per tissue
#5. Compare PCA to tSNE. Unlike the previous example, these should look very different now!
#6. How would you cluster the data? How many clusters would you choose?
#7. optional - according to each cluster, if familiar with differential expression analysis, find the most highly differentially expressed genes and for each such gene colour tSNE, PCA plots with that gene's expression values 




