a <- read.csv("sample_info.txt", sep='\t', header=T)
head(a$sample_title)

unlist(strsplit(as.character(a$sample_title),"_", fixed=TRUE))
matrix(unlist(strsplit(as.character(a$sample_title),"_", fixed=TRUE)),byrow=T,ncol=4)

## Split Sample_accession into different parts and create columns for each of them in the data frame

m = matrix(unlist(strsplit(as.character(a$sample_title),"_", fixed=TRUE)),byrow=T,ncol=4)

##double knockout DKO mice, knockin Kin mice, TRF transcriptional regulation factor
a$TRF = m[,1]
a
m
#a$Mice = m[,2]
#a$Marker = m[,3]
#a$Rep = m[,4]