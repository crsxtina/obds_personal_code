our_data <- read.csv("sample_info.txt", sep='\t', header=T)
class(our_data)
str(our_data)
dim(our_data)
summary(our_data)
names(our_data)
colnames(our_data)
row5 <- our_data[5,]
entry <- our_data[4, "sample_title"]
entry

entry_samples <- as.data.frame(our_data$sample_title)
head(entry_samples)
class(entry_samples)
entry4_samples <- entry_samples[4,]
entry4_samples

new_col <- our_data$mill_reads <- our_data$read_count/10000000
head(new_col)

# main is for title
# xlab or names.
# las=2 to put the x labels vertically

barplot(our_data$mill_reads, main="mill_reads", names.arg = our_data$sample_title, las=2)
?barplot

athena_data <- read.csv("joined_final.txt", sep = ' ', header = T)
str(athena_data)
summary(athena_data)
class(athena_data)
dim(athena_data)
colnames(athena_data)
head(athena_data)

common_genes <- athena_data[rowSums(athena_data[-1]) > 0,] #-1 means ignore first row
head(common_genes)

max_index <- common_genes[order(rowSums(common_genes[-1])),]
head(max_index)
tail(max_index)
