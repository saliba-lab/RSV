### the scenic results are loaded here
### this data can be accessed in the cloned directory
binary <- readRDS("int/4.1_binaryRegulonActivity.Rds")
dim(binary)
binary_NonDupl <- readRDS("int/4.2_binaryRegulonActivity_nonDupl.Rds")
dim(binary_NonDupl)
selection <- readRDS("int/4.3_regulonSelections.Rds")
selection_corr <- selection$corr
length(selection_corr)


### Only regulons with positive correlation with other regulons are selected 
selection_corr <- binary_NonDupl[selection_corr, ]
selection_corr <- as.data.frame(selection_corr)
dim(selection_corr)


### Naive, bystander and infected are seperated
cells <- colnames(selection_corr)
Sample <- gsub("[[:upper:]]","", cells)
Sample <- gsub("\\.", "", Sample)
Sample <- gsub("1", "naive", Sample)
Sample <- gsub("2", "bystander", Sample)
Sample <- gsub("3", "infected", Sample)
Sample <- gsub("4", "naive", Sample)
Sample <- gsub("5", "bystander", Sample)
Sample <- gsub("6", "infected", Sample)
Sample <- factor(Sample)



DataSet <- gsub("[[:upper:]]","", cells)
DataSet <- gsub("\\.", "", DataSet)
DataSet <- gsub("1", "ExpSecond", DataSet)
DataSet <- gsub("2", "ExpSecond", DataSet)
DataSet <- gsub("3", "ExpSecond", DataSet)
DataSet <- gsub("4", "ExpFirst", DataSet)
DataSet <- gsub("5", "ExpFirst", DataSet)
DataSet <- gsub("6", "ExpFirst", DataSet)
DataSet <- factor(DataSet)



pdata <- data.frame(cells, Sample, DataSet)



### A heatmap of the whole data is created... In this heatmap, the viral load is not used.
library(ComplexHeatmap)
colors = structure(c("black", "white"), names = c("1", "0"))
WholeHeatMap <- Heatmap(selection_corr, col = colors, cluster_rows = TRUE, cluster_columns = TRUE, name = "legend", show_column_names = FALSE)
WholeHeatMap



t <- row_order(WholeHeatMap)
t <- as.data.frame(t)
t <- t$c.36L..27L..9L..53L..6L..21L..34L..51L..22L..66L..61L..44L..41L..
row <- rownames(selection_corr)

order <- c(1:68)
for (i in order) {
  z <- t[i]
  order[i] <- row[z]
}

order <- rev(order)


### Next a seperate heatmap for naive, bystander and infected are created. 
### Here the viral load is effective,
### at the end, three heatmaps are attached together

### Infected
###############################################################################
### retrieving the percentage of virus in each cell
#subsetting the cells
### This object can be accessed in the cloned directory
load("/home/ehsan/Desktop/Sibylle/Sibylle/E8E9E10B9B10B11/E8E9E10B9B10B11Analysed.Robj")
E8E9E10B9B10B11sub <- SubsetData(E8E9E10B9B10B11, ident.use = c(1,2,3,10,7,9))
cells <- colnames(as.matrix(E8E9E10B9B10B11sub@data))

### Note that the way that viral RNA is calculated here is different from the way it is calculated in the new data sets

data <- as.matrix(E8E9E10B9B10B11sub@data)
data <- exp(data) - 1
dim(data)
## Calculating total viral mRNA
virus <- data[18693:18704, ]
virus.tot <- as.data.frame(colSums(virus))
cell.tot <- as.data.frame(colSums(data))
cell.percent <- (virus.tot/cell.tot)*100
colnames(virus.tot) <- "Total_viral_mRNA"
colnames(cell.percent) <- "Virus_percentage"
t <- gsub("-", ".", rownames(cell.percent))
rownames(cell.percent) <- t


###############################################################################


infected <- pdata[pdata$Sample == "infected", ]
dim(infected)
infMatrix <- selection_corr[order ,as.character(infected$cells)]
dim(infMatrix)



## producing the heatmap
library(ComplexHeatmap)
library(circlize) ## necessary for colour assignment
#############################
#Retrieving viral percentage information
#f <- as.character(infected$cells)
#cell.percent <- cell.percent[f,]
#############################
## Producing the annotation files
df = data.frame(infected$DataSet, infected$Sample)
ha = HeatmapAnnotation(df = df, col = list(infected.DataSet = c("ExpFirst" =  "red", "ExpSecond" = "blue"), infected.Sample = c("naive" =  "green", "bystander" = "steelblue4", "infected" = "brown")))

#############################
colors = structure(c("black", "white"), names = c("1", "0"))


png(file = "Infected.png", width = 434.5, height = 1080)
InfHeatMap <- Heatmap(infMatrix, col = colors, cluster_rows = FALSE, cluster_columns = TRUE, name = "legend", show_column_names = FALSE, top_annotation = ha)
print(InfHeatMap)
dev.off()


### bystander

##############################################################################



bystander <- pdata[pdata$Sample == "bystander", ]
dim(bystander)
byMatrix <- selection_corr[ order,as.character(bystander$cells)]
dim(byMatrix)



## producing the heatmap
library(ComplexHeatmap)
library(circlize) ## necessary for colour assignment
#############################
## Producing the annotation files
df = data.frame(bystander$DataSet, bystander$Sample)
ha = HeatmapAnnotation(df = df, col = list(bystander.DataSet = c("ExpFirst" =  "red", "ExpSecond" = "blue"), bystander.Sample = c("naive" =  "green", "bystander" = "steelblue4", "infected" = "brown")))

#############################
colors = structure(c("black", "white"), names = c("1", "0"))


png(file = "bystander.png", width = 544, height = 1080)
ByHeatMap <- Heatmap(byMatrix, col = colors, cluster_rows = FALSE, cluster_columns = TRUE, name = "legend", show_column_names = FALSE, show_row_names = FALSE ,top_annotation = ha)
print(ByHeatMap)
dev.off()


### Naive


naive <- pdata[pdata$Sample == "naive", ]
dim(naive)
naiveMatrix <- selection_corr[ order,as.character(naive$cells)]
dim(naiveMatrix)



## producing the heatmap
library(ComplexHeatmap)
library(circlize) ## necessary for colour assignment
#############################
## Producing the annotation files
df = data.frame(naive$DataSet, naive$Sample)
ha = HeatmapAnnotation(df = df, col = list(naive.DataSet = c("ExpFirst" =  "red", "ExpSecond" = "blue"), naive.Sample = c("naive" =  "green", "bystander" = "steelblue4", "infected" = "brown")))

#############################
colors = structure(c("black", "white"), names = c("1", "0"))


png(file = "naive.png", width = 562.5, height = 1080)
NaiveHeatMap <- Heatmap(naiveMatrix, col = colors, cluster_rows = FALSE, cluster_columns = TRUE, name = "legend", show_column_names = FALSE, show_row_names = FALSE , top_annotation = ha)
print(NaiveHeatMap)
dev.off()



### combining three heatmaps 
class(InfHeatMap)
class(ByHeatMap)
class(NaiveHeatMap)


#### In this heatmap infected cells are not ordered based on viral load
NaiveHeatMap + ByHeatMap + InfHeatMap

### Next we produce a heatmap that infected cells are ordered based on viral load
### Infected heatmap with new criteria


infected <- pdata[pdata$Sample == "infected", ]
dim(infected)
infMatrix <- selection_corr[order ,as.character(infected$cells)]
dim(infMatrix)



## producing the heatmap
library(ComplexHeatmap)
library(circlize) ## necessary for colour assignment
#############################
#Retrieving viral percentage information
f <- as.character(infected$cells)
cell.percent2 <- cell.percent[f,]
#############################
## Producing the annotation files
df = data.frame(infected$DataSet, infected$Sample, cell.percent2)
ha = HeatmapAnnotation(df = df, col = list(infected.DataSet = c("ExpFirst" =  "red", "ExpSecond" = "blue"), infected.Sample = c("naive" =  "green", "bystander" = "steelblue4", "infected" = "brown"), cell.percent2 = colorRamp2(c(0, max(cell.percent2)), c("yellow", "red"))))

#############################
colors = structure(c("black", "white"), names = c("1", "0"))


#png(file = "Infected.png", width = 434.5, height = 1080)
InfHeatMap <- Heatmap(infMatrix, col = colors, cluster_rows = FALSE, cluster_columns = TRUE, name = "legend", show_column_names = FALSE, top_annotation = ha)
print(InfHeatMap)
#dev.off()










########################################################
########################################################


bystander <- pdata[pdata$Sample == "bystander", ]
dim(bystander)
infMatrix <- selection_corr[order ,as.character(bystander$cells)]
dim(infMatrix)



## producing the heatmap
library(ComplexHeatmap)
library(circlize) ## necessary for colour assignment
#############################
#Retrieving viral percentage information
f <- as.character(bystander$cells)
cell.percent2 <- cell.percent[f,]
#############################
## Producing the annotation files
df = data.frame(bystander$DataSet, bystander$Sample, cell.percent2)
ha = HeatmapAnnotation(df = df, col = list(bystander.DataSet = c("ExpFirst" =  "red", "ExpSecond" = "blue"), bystander.Sample = c("naive" =  "green", "bystander" = "steelblue4", "infected" = "brown"), cell.percent2 = colorRamp2(c(0, max(cell.percent2)), c("yellow", "red"))))

#############################
colors = structure(c("black", "white"), names = c("1", "0"))


#png(file = "bystander.png", width = 434.5, height = 1080)
ByHeatMap <- Heatmap(infMatrix, col = colors, cluster_rows = FALSE, cluster_columns = TRUE, name = "legend", show_column_names = FALSE, show_row_names = FALSE, top_annotation = ha)
print(ByHeatMap)
#dev.off()


########################################################
# Here the final heatmap is created
########################################################


NaiveHeatMap + ByHeatMap + InfHeatMap


