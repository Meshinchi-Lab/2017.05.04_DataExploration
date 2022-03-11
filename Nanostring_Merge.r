#Jenny Smith

#April 21, 2017

#Purpose: Merge dataframes for Nanostring


# setwd("H:/Nanostring_Analysis/data/")
setwd(file.path(TARGET,"RNA/Nanostring/2017.05.04_DataExploration"))

#Reference mapping file
ref <- read.csv("H:/reference_mapping-files/TARGET_AML_current_asof_june30_2016_FINAL.csv", stringsAsFactors = FALSE)

head(ref[,1:10])

ID_map <- ref[,1:2]

#Read in the data
Apr04 <- read.csv("AML_EXPR_C5117_170406_01-03.csv", stringsAsFactors = FALSE)

Apr05 <- read.csv("AML_EXPR_C5117_170410_04-06.csv", stringsAsFactors = FALSE)

dim(Apr04) #272  39
dim(Apr05) #272  39

head(Apr04)
head(Apr05)

IDs_04 <- Apr04[2,][-(1:3)] #$remove rownames filename, ID, and X


IDs_05 <- Apr05[2,][-(1:3)]

#How many IDs are in both sets
length(which(IDs_04 %in% IDs_05 ))
length(which(IDs_05 %in% IDs_04 ))

#genes in same order?
identical(Apr04$File.Name, Apr05$File.Name) #TRUE
identical(Apr04$X, Apr05$X) #True
identical(Apr04$X.1, Apr05$X.1) #True

#Remove the extra repetitiv info for now


#merge the sets
colnames(Apr04)[1:10]
merged <- cbind(Apr04, Apr05) #will have duplicates of teh filename,X,and X.1 

#characterize the merge
dim(merged) #272 78 
head(merged)

#remove the duplicated cols. Did not use merge, becuase that caused issues with formatting due to top 15 rows.
merged <- merged[,-(39:41)] 

dim(merged) #272 75




write.csv(merged, file = "TARGET_AML_nanostring_merged_21April2017.csv", row.names = FALSE)
save(merged, file = "TARGET_AML_nanostring_merged_21April2017.RData")



