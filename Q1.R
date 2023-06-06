library(vegan)
library(ecodist)
#install.packages("dplyr")                             # Install dplyr package
library("dplyr") 
library("ggplot2")

#You have to investigate the following:
#1. a. Is there a significant variation in the Evenness (Pielou Evenness) and
#overall alpha diversity (both Richness and Evenness quantified by the
#Shannon Index) between the patients across the different severity
#groups? Show the results as boxplots and indicate the appropriate tests used.

##Reading the dataset files
data1_rawcount <- read.csv("Assignment_RawCount_Species - Sheet1.csv")
data2_metadata<-read.csv("Assignment1_Metadata - Sheet1.csv")

###Question 1a ###

c<-data2_metadata$X #Retrieving first column from metadata file

df_merged <- data1_rawcount[data1_rawcount$X %in% c,] #Rows from rawcount are copied based on column X in rawcount and metadata files
#data1$WHO_severity <- data2$WHO_severity[match(data1$X, data2$X)]

#Calculating Shannon Index
df_shannon=(diversity(df_merged[-1], index="shannon")) 
df_merged['Shannon Index']=df_shannon             #Adding Shannon Index values to df_merged


#Adding WHO_severity column from metadata to df_merged
df_merged$WHO_severity <- data2_metadata$WHO_severity[match(df_merged$X, data2_metadata$X)] 

#Converting mild,moderate,critical_severe to 1,2 and 3
df_merged$WHO_severity[which(df_merged$WHO_severity=="mild")] <- 1
df_merged$WHO_severity[which(df_merged$WHO_severity=="moderate")]<-2
df_merged$WHO_severity[which(df_merged$WHO_severity=="critical_severe")] <- 3

#Calculating Evenness

S <- apply(df_merged[,-1]>0,1,sum)
df_evenness<-df_shannon/log(S)
df_merged['Evenness']<-df_evenness #Adding evenness column to df_merged

#Plotting Shannon Index vs Covid Severity
x<- df_merged$WHO_severity
y<-df_merged$`Shannon Index`
boxplot(y~ x, data = df_merged,names=c("mild", "moderate","severe"),xlab = "Covid Severity",ylab = "Shannon Index",main = "Shannon Index vs Covid Severity",col = c("green", "red", "blue"))

#Plotting Evenness vs Covid Severity
x<- df_merged$WHO_severity
y1<-df_merged$Evenness
boxplot(y1~ x, data = df_merged,names=c("mild", "moderate","severe"),xlab = "Covid Severity",ylab = "Evenness",main = "Evenness vs Covid Severity",col=c("green","red","blue"))

#After visualizing the boxplots, we can perform statistical tests to determine if there is a significant difference in the alpha diversity and evenness measures between the Covid severity groups.
#Since the data may not follow a normal distribution, we can see outliers in box plot, we can use non-parametric tests such as the Kruskal-Wallis test to compare the distributions 
#of alpha diversity and evenness measures across the different Covid severity groups.
#We can perform these tests using the stats package in R. The code for performing the tests is as follows:

# Kruskal-Wallis test for Shannon Index - overall alpha diversity 
kw_diversity <- kruskal.test(y ~ x, data = df_merged)
print(kw_diversity)

# Kruskal-Wallis test for Pielou Evenness
kw_evenness <- kruskal.test(y1 ~ x, data = df_merged)
print(kw_evenness)

#Based on the Kruskal-Wallis rank sum test results, there is no significant variation in the evenness or overall alpha diversity (both richness and evenness quantified by the Shannon index) between the patients across the different severity groups. 
#The p-values for both tests are greater than 0.05(kw_diversity=0.2486 and kw_evenness=0.1732), indicating that we cannot reject the null hypothesis that there is no significant difference between the groups.



##Question 1 b- Beta Diversity##
#b. Is there a difference in the overall abundance pattern of the
#microbiomes (beta-diversity) across the different severity groups. Show
#the results as Prinicipal Component Analysis or Principal Coordinate
#Analysis plots and indicate the significance of separation. Which test
#would you use for this?
#####################################################
##Creating a dataframe to calculate betadiversity by removing "X","Shannon Index", "WHO_severity", "Evenness" columns
data2_metadata<-read.csv("Assignment1_Metadata - Sheet1.csv")
data3_clrcount <- read.csv("Assignment1_ClrTrans_Species - Sheet1.csv")


df_merged_clr <- data3_clrcount[data3_clrcount$X %in% data2_metadata$X,]
#df_beta_div<- df_merged[,!names(df_merged) %in% c("X","Shannon Index", "WHO_severity", "Evenness")]
mat_df.mdf <- as.matrix.data.frame(df_merged_clr[-1])

###Calculating the distance matrix
div_beta=vegdist(mat_df.mdf,method="bray")
res<-pco(div_beta)

#Adding WHO_severity column from metadata to df_merged
df_merged_clr$WHO_severity <- data2_metadata$WHO_severity[match(df_merged_clr$X, data2_metadata$X)] 

#Converting mild,moderate,critical_severe to 1,2 and 3
df_merged_clr$WHO_severity[which(df_merged_clr$WHO_severity=="mild")] <- 1
df_merged_clr$WHO_severity[which(df_merged_clr$WHO_severity=="moderate")]<-2
df_merged_clr$WHO_severity[which(df_merged_clr$WHO_severity=="critical_severe")] <- 3
###Adding PC1,PC2 and WHO Severity to res_pcoa_df
res_pcoa_df <- data.frame(pcoa1 = res$vectors[,1], pcoa2 = res$vectors[,2])

res_pcoa_df <- cbind(res_pcoa_df,severity=df_merged_clr$WHO_severity)

###Plotting PCOA using data=res_pcoa_df

braycurtis_PCOA_plot <- ggplot(data = res_pcoa_df, aes(x=pcoa1, y=pcoa2,
                                                       color = severity)) + geom_point() +labs(x = "PC1",y = "PC2",title = "PCoA with BrayCurtis distances") +theme(title = element_text(size = 12)) 

braycurtis_PCOA_plot

##This will produce a PCoA plot of the gut microbiome beta-diversity, with the points colored by the severity of Covid-19. 
##We can use PERMANOVA (permutational multivariate analysis of variance) to test whether there is a 
##significant difference in the microbiome composition across the different severity groups:

pcoa_adonis <- adonis2(div_beta ~ WHO_severity, data = data2_metadata, permutations = 999)
print(pcoa_adonis)

#The p-value of 0.043 suggests that there is a significant difference 
#in the gut microbiome composition across the different severity groups 
#of COVID-19 in dataset.


#attempt2 for both question 1a and 1b

#1a -----------------------------------------------------------------------------------

library(readxl)
library(dplyr)

#Load the metadata and abundance profiles into R using the read_excel() function from the "readxl" package:

# Load metadata
metadata <- read_excel("F:/IIIT-D/HMDS/Assignment1_Metadata.xlsx")

# Load abundance profiles
#abundance <- read_excel("F:/IIIT-D/HMDS/Assignment1_ClrTrans_Species.xlsx", sheet = 1)
abundance <- read_excel("F:/IIIT-D/HMDS/Assignment_RawCount_Species.xlsx", sheet = 1)

# merge the two data frames on the 'id' column
merged_df <- merge(abundance, metadata, by = "SampleID")

# # remove rows from abundance that are not present in metadata
merged_df  <- merged_df [!is.na(merged_df $Sex),]

# print the resulting data frame
merged_df
abundance <- select(merged_df , -Sex,-Age,-BMI,-HTN,-Diabetes,-Respiratory_disease,-Heart_disease,-Renal_Disease,-Liver_Disease,-Obesity,-Malignancy,-Immunosuppressive_Disease,-Neurological_disease,-Metabolic_Disease,-Cardiovascular_Disease,-comorbidities_total,-WHO_severity)

#remove 1st column i.e sample-id
#rownames(abundance) <- abundance$SampleID
abundance <- abundance[, -1]
abundance

#Calculate alpha diversity and evenness metrics using the "vegan" package:
#install.packages("vegan")
#install.packages("permute")
library(permute)
library(vegan)

# Calculate richness
richness <- specnumber(abundance)
# Calculate Shannon diversity index
diversity <- diversity(abundance, index = "shannon")
# Calculate Pielou's evenness index
evenness <- diversity / log(richness)


#Create boxplots for each metric using the "ggplot2" package:


library(ggplot2)
# Combine metadata and diversity/evenness values
data <- data.frame(Sample_ID = rownames(metadata), metadata[, c("Sex","Age","BMI","HTN","Diabetes","Respiratory_disease","Heart_disease","Renal_Disease","Liver_Disease","Obesity","Malignancy","Immunosuppressive_Disease","Neurological_disease","Metabolic_Disease","Cardiovascular_Disease","comorbidities_total","WHO_severity")], richness, diversity, evenness)
data$WHO_severity <- factor(data$WHO_severity, levels = c("mild", "moderate", "critical_severe"))

# Plot richness
ggplot(data, aes(x = WHO_severity, y = richness, fill = WHO_severity)) +
  geom_boxplot() +
  xlab("Covid Severity") +
  ylab("Richness") +
  ggtitle("Alpha Diversity: Richness")

# Plot Shannon diversity index
ggplot(data, aes(x = WHO_severity, y = diversity, fill = WHO_severity)) +
  geom_boxplot() +
  xlab("Covid Severity") +
  ylab("Shannon Diversity Index") +
  ggtitle("Alpha Diversity: Shannon Index")

# Plot Pielou's evenness index
ggplot(data, aes(x = WHO_severity, y = evenness, fill = WHO_severity)) +
  geom_boxplot() +
  xlab("Covid Severity") +
  ylab("Pielou's Evenness Index") +
  ggtitle("Evenness")


#After visualizing the boxplots, we can perform statistical tests to determine if there is a significant difference in the alpha diversity and evenness measures between the Covid severity groups.
#Since the data may not follow a normal distribution, we can see outliers in box plot, we can use non-parametric tests such as the Kruskal-Wallis test to compare the distributions of alpha diversity and evenness measures across the different Covid severity groups. We can perform these tests using the stats package in R. The code for performing the tests is as follows:
data <- subset(data, select = -richness)
metadata_diversity <- subset(data, select = -evenness)
metadata_evenness <- subset(data, select = -diversity)

# Kruskal-Wallis test for Shannon Index - overall alpha diversity 
kw_diversity <- kruskal.test(diversity ~ WHO_severity, data = metadata_diversity)

# Kruskal-Wallis test for Pielou Evenness
kw_evenness <- kruskal.test(evenness ~ WHO_severity, data = metadata_evenness)

# Print the test results
print(kw_diversity)
print(kw_evenness)

#Based on the Kruskal-Wallis rank sum test results, there is no significant variation in the evenness or overall alpha diversity (both richness and evenness quantified by the Shannon index) between the patients across the different severity groups. 
#The p-values for both tests are greater than 0.05, indicating that we cannot reject the null hypothesis that there is no significant difference between the groups.

#1b -----------------------------------------------------------------------------------

library(vegan)
library(ggplot2)
library(readxl)

#metadata <- read_excel("F:/IIIT-D/HMDS/Assignment1_Metadata.xlsx")
clr_data <-read_excel("F:/IIIT-D/HMDS/Assignment1_ClrTrans_Species.xlsx", sheet = 1)

#Next, we will create a distance matrix using the Bray-Curtis dissimilarity index:
dist_matrix <- vegdist(abundance, method = "bray")
#We can then perform the PCoA and plot the results:
pcoa <- cmdscale(dist_matrix, k = 3, eig = TRUE, add = TRUE)
pcoa_df <- data.frame(pcoa$points, metadata)
colnames(pcoa_df)[1:3] <- c("PC1", "PC2", "PC3")

pcoa_df$WHO_severity <- factor(pcoa_df$WHO_severity)

ggplot(pcoa_df, aes(x = PC1, y = PC2, color = WHO_severity)) +
  geom_point(size = 3) +
  theme_bw() +
  ggtitle("PCoA plot of gut microbiome beta-diversity") +
  xlab(paste0("PC1 (", round(pcoa$eig[1]/sum(pcoa$eig)*100, 2), "%)")) +
  ylab(paste0("PC2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100, 2), "%)")) +
  scale_color_manual(name = "WHO_severity", values = c("green", "red", "blue")) +
  guides(color = guide_legend(reverse = TRUE))

#This will produce a PCoA plot of the gut microbiome beta-diversity, with the points colored by the severity of Covid-19. We can use PERMANOVA (permutational multivariate analysis of variance) to test whether there is a significant difference in the microbiome composition across the different severity groups:
pcoa_adonis <- adonis2(dist_matrix ~ WHO_severity, data = metadata, permutations = 999)
pcoa_adonis
#This will give us the results of the PERMANOVA, including the R2 value, F statistic, and p-value. If the p-value is less than 0.05, we can conclude that there is a significant difference in the microbiome composition across the different severity groups.
#The p-value of 0.841 suggests that there is no significant difference 
#in the gut microbiome composition across the different severity groups 
#of COVID-19 in dataset.
