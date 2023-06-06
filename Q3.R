## Installing Packages
#install.packages("vegan") 
#install.packages("ecodist") 
#install.packages("ggplot2")
#install.packages("dplyr")                             # Install dplyr package
#install.packages("readxl")

library("readxl")
library(vegan)
library(ecodist)
library("dplyr") 
library("ggplot2")
setwd("C:\\Users\\anany\\Desktop\\semester8\\hmds\\Ass1")

##Question 3##
meth<-"bray"
data1 <- read_excel("Assignment1_ClrTrans_Species.xlsx")
data1 <- data1[order(data1$X),]

data2<-read_excel("Assignment1_Metadata.xlsx")
data2 <- data2[order(data2$X),]

#Converting to numeric
data2$WHO_severity[which(data2$WHO_severity=="mild")] <- 1
data2$WHO_severity[which(data2$WHO_severity=="moderate")]<-2
data2$WHO_severity[which(data2$WHO_severity=="critical_severe")] <- 3

data2$Sex[which(data2$Sex=="M")]<-1
data2$Sex[which(data2$Sex=="F")] <- 2

data2[data2=="y"] <-"1"
data2[data2=="n"] <-"0"

data2_new<-data2[,-1] ###Metadata
data1<-data1[,-1] ###data

data2_new <- mutate_all(data2_new, function(x) as.numeric(as.character(x)))
newdf <- mutate_all(data1, function(x) as.numeric(as.character(x)))

data2_sex.div<-adonis2(newdf ~ Sex, data = data2_new, permutations = 1000, method=meth)
df_data2_sex <- data.frame(data2_sex.div)

data2_age.div<-adonis2(newdf ~ Age, data = data2_new, permutations = 1000, method=meth)
df_data2_age<-data.frame(data2_age.div)

#data2_age.div['Pr(>F)']
data2_BMI.div<-adonis2(newdf ~ BMI, data = data2_new, permutations = 1000, method=meth)
df_data2_BMI<-data.frame(data2_BMI.div)

data2_WHO_severity.div<-adonis2(newdf ~ WHO_severity, data = data2_new, permutations = 1000, method=meth)
df_data2_WHO_severity<-data.frame(data2_WHO_severity.div)

data2_comorbidities_total.div<-adonis2(newdf ~ comorbidities_total, data = data2_new, permutations = 1000, method=meth)
df_data2_comorbidities_total<-data.frame(data2_comorbidities_total.div)

data2_HTN.div<-adonis2(newdf ~ HTN, data = data2_new, permutations = 1000, method=meth)
df_data2_HTN<-data.frame(data2_HTN.div)

data2_Diabetes.div<-adonis2(newdf ~ Diabetes, data = data2_new, permutations = 1000, method=meth)
df_data2_Diabetes<-data.frame(data2_Diabetes.div)

data2_Respiratory_disease.div<-adonis2(newdf ~ Respiratory_disease, data = data2_new, permutations = 1000, method=meth)
df_data2_Respiratory_disease<-data.frame(data2_Respiratory_disease.div)

data2_Heart_disease.div<-adonis2(newdf ~ Heart_disease, data = data2_new, permutations = 1000, method=meth)
df_data2_Heart_disease<-data.frame(data2_Heart_disease.div)

data2_Renal_Disease.div<-adonis2(newdf ~ Renal_Disease, data = data2_new, permutations = 1000, method=meth)
df_data2_Renal_disease<-data.frame(data2_Renal_Disease.div)

data2_Liver_Disease.div<-adonis2(newdf ~ Liver_Disease, data = data2_new, permutations = 1000, method=meth)
df_data2_Liver_disease<-data.frame(data2_Liver_Disease.div)

data2_Obesity.div<-adonis2(newdf ~ Obesity, data = data2_new, permutations = 1000, method=meth)
df_data2_Obesity<-data.frame(data2_Obesity.div)

data2_Malignancy.div<-adonis2(newdf ~ Malignancy, data = data2_new, permutations = 1000, method=meth)
df_data2_Malignancy<-data.frame(data2_Malignancy.div)

data2_Immunosuppressive_Disease.div<-adonis2(newdf ~ Immunosuppressive_Disease, data = data2_new, permutations = 1000, method=meth)
df_data2_Immunosuppressive_disease<-data.frame(data2_Immunosuppressive_Disease.div)


data2_Neurological_disease.div<-adonis2(newdf ~ Neurological_disease, data = data2_new, permutations = 1000, method=meth)
df_data2_Neurological_disease<-data.frame(data2_Neurological_disease.div)

data2_Metabolic_Disease.div<-adonis2(newdf ~ Metabolic_Disease, data = data2_new, permutations = 1000, method=meth)
df_data2_Metabolic_disease<-data.frame(data2_Metabolic_Disease.div)

data2_Cardiovascular_Disease.div<-adonis2(newdf ~ Cardiovascular_Disease, data = data2_new, permutations = 1000, method=meth)
df_data2_Cardiovascular_disease<-data.frame(data2_Cardiovascular_Disease.div)

df_all<-rbind(df_data2_sex[1,], df_data2_age[1,],df_data2_BMI[1,],df_data2_WHO_severity[1,],df_data2_comorbidities_total[1,],df_data2_HTN[1,],df_data2_Diabetes[1,],df_data2_Respiratory_disease[1,],df_data2_Heart_disease[1,],df_data2_Renal_disease[1,],df_data2_Liver_disease[1,],df_data2_Obesity[1,],df_data2_Malignancy[1,],df_data2_Immunosuppressive_disease[1,],df_data2_Neurological_disease[1,],df_data2_Metabolic_disease[1,],df_data2_Cardiovascular_disease[1,])
df_all_desc<-arrange(df_all,desc(df_all$Pr..F.))
df_all_desc
write.csv(df_all_desc,"C:\\Users\\anany\\Desktop\\semester8\\hmds\\Ass1\\Metadata_MicrobiomeCompostion.csv") 

####Permanova is used to find the association between Microbiome Composition and Metadata.
## WHO_severity , Respiratory_disease , HTN and Malignancy were found to have significant associations ( <= 0.05)

