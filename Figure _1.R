
#################################### Complex Metadata heatmap of pre-ICI cohort ############################################################
# load librarries
library(openxlsx)
library(ComplexHeatmap)

# load data
data <- read.xlsx("/Users......./Supp_Table_1.xlsx")

# clean data frame 
df <- data # rename data frame

# remove CIBERSORTx columns 
df2 <- df[, -9: -30] # remove addtional columns that will not appear in heatmap

# add numeric columns of 0 and 1 
df2['column1'] = 1
df2['column2'] = 0

# Assuming your dataframe is named 'df'

# Select annotation columns
annotation_df <- df2[, c("Cohort", "Response", "Detailed.Response", "Therapy", "Overall_survival", "Alive", "Imm.Subtype")]

# Create the Heatmap object
Heatmap(
  t(df2[, c("column1", "column2")]),
  name = "Value",
  bottom_annotation = HeatmapAnnotation(df = annotation_df, which = "column"),
  cluster_rows = FALSE
)

annotation_df.sorted <- annotation_df[order(annotation_df$Cohort,
                                            annotation_df$Response,
                                            annotation_df$Detailed.Response),]
# Create the Heatmap object
Heatmap(
  t(df2[, c("column1", "column2")]),
  name = "Value",
  bottom_annotation = HeatmapAnnotation(df = annotation_df.sorted, which = "column"),
  cluster_rows = FALSE
)

########################################### Imm Subtypes and Response plots #################################################################
library(openxlsx)
#load data
data <- read.xlsx("/Users......./Supp_Table_1.xlsx", sheet = 2)

#data filter
#remove survival as continuous variation

#data <- data[, -6]


#load libraray
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggplot2)
library(readxl)

library(RColorBrewer)
par(mar=c(3,4,2,2))
display.brewer.all()

nb.cols <- 22
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)


# convert survival numeric data to character
data$Overall_survival <- as.character(data$Overall_survival)

# melt data longer format
data_melt <- melt(data)
data_melt <- as.data.frame(data_melt)

#check value coloumn  

############################# Imm Subtype and Response ######################
p <- ggplot(data_melt, aes(Imm.Subtype, fill= Response)) +
  geom_bar( position='fill')+
  ylab("fraction") + xlab("Immune Subtype")

p

countSubtype <- 
  
  # cohort as facet 
  plot <- p + facet_grid( ~ Cohort, scales = "free", space = "free_x")

plot 

################################### End ########################################

############################# Patient and cohort specific######################
p <-  ggplot(data_melt, aes(x= Mixture, y=value, fill=variable)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = mycolors) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("fraction") + xlab("")
p

# cohort as facet 
plot <- p + facet_grid( ~ Cohort, scales = "free", space = "free_x")

plot 

################################### End ########################################

#################### Patient and Response specific##############################
p <-  ggplot(data_melt, aes(x= Mixture, y=value, fill=variable)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = mycolors) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("fraction") + xlab("")
p


# cohort as facet 
plot <- p + facet_grid( ~ Response, scales = "free", space = "free_x")  


plot

###################################End##########################################

### Avg fraction cohortwise
p <-  ggplot(data_melt, aes(x= Cohort, y=value, fill=variable)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = mycolors) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("fraction") + xlab("")
p

# cohort as facet 
plot <- p + facet_grid( ~ Response, scales = "free", space = "free_x")

plot

############################## STAT ############################################
library(rstatix)
library(ggpubr)


stat.test <- data_melt %>%
  group_by(variable) %>%
  anova_test(value ~ Cohort) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

write.xlsx(stat.test, "LM22_Cohort_anova_test.xlsx")


#################################### End #######################################

# Clinical data plots

# Butterfly plot of Response and Cohorts

# load Cohort, Response, Therapy Count data

count <- read.xlsx("/Users......./Supp_Table_1.xlsx.xlsx")
# Load libraries
library(tidyverse)
install.packages("conmat")
library(conmat)

# Adding negative value when Response category is "NR" 
count_pyramid <- count %>% 
  mutate (Value = case_when(
    Response == "NR" ~ -Value,
    TRUE ~ Value))
count_pyramid

# Check range of X-axis
pop_range <- range(count_pyramid$Value)
pop_range

# Making negative values to positive
pop_range_seq <- seq(pop_range[1], pop_range[2], by = 10)
pop_range_seq
pop_range_breaks <- pretty(pop_range, n = 7)

#ggplot
plot <- ggplot(count_pyramid, aes(x = Value, y = Cohort, fill = Response)) +
  geom_col() +
  facet_wrap(vars(Therapy), ncol = 2, scales = "free_x")
scale_x_continuous(breaks  = pop_range_breaks,
                   labels = abs(pop_range_breaks))

plot

################################################################################

#################### Patient and Response and Subtype specific##############################

p <-  ggplot(data_melt, aes(x= Response, y=value, fill=variable)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = mycolors) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("fraction") + xlab("")

p

# cohort as facet 
plot <- p + facet_grid( ~ Imm.Subtype, scales = "free", space = "free_x")  

plot

#STAT
#Subtype 1
Sub1 <- data_melt[data_melt$Imm.Subtype == "Subtype1", ]
Sub1

#Subtype 2
Sub2 <- data_melt[data_melt$Imm.Subtype == "Subtype2", ]
Sub2

#Subtype 3
Sub3 <- data_melt[data_melt$Imm.Subtype == "Subtype3", ]
Sub3

#Subtype 4
Sub4 <- data_melt[data_melt$Imm.Subtype == "Subtype4", ]
Sub4


#Subtype 6
Sub6 <- data_melt[data_melt$Imm.Subtype == "Subtype6", ]
Sub6

stat.test <- Sub4 %>%
  group_by(variable) %>%
  wilcox_test(value ~ Response) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

write.xlsx(stat.test, "Sub4_Response_M2_M1_wilcox_test.xlsx")

###################################End##########################################

#M2 Violin plot

col_sub <- c("#F55050", "#FFEA20", "#ADE792", "#00FFF6", "#C539B4")

v <- ggplot(data, aes(x = Response, y = Macrophages.M2, fill = Imm.Subtype)) +
  geom_jitter(position=position_jitter(0.2), size = 1, alpha = 0.3, colour = "black") +
  scale_fill_manual(values=c("#F55050", "#FFEA20", "#ADE792", "#00FFF6", "#C539B4"))+
  facet_wrap( Imm.Subtype ~ . ) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black")+ 
  theme_bw()

v

################################################################################
#M2/M1 ratio plot

col_sub <- c("#F55050", "#FFEA20", "#ADE792", "#00FFF6", "#C539B4")

v <- ggplot(data, aes(x = Response, y = M2_M1_Ratio, fill = Imm.Subtype)) +
  geom_jitter(position=position_jitter(0.1), size = 1.5, alpha = 0.2) +
  scale_fill_manual(values=c("#F55050", "#FFEA20", "#ADE792", "#00FFF6", "#C539B4"))+
  facet_wrap( Imm.Subtype ~ . ) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black")+ 
  theme_bw()

v


