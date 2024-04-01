# Function: summary_analysis. 

# Description: This function calculates summary statistics from a given subset of data.
# It computes various metrics including:
# - Number of patients in the subset.
# - Sex information availability and the distribution of male and female patients.
# - Age information availability and the distribution of patients across different age groups.
# - Mutation type availability and the distribution of different mutation types.
# The function returns a data frame with the computed summary statistics.

phenotypes <- c("VI","MT","HL","HRT","LIV","REN","PUL","SHS","REP","TYD","MEND","ABFING","INT","SCO","NER","ALO")

summary_analysis <- function(subset_x) {
  
  # number of patients
  df_x = data.frame(count = nrow(subset_x))
  # sex
  
  df_x$sex = sum(is.na(subset_x$Sex)==FALSE)
  df_x$male = sum(subset_x$Sex=="M", na.rm = TRUE)
  df_x$female = sum(subset_x$Sex=="F", na.rm = TRUE)
  df_x$per_male = df_x$male/df_x$sex*100
  df_x$per_female = df_x$female/df_x$sex*100
  
  # age
  df_x$age = sum(is.na(subset_x$Age)==FALSE)
  df_x$age1 = sum(subset_x$Age=="0-9", na.rm = TRUE)
  df_x$age2 = sum(subset_x$Age=="10-19", na.rm = TRUE)
  df_x$age3 = sum(subset_x$Age=="20-29", na.rm = TRUE)
  df_x$age4 = sum(subset_x$Age=="30-39", na.rm = TRUE)
  df_x$age5 = sum(subset_x$Age=="40-49", na.rm = TRUE)
  df_x$age6 = sum(subset_x$Age=="50-59", na.rm = TRUE)
  df_x$per_age1 = df_x$age1/df_x$age*100
  df_x$per_age2 = df_x$age2/df_x$age*100
  df_x$per_age3 = df_x$age3/df_x$age*100
  df_x$per_age4 = df_x$age4/df_x$age*100
  df_x$per_age5 = df_x$age5/df_x$age*100
  df_x$per_age6 = df_x$age6/df_x$age*100
  
  
  # mutations
  df_x$mut = sum(is.na(subset_x$TM_1)==FALSE) + sum(is.na(subset_x$TM_2)==FALSE)
  df_x$cLOF = sum(subset_x$TM_1=="cLOF", na.rm = TRUE) + sum(subset_x$TM_2=="cLOF", na.rm = TRUE)
  df_x$MS = sum(subset_x$TM_1=="Ms", na.rm = TRUE) + sum(subset_x$TM_2=="Ms", na.rm = TRUE)
  df_x$intronic = sum(subset_x$TM_1=="intronic", na.rm = TRUE) + sum(subset_x$TM_2=="intronic", na.rm = TRUE)
  df_x$per_cLOF = df_x$cLOF/df_x$mut*100
  df_x$per_MS = df_x$MS/df_x$mut*100
  df_x$per_intronic = df_x$intronic/df_x$mut*100
  
  
  colnames(df_x) <- c("Number of patients",
                      "Sex information available",
                      "Number of male patients",
                      "Number of female patients",
                      "% of male",
                      "% of female",
                      "Age information available",
                      "Number of age 0-9 years",
                      "Number of age 10-19 years",
                      "Number of age 20-29 years",
                      "Number of age 30-39 years",
                      "Number of age 40-49 years",
                      "Number of age 50-59 years",
                      "% of age 0-9 years",
                      "% of age 10-19 years",
                      "% of age 20-29 years",
                      "% of age 30-39 years",
                      "% of age 40-49 years",
                      "% of age 50-59 years",
                      "Mutation type avaiable",
                      "Number of cLOF",
                      "Number of missense",
                      "Number of intronic",
                      "% of cLOF mutations",
                      "% of missense mutations",
                      "% of intronic mutations")
  
  df_x
  
}

# Function: phenotype_analysis
# 
# Description: This function calculates various statistics related to phenotype analysis based on the input data frame `df`. 
# It calculates counts and prevalence of phenotypes in the group. The function returns a data frame `df_x` containing 
# these statistics including the count of non-missing values, the count of "yes" values, the count of "no" values, 
# the count of missing values, the percentage of "yes" values, and confidence intervals for the prevalence of "yes" values.
# 
# Args:
#   df: A data frame containing information about phenotypes.
# 
# Returns:
#   df_x: A data frame with counts and prevalence of phenotypes.

phenotype_analysis <- function(df) {
  df_x <- data.frame( count = sum(is.na(df)==FALSE),
                      yes = sum(df, na.rm = TRUE),
                      no = sum(df=="0", na.rm = TRUE),
                      na = sum(is.na(df)))
  
  # prevalence of symptoms in group
  df_x$percent = round(df_x$yes/df_x$count*100,digits=2)
  df_x$upper = round(qbeta(0.975, df_x$yes + 1, df_x$no + 1)*100, digits = 2)
  df_x$lower = round(qbeta(0.025,df_x$yes + 1, df_x$no + 1)*100, digits = 2)
  df_x
}

# Function: prevalence_symptons_global
# 
# Description: This function calculates the prevalence of symptoms globally across different phenotypes based on 
# the input data frame `update8forR`. It calls the `phenotype_analysis` function for each phenotype to calculate
# the percentage of patients exhibiting each symptom. The function then compiles these percentages into a data frame 
# `prevalence` with each row corresponding to a phenotype and its corresponding prevalence percentage. The function 
# returns the `prevalence` data frame sorted in descending order of prevalence percentages.
# 
# Args:
#   update8forR: A list or data frame containing information about different phenotypes, each identified by a 
#                specific key.
# 
# Returns:
#   prevalence: A data frame containing the prevalence of symptoms across different phenotypes.
#

prevalence_symptons_global <- function(update8forR){
  prevalence <- data.frame( c( phenotype_analysis(update8forR$VI)$percent,
                               phenotype_analysis(update8forR$MT)$percent,
                               phenotype_analysis(update8forR$HL)$percent,
                               phenotype_analysis(update8forR$HRT)$percent,
                               phenotype_analysis(update8forR$LIV)$percent,
                               phenotype_analysis(update8forR$REN)$percent,
                               phenotype_analysis(update8forR$PUL)$percent,
                               phenotype_analysis(update8forR$SHS)$percent,
                               phenotype_analysis(update8forR$REP)$percent,
                               phenotype_analysis(update8forR$TYD)$percent,
                               phenotype_analysis(update8forR$MEND)$percent,
                               phenotype_analysis(update8forR$ABFING)$percent,
                               phenotype_analysis(update8forR$INT)$percent,
                               phenotype_analysis(update8forR$SCO)$percent,
                               phenotype_analysis(update8forR$NER)$percent,
                               phenotype_analysis(update8forR$ALO)$percent)
  )
  colnames(prevalence)<- "Percent"
  rownames(prevalence)<- phenotypes
  prevalence <- prevalence %>% arrange(desc(prevalence$Percent))
  prevalence
}

# Function: prevalence_symptons_top
# 
# Description: This function calculates the prevalence of symptoms for the top five phenotypes based on the input 
# data frame `update8forR`. It calls the `phenotype_analysis` function for each phenotype to calculate the percentage
# of patients exhibiting each symptom. The function then compiles these percentages into a data frame `prevalence` 
# with each row corresponding to a phenotype and its corresponding prevalence percentage. The function returns the 
# `prevalence` data frame sorted in descending order of prevalence percentages for the top five phenotypes.
# 
# Args:
#   update8forR: A list or data frame containing information about different phenotypes, each identified by a 
#                specific key.
# 
# Returns:
#   prevalence: A data frame containing the prevalence of symptoms for the top five phenotypes.
#

prevalence_symptons_top <- function(update8forR){
  prevalence <- data.frame( c( phenotype_analysis(update8forR$VI)$percent,
                               phenotype_analysis(update8forR$MT)$percent,
                               phenotype_analysis(update8forR$HL)$percent,
                               phenotype_analysis(update8forR$HRT)$percent,
                               phenotype_analysis(update8forR$LIV)$percent
  ))
  colnames(prevalence)<- "Percent"
  rownames(prevalence)<- c("VI","MT","HL","HRT","LIV")
  prevalence <- prevalence %>% arrange(desc(prevalence$Percent))
  prevalence
}

# Function: prevalence_statistics_chisq_global
# 
# Description: This function performs a chi-square test of significance to determine if there are significant 
# differences in the prevalence of symptoms between them. It constructs a contingency table 
# `prevalence_table` containing the counts of patients with and without each kind of symptoms. Then, it 
# conducts a chi-square test using the contingency table to assess the significance of differences in symptom 
# prevalence. The function returns the result of the chi-square test and the associated p-value.
# 
# Args:
#   update8forR: A list or data frame containing information about different phenotypes, each identified by a 
#                specific key.
# 
# Returns:
#   result: A list containing the result of the chi-square test and the associated p-value.
# 

prevalence_statistics_chisq_global <- function(update8forR){
  ## test of significance
  prevalence_table = matrix(c(
    phenotype_analysis(update8forR$VI)$no,
    phenotype_analysis(update8forR$MT)$no,
    phenotype_analysis(update8forR$HL)$no,
    phenotype_analysis(update8forR$HRT)$no,
    phenotype_analysis(update8forR$LIV)$no,
    phenotype_analysis(update8forR$REN)$no,
    phenotype_analysis(update8forR$PUL)$no,
    phenotype_analysis(update8forR$SHS)$no,
    phenotype_analysis(update8forR$REP)$no,
    phenotype_analysis(update8forR$TYD)$no,
    phenotype_analysis(update8forR$MEND)$no,
    phenotype_analysis(update8forR$ABFING)$no,
    phenotype_analysis(update8forR$INT)$no,
    phenotype_analysis(update8forR$SCO)$no,
    phenotype_analysis(update8forR$NER)$no,
    phenotype_analysis(update8forR$ALO)$no,
    phenotype_analysis(update8forR$VI)$yes,
    phenotype_analysis(update8forR$MT)$yes,
    phenotype_analysis(update8forR$HL)$yes,
    phenotype_analysis(update8forR$HRT)$yes,
    phenotype_analysis(update8forR$LIV)$yes,
    phenotype_analysis(update8forR$REN)$yes,
    phenotype_analysis(update8forR$PUL)$yes,
    phenotype_analysis(update8forR$SHS)$yes,
    phenotype_analysis(update8forR$REP)$yes,
    phenotype_analysis(update8forR$TYD)$yes,
    phenotype_analysis(update8forR$MEND)$yes,
    phenotype_analysis(update8forR$ABFING)$yes,
    phenotype_analysis(update8forR$INT)$yes,
    phenotype_analysis(update8forR$SCO)$yes,
    phenotype_analysis(update8forR$NER)$yes,
    phenotype_analysis(update8forR$ALO)$yes
  ), ncol = 16, byrow = TRUE)
  
  chisq.test(prevalence_table)
  chisq.test(prevalence_table)$p.value
}

# Function: prevalence_statistics_chisq_top
# Description: This function performs a chi-square test of significance for the top 5 prevalent symptoms.
# Input:
#   - update8forR: A dataframe containing binary data indicating the presence (1) or absence (0) of symptoms for each patient.
# Output:
#   - p-value: A numeric value representing the result of the chi-square test.
# Details: This function calculates the prevalence of the top 5 symptoms and performs a chi-square test of significance based on the presence or absence of these symptoms.

prevalence_statistics_chisq_top <- function(update8forR){
  ## test of significance
  prevalence_table = matrix(c(
    phenotype_analysis(update8forR$VI)$no,
    phenotype_analysis(update8forR$MT)$no,
    phenotype_analysis(update8forR$HL)$no,
    phenotype_analysis(update8forR$HRT)$no,
    phenotype_analysis(update8forR$LIV)$no,
    phenotype_analysis(update8forR$VI)$yes,
    phenotype_analysis(update8forR$MT)$yes,
    phenotype_analysis(update8forR$HL)$yes,
    phenotype_analysis(update8forR$HRT)$yes,
    phenotype_analysis(update8forR$LIV)$yes
  ), ncol = 5, byrow = TRUE)
  
  chisq.test(prevalence_table)
  chisq.test(prevalence_table)$p.value
}

# Function: plot_syndromic_score
# Description: This function generates a bar plot to visualize the distribution of syndromic scores.
# Input:
#   - df: A dataframe containing the syndromic scores.
#   - x_title: A character string specifying the title for the x-axis.
# Output:
#   - A bar plot visualizing the distribution of syndromic scores.
# Details: This function creates a bar plot using ggplot2 to display the distribution of syndromic scores. 
#          It adds a dashed vertical line at the mean of the syndromic scores for reference.

plot_syndromic_score <- function(df,x_title){
  ggplot(df, aes(x= SS))+
    geom_bar()+
    geom_vline(aes(xintercept = mean(SS)),col='black',size=1, linetype = "dashed")+
    scale_x_continuous(breaks=seq(0,1,0.2), limits = c(0,1.2))+
    scale_y_continuous(expand = c(0, 0))+
    labs(title = "Distribution Syndromic score",
         x = "Syndromic Score",
         y = x_title)+
    coord_flip()+
    theme_classic() +
    theme()
}

# Function: plot_syndromic_score_box
# Description: This function generates a box plot to visualize the distribution of syndromic scores across different groups.
# Input:
#   - df: A dataframe containing syndromic scores and grouping variable.
#   - x_var: A character string specifying the grouping variable.
#   - x_title: A character string specifying the title for the x-axis.
# Output:
#   - A box plot visualizing the distribution of syndromic scores across different groups.
# Details: This function creates a box plot using ggplot2 to display the distribution of syndromic scores 
#          across different groups specified by the x_var variable. It also includes jittered points for better visualization.

plot_syndromic_score_box <- function(df, x_var, x_title = "All patients n=227"){
  ggplot(df, aes(x= x_var, y=SS))+
    geom_boxplot(aes(fill = x_var, alpha = 0.1), outlier.shape = NA)+
    geom_jitter(aes(colour = x_var), height = 0.05, width = 0.2)+
    scale_color_brewer(palette = "Set2")+
    scale_fill_brewer(palette = "Set2")+
    scale_y_continuous(breaks=seq(0,1.4,0.2), limits = c(0,1.4))+
    labs(title = "Syndromic score distribution",
         x = x_title,
         y = "Syndromic score")+
    theme_bw() +
    theme(legend.position = "none")
}

# Function: plot_prevalence
# Description: This function generates a bar plot to visualize the prevalence of ALMS symptoms based on the specified mode of analysis.
# Input:
#   - update8forR: A list containing data frames of ALMS symptoms.
#   - mode_analysis: A character string specifying the mode of analysis ("global" or "top").
# Output:
#   - A bar plot visualizing the prevalence of ALMS symptoms.
# Details: This function first calculates the prevalence of ALMS symptoms and performs a test of significance based on the specified mode of analysis. 
#          It then generates a bar plot displaying the prevalence of symptoms along with the p-value from the chi-squared test.

plot_prevalence <- function(update8forR, mode_analysis){
  
  if(mode_analysis=="global"){
    ## prevalence of phenotypes 
    prevalence <- prevalence_symptons_global(update8forR)
    
    ## test of significance
    chisq_result  <- prevalence_statistics_chisq_global(update8forR)
    
    ## Plot results
    barplot(prevalence$Percent, ylim = c(0,100), ylab = "Prevalence of symptoms (%)", las = 1, names.arg = rownames(prevalence), main = "Prevalence of ALMS symptoms",
            xlab = paste("chisq.test p-value = ", format(chisq_result,digits = 4)))
    
    abline(h=15, lty = "longdash")
    
  }
  if(mode_analysis=="top"){
    ## prevalence of phenotypes 
    prevalence <- prevalence_symptons_top(update8forR)
    
    ## test of significance
    chisq_result  <- prevalence_statistics_chisq_top(update8forR)
    
    ## Plot results
    barplot(prevalence$Percent, ylim = c(0,110), ylab = "Prevalence of symptoms (%)", las = 1, names.arg = rownames(prevalence), main = "Prevalence of ALMS symptoms",
            xlab = paste("chisq.test p-value = ", format(chisq_result,digits = 4)))
    
    prevalence_top <- recordPlot()
  }
  
}

save_plot <- function(plot,name, w, h){
  pdf(file=paste0("./results/figures/",name), width=w, height=h)
  plot
}
