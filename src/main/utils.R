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
    barplot(prevalence$Percent, ylim = c(0,110), ylab = "Prevalence of symptoms (%)", las = 1, names.arg = rownames(prevalence), main = "prevalence of ALMS symptoms",
            xlab = paste("chisq.test p-value = ", format(chisq_result,digits = 4)))
    
    prevalence_top <- recordPlot()
  }
  
}

save_plot <- function(plot,name, w, h){
  pdf(file=paste0("./output/figures/",name), width=w, height=h)
  plot
}
