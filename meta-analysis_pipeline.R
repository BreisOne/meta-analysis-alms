#---
#title: "Analysis of genotype-phenotype associations in ALMS"
#author: "Brais Bea Mascato"
#date: '`r format(Sys.time(), "%d %m, %Y")`'
#output:
#---
  
libraries <- c('InteractionSet','trackViewer','org.Hs.eg.db','TxDb.Hsapiens.UCSC.hg19.knownGene','kableExtra','gridGraphics','ggpubr','cowplot','rstatix','readxl','dendsort','apeglm','VennDiagram', 'RColorBrewer', 'pheatmap', 'tidyverse','scales','ggrepel')
lapply(libraries,library, character.only = TRUE)

setwd("C:/Users/Brise/OneDrive - Universidade de Vigo/Tesis - Universidad de Vigo/Meta analisis ALMS1")

Curated_ALMS1_DB <- read_xlsx("Curated ALMS1 DB.xlsx",
                               sheet = "Genotype-Phenotype")

ALMS1_cohort <- Curated_ALMS1_DB[complete.cases(Curated_ALMS1_DB[,c(5:18,20:36)]),c(1:36)]
colnames(ALMS1_cohort)[6:15] <- c("PV1_DNA", "EA1", "PV1_PROTEIN","TM_1","PV2_DNA", "EA2", "PV2_PROTEIN","TM_2","GS","LT")
colnames(ALMS1_cohort)[36] <- "SS"

ALMS1_cohort <- ALMS1_cohort %>% mutate(group = as.factor(ifelse(LT < 9, 1,ifelse(LT >= 9 & LT < 14,2,3))))

phenotypes <- c("VI","MT","HRT","HL","LIV","REN","RES","SHS","REP","TYD","MEND","ABFING","INT","SCO","NER","ALO")


####Functions####
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
  
  # penetrance of symptoms in group
  df_x$percent = round(df_x$yes/df_x$count*100,digits=2)
  df_x$upper = round(qbeta(0.975, df_x$yes + 1, df_x$no + 1)*100, digits = 2)
  df_x$lower = round(qbeta(0.025,df_x$yes + 1, df_x$no + 1)*100, digits = 2)
  df_x
}

penetrance_symptons_global <- function(update8forR){
  penetrance <- data.frame( c( phenotype_analysis(update8forR$VI)$percent,
                               phenotype_analysis(update8forR$MT)$percent,
                               phenotype_analysis(update8forR$HRT)$percent,
                               phenotype_analysis(update8forR$HL)$percent,
                               phenotype_analysis(update8forR$LIV)$percent,
                               phenotype_analysis(update8forR$REN)$percent,
                               phenotype_analysis(update8forR$RES)$percent,
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
  colnames(penetrance)<- "Percent"
  rownames(penetrance)<- phenotypes
  penetrance <- penetrance %>% arrange(desc(penetrance$Percent))
  penetrance
}

penetrance_symptons_top <- function(update8forR){
  penetrance <- data.frame( c( phenotype_analysis(update8forR$VI)$percent,
                               phenotype_analysis(update8forR$MT)$percent,
                               phenotype_analysis(update8forR$HRT)$percent,
                               phenotype_analysis(update8forR$HL)$percent,
                               phenotype_analysis(update8forR$LIV)$percent
  ))
  colnames(penetrance)<- "Percent"
  rownames(penetrance)<- c("VI","MT","HRT","HL","LIV")
  penetrance <- penetrance %>% arrange(desc(penetrance$Percent))
  penetrance
}

penetrance_statistics_chisq_global <- function(update8forR){
        ## test of significance
        penetrance_table = matrix(c(
          phenotype_analysis(update8forR$VI)$no,
          phenotype_analysis(update8forR$MT)$no,
          phenotype_analysis(update8forR$HRT)$no,
          phenotype_analysis(update8forR$HL)$no,
          phenotype_analysis(update8forR$LIV)$no,
          phenotype_analysis(update8forR$REN)$no,
          phenotype_analysis(update8forR$RES)$no,
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
          phenotype_analysis(update8forR$HRT)$yes,
          phenotype_analysis(update8forR$HL)$yes,
          phenotype_analysis(update8forR$LIV)$yes,
          phenotype_analysis(update8forR$REN)$yes,
          phenotype_analysis(update8forR$RES)$yes,
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
        
        chisq.test(penetrance_table)
        chisq.test(penetrance_table)$p.value
}

penetrance_statistics_chisq_top <- function(update8forR){
  ## test of significance
  penetrance_table = matrix(c(
    phenotype_analysis(update8forR$VI)$no,
    phenotype_analysis(update8forR$MT)$no,
    phenotype_analysis(update8forR$HRT)$no,
    phenotype_analysis(update8forR$HL)$no,
    phenotype_analysis(update8forR$LIV)$no,
    phenotype_analysis(update8forR$VI)$yes,
    phenotype_analysis(update8forR$MT)$yes,
    phenotype_analysis(update8forR$HRT)$yes,
    phenotype_analysis(update8forR$HL)$yes,
    phenotype_analysis(update8forR$LIV)$yes
  ), ncol = 5, byrow = TRUE)
  
  chisq.test(penetrance_table)
  chisq.test(penetrance_table)$p.value
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
    labs(title = "Distribution Syndromic score",
         x = x_title,
         y = "Syndromic score")+
    theme_bw() +
    theme(legend.position = "none")
}

plot_penetrance <- function(update8forR, mode_analysis){
  
  if(mode_analysis=="global"){
    ## penetrance of phenotypes 
    penetrance <- penetrance_symptons_global(update8forR)
    
    ## test of significance
    chisq_result  <- penetrance_statistics_chisq_global(update8forR)
    
    ## Plot results
    barplot(penetrance$Percent, ylim = c(0,100), ylab = "Penetrance of symptoms (%)", las = 1, names.arg = rownames(penetrance), main = "Penetrance of ALMS symptoms",
            xlab = paste("chisq.test p-value = ", format(chisq_result,digits = 4)))
    
    abline(h=15, lty = "longdash")
    
  }
  if(mode_analysis=="top"){
    ## penetrance of phenotypes 
    penetrance <- penetrance_symptons_top(update8forR)
    
    ## test of significance
    chisq_result  <- penetrance_statistics_chisq_top(update8forR)
    
    ## Plot results
    barplot(penetrance$Percent, ylim = c(0,110), ylab = "Penetrance of symptoms (%)", las = 1, names.arg = rownames(penetrance), main = "Penetrance of ALMS symptoms",
                        xlab = paste("chisq.test p-value = ", format(chisq_result,digits = 4)))
    
    penetrance_top <- recordPlot()
  }

}

save_plot <- function(plot,name, w, h){
  pdf(file=paste0("./Figures/",name), width=w, height=h)
  plot
}

####Analysis####

###Summary###

Cohort_summary <- t(summary_analysis(ALMS1_cohort))
colnames(Cohort_summary) <- "Cohort of study"


kable(Cohort_summary, format = "html") %>%
  kable_styling(full_width = F, font_size = 9,bootstrap_options = c("striped", "hover", "condensed", "responsive"))


##Mutations in the cohort
ALMS1_cohort <- ALMS1_cohort %>% mutate( concat_1 = paste0(ALMS1_cohort$PV1_PROTEIN,ALMS1_cohort$PV1_DNA), concat_2= paste0(ALMS1_cohort$PV2_PROTEIN, ALMS1_cohort$PV2_DNA))

df_mutations <- data.frame(table(c(str_squish(ALMS1_cohort$PV1_PROTEIN),str_squish(ALMS1_cohort$PV2_PROTEIN))))
colnames(df_mutations) <- c("mutation_protein","count")
df_mutations <- df_mutations[df_mutations$mutation_protein!="splice site",]

# df_mutations_concat<- ALMS1_cohort[ALMS1_cohort$PV1_PROTEIN!="splice site"&ALMS1_cohort$PV2_PROTEIN!="splice site",]
# 
# df_mutations_concat<- data.frame(table(c(str_squish(ALMS1_cohort$concat_1),str_squish(ALMS1_cohort$concat_2))))

df_mutations_cDNA<- ALMS1_cohort[ALMS1_cohort$PV1_PROTEIN!="splice site"&ALMS1_cohort$PV2_PROTEIN!="splice site",]

df_mutations_cDNA<- data.frame(table(c(str_squish(ALMS1_cohort$PV1_DNA),str_squish(ALMS1_cohort$PV2_DNA))))
colnames(df_mutations_cDNA) <- c("mutation_cDNA","count")

mutations_cohort <- ggplot(df_mutations, aes(x = mutation_protein , y = count) )+
                        geom_bar(
                          # fill= "#8EA0CB",
                          stat="identity")+
                        labs(title = "Count of mutations in ALMS1 cohort",
                             x = "Mutations",
                             y = "Count")+
                        # scale_fill_grey()+
                        theme_classic2() +
                        coord_flip()+
                        theme(legend.position = "none")
mutations_cohort

save_plot(mutations_cohort,"_FigS1_mutations_cohort.pdf",5,30)
dev.off()

mutations_cohort_top <- ggplot(df_mutations[which(df_mutations$count>2),], aes(x = mutation_protein , y = count) )+
                            geom_bar(stat="identity")+
                            labs(title = "Count of mutations in ALMS1 cohort",
                                 x = "Mutations",
                                 y = "Count")+
                            scale_y_continuous(breaks=seq(0,40,5), limits = c(0,40))+
                            # scale_fill_grey()+
                            theme_classic2() +
                            theme(legend.position = "none",
                                  axis.text.x = element_text(angle = 45, hjust=1),
                                  plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
mutations_cohort_top

save_plot(mutations_cohort_top,"_Fig1_mutations_cohort_top.pdf",9,5)
dev.off()

## Penetrance of mutations in ALMS1 exons genes

df <- data.frame(table(c(ALMS1_cohort$EA1,ALMS1_cohort$EA2)))
colnames(df) <- c("exon","count")
df <- df[df$exon!="intron 17",]
df$exon <- factor(df$exon, levels = c(1,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20))
df$percentage = round(df$count/length(c(ALMS1_cohort$EA1,ALMS1_cohort$EA2))*100,digits = 2)

Penetrance_alleles <- ggplot(df, aes(x = exon, y = count) )+
                      geom_bar(stat="identity")+
                      labs(title = "Penetrance of mutated alelles in ALMS1 exons",
                           x = "Exon",
                           y = "Alleles")+
                      # scale_fill_grey()+
                      theme_classic2() +
                      theme(legend.position = "none")
Penetrance_alleles

save_plot(Penetrance_alleles,"_Fig1_Penetrance_alleles.pdf",7,5)
dev.off()

#Groups by Longest transcript
groups_LT <-    ggplot(ALMS1_cohort, aes(x=LT))+
                geom_bar() +
                labs(title = "Groups by longest transcript",
                     x = "Exon",
                     y = "Number of patients")+
                scale_x_continuous(breaks=seq(5,20,1))+
                theme_classic2() +
                theme(legend.position = "none")

groups_LT

save_plot(groups_LT,"_Fig1_groups_LT.pdf",7,5)
dev.off()

#Agregated groups of Longest transcripts with age stratification

agregated_groups_LT <-  ggplot(ALMS1_cohort, aes(x= group, fill = Age))+
                            geom_bar() +
                            labs(title = "Groups by clusterised Longest trancript",
                                 x = "Group",
                                 y = "Number of patients")+
                            scale_x_discrete(labels= c("1"="E8","2"="E10","3"="E16"))+
                            scale_fill_grey()+
                            theme_classic2()
  
agregated_groups_LT

save_plot(agregated_groups_LT,"_Fig1_agregated_groups_LT.pdf",7,5)
dev.off()

#Groups by Sex bar_plot
groups_by_sex <-ggplot(ALMS1_cohort, aes(x= Sex, fill=Sex))+
                geom_bar() +
                labs(title = "Cohort Patients",
                     y = "Number of patients",
                     x = "Sex")+
                scale_fill_grey()+
                theme_classic2() +
                theme(legend.position = "none")

groups_by_sex

save_plot(groups_by_sex,"_Fig1_groups_by_sex.pdf",7,5)
dev.off()

#Penetrance of symptoms

pdf(file = paste0("./Figures/", "_FigS1_penetrance_global.pdf"), width = 12, height = 5)
plot_penetrance(ALMS1_cohort, "global")
dev.off()


penetrance_top <- plot_penetrance(ALMS1_cohort, "top")
save_plot(penetrance_top,"_Fig2_penetrance_top.pdf",5,5)
dev.off()

#Distribution syndromic score complete cohort

allpatients_ss <- plot_syndromic_score(ALMS1_cohort, "All patients \n n=227")
allpatients_ss

save_plot(allpatients_ss,"_Fig2_allpatients_ss.pdf",5,7)
dev.off()

#Distribution syndromic score by sex

sexgroups_ss_box_plot <- plot_syndromic_score_box(ALMS1_cohort, ALMS1_cohort$Sex)+
                          stat_compare_means(method = "wilcox.test", label.y = 1.3 )

sexgroups_ss_box_plot

save_plot(sexgroups_ss_box_plot,"_Fig2_sexgroups_ss_boxplot.pdf",6,6)
dev.off()

#Distribution syndromic score by ages

agegroups_ss_box_plot <- plot_syndromic_score_box(ALMS1_cohort, ALMS1_cohort$Age2)+
                          stat_compare_means(method = "kruskal.test", label.y = 1.3)

agegroups_ss_box_plot

save_plot(agegroups_ss_box_plot,"_Fig2_agegroups_ss_box_plot.pdf",6,6)
dev.off()

#Distribution syndromic score in subgroups (E8, E10, E16)

#Distribution syndromic score in subgroups
stat.test <- compare_means(SS ~ group,  data = ALMS1_cohort, p.adjust.method ="BH")


subgroups_ss_box_plot <- plot_syndromic_score_box(ALMS1_cohort, ALMS1_cohort$group)+
                          scale_x_discrete(labels=c("1" = "E8", "2" = "E10", "3" = "E16"))+
                                stat_pvalue_manual(
                                  stat.test, 
                                  y.position = 1.1, step.increase = 0.1,
                                  label = "p.adj"
                                ) + # Add adj.p-value
                          stat_compare_means(method = "kruskal.test", label.y = 1.3) # Add global p-value
  
subgroups_ss_box_plot

save_plot(subgroups_ss_box_plot,"_Fig2_subgroups_ss_boxplot.pdf",6,6)
dev.off()

#Distribution syndromic score in subgroups by ages

stat.test <- ALMS1_cohort %>%
  group_by(Age2)%>%
  pairwise_wilcox_test(SS ~ group, p.adjust.method ="BH")

stat.test

subgroups_ss_ages_box_plot <- plot_syndromic_score_box(ALMS1_cohort, ALMS1_cohort$group)+
                                  scale_x_discrete(labels=c("1" = "E8", "2" = "E10", "3" = "E16"))+
                                  facet_grid(cols = vars(Age2))+
                                  stat_pvalue_manual(
                                    stat.test, 
                                    y.position = 1.1, 
                                    step.increase = 0.08,
                                    hide.ns = TRUE,
                                    label = "p.adj")

subgroups_ss_ages_box_plot

save_plot(subgroups_ss_ages_box_plot,"_Fig3_subgroups_ss_ages_boxplot.pdf",9,6)
dev.off()

#Penetrance of mains phenotypes features in subgroups (E8,E10,E16)

group_exon8 <- ALMS1_cohort[ALMS1_cohort$group==1,]
group_exon10 <- ALMS1_cohort[ALMS1_cohort$group==2,]
group_exon16 <- ALMS1_cohort[ALMS1_cohort$group==3,]

Phenotypes_groups <- data.frame( "Penetrance"=c(phenotype_analysis(group_exon8$VI)$percent,
                                                phenotype_analysis(group_exon8$MT)$percent,
                                                phenotype_analysis(group_exon8$HL)$percent,
                                                phenotype_analysis(group_exon8$HRT)$percent,
                                                phenotype_analysis(group_exon8$LIV)$percent,
                                                phenotype_analysis(group_exon8$REN)$percent,
                                                phenotype_analysis(group_exon8$MEND)$percent,
                                                phenotype_analysis(group_exon8$RES)$percent,
                                                phenotype_analysis(group_exon8$REP)$percent,
                                               phenotype_analysis(group_exon10$VI)$percent,
                                               phenotype_analysis(group_exon10$MT)$percent,
                                               phenotype_analysis(group_exon10$HL)$percent,
                                               phenotype_analysis(group_exon10$HRT)$percent,
                                               phenotype_analysis(group_exon10$LIV)$percent,
                                               phenotype_analysis(group_exon10$REN)$percent,
                                               phenotype_analysis(group_exon10$MEND)$percent,
                                               phenotype_analysis(group_exon10$RES)$percent,
                                               phenotype_analysis(group_exon10$REP)$percent,
                                               phenotype_analysis(group_exon16$VI)$percent,
                                               phenotype_analysis(group_exon16$MT)$percent,
                                               phenotype_analysis(group_exon16$HL)$percent,
                                               phenotype_analysis(group_exon16$HRT)$percent,
                                               phenotype_analysis(group_exon16$LIV)$percent,
                                               phenotype_analysis(group_exon16$REN)$percent,
                                               phenotype_analysis(group_exon16$MEND)$percent,
                                               phenotype_analysis(group_exon16$RES)$percent,
                                               phenotype_analysis(group_exon16$REP)$percent),
                                    
                                  "Symptom" =c(rep(c("VI","MT","HRT","HL","LIV","REN","MEND","RES","REP"),3)),
                                    
                                   "group" =as.factor(c(rep(1,9),rep(2,9),rep(3,9)))
                                 ) 

Phenotypes_groups$Symptom <- factor(Phenotypes_groups$Symptom, levels=c("VI","MT","HRT","HL","LIV","REN","MEND","RES","REP"))


#Plot penetrance of mains phenotypes features in subgroups (E8,E10,E16)
  
# phenotypes_by_group <- ggplot(Phenotypes_groups, aes(x= group, y= Penetrance, fill=group),stat = "summary", fun.y = "mean")+
#                                 geom_bar(stat = "identity")+
#                                 facet_grid(cols = vars(Symptom))+
#                                 scale_x_discrete(labels=c("1" = "E8", "2" = "E10", "3" = "E16"))+
#                                 scale_fill_brewer(palette = "Set2")+
#                                 labs(title = "Penetrance by symptom in each group",
#                                            x="All patietns n = 227",
#                                            y="Penetrance of the symptom(%)")+
#                                 
#                                 theme_bw()+
#                                 theme(legend.position = "none")
# 
# phenotypes_by_group
# save_plot(phenotypes_by_group,"phenotypes_by_group.pdf",7,5)
# dev.off()

#Plot penetrance of mains phenotypes features in subgroups (E8,E10,E16) using beta-distributions

Phenotypes_groups_qvalues <- data.frame( "Penetrance"=c(unlist(phenotype_analysis(group_exon8$VI)[,5:7]),
                                                unlist(phenotype_analysis(group_exon8$MT)[,5:7]),
                                                unlist(phenotype_analysis(group_exon8$HL)[,5:7]),
                                                unlist(phenotype_analysis(group_exon8$HRT)[,5:7]),
                                                unlist(phenotype_analysis(group_exon8$LIV)[,5:7]),
                                                unlist(phenotype_analysis(group_exon8$REN)[,5:7]),
                                                unlist(phenotype_analysis(group_exon8$MEND)[,5:7]),
                                                unlist(phenotype_analysis(group_exon8$RES)[,5:7]),
                                                unlist(phenotype_analysis(group_exon8$REP)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$VI)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$MT)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$HL)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$HRT)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$LIV)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$REN)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$MEND)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$RES)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$REP)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$VI)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$MT)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$HL)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$HRT)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$LIV)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$REN)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$MEND)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$RES)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$REP)[,5:7])),
                                 
                                 "Symptom" =c(rep(c(rep("VI",3),rep("MT",3),rep("HRT",3),rep("HL",3),rep("LIV",3),rep("REN",3),rep("MEND",3),rep("RES",3),rep("REP",3)),3)),
                                 
                                 "group" =as.factor(c(rep(1,27),rep(2,27),rep(3,27)))
) 

Phenotypes_groups_qvalues$Symptom <- factor(Phenotypes_groups_qvalues$Symptom, levels=c("VI","MT","HRT","HL","LIV","REN","MEND","RES","REP"))

df.summary <- Phenotypes_groups_qvalues %>%
                  group_by(group, Symptom) %>%
                  summarise(
                    sd = sd(Penetrance, na.rm = TRUE),
                    Penetrance = mean(Penetrance)
                  )
df.summary

for(i in c(20:26,28,30)){
  
  if(i==20){df <- data.frame()}
  
  df <- rbind(df,data.frame(list(Symptom=phenotypes[i-19], 
                                 .y.= rep("Penetrance",3), 
                                 pairwise_fisher_test(matrix( c(
                                   phenotype_analysis(group_exon8[i])$no,
                                   phenotype_analysis(group_exon10[i])$no,
                                   phenotype_analysis(group_exon16[i])$no,
                                   phenotype_analysis(group_exon8[i])$yes,
                                   phenotype_analysis(group_exon10[i])$yes,
                                   phenotype_analysis(group_exon16[i])$yes),
                                   byrow = TRUE, 
                                   ncol = 3,
                                   dimnames = list(c("No","Yes"), c("1","2","3"))),
                                   p.adjust.method = "BH"))))
  df                         
}

stat.test <- as_tibble(df)
stat.test$Symptom <- factor(stat.test$Symptom, levels=c("VI","MT","HRT","HL","LIV","REN","MEND","RES","REP"))

kable(stat.test, format = "html") %>%
  kable_styling(full_width = F, font_size = 9,bootstrap_options = c("striped", "hover", "condensed", "responsive"))

phenotypes_by_group_qvalues_plot<- ggplot(Phenotypes_groups_qvalues, aes(x = group, y = Penetrance))+
                                geom_bar(data = Phenotypes_groups, aes(fill = group), stat = "identity")+
                                geom_jitter(position = position_jitter(0.1))+
                                geom_errorbar(data = df.summary, aes(ymin = Penetrance-sd, ymax = Penetrance+sd), width = 0.5)+
                                facet_grid(cols = vars(Symptom))+
                                scale_y_continuous(breaks=seq(0,120,20), limits = c(0,120))+
                                scale_x_discrete(labels=c("1" = "E8", "2" = "E10", "3" = "E16"))+
                                scale_fill_brewer(palette = "Set2")+
                                labs(title = "Penetrance by symptom in each group",
                                     x="All patients n = 227",
                                     y="Penetrance of the symptom(%)")+
                                theme_bw()+
                                theme(legend.position = "none")+
                                stat_pvalue_manual(
                                  stat.test,
                                  y.position = 75,
                                  step.increase = 0.12,
                                  hide.ns = TRUE,
                                  label = "p.adj")# Add adj.p-value

phenotypes_by_group_qvalues_plot

save_plot(phenotypes_by_group_qvalues_plot,"_Fig3_phenotypes_by_group_qvalue.pdf",9,5)
dev.off()


Fig1 <- plot_grid(mutations_cohort_top,
                  Penetrance_alleles,
                  groups_LT,
                  agregated_groups_LT,
                  ncol=2,
                  labels = c("A","B","C","D"),
                  label_size = 22)
Fig1

save_plot(Fig1,"Fig1_cohort_description.pdf",17,10)
dev.off()

Fig2 <- plot_grid(penetrance_top,
                  allpatients_ss,
                  sexgroups_ss_box_plot,
                  agegroups_ss_box_plot,
                  subgroups_ss_box_plot,
                  ncol=2,
                  labels = c("A","B","C","D","E"),
                  label_size = 22)
Fig2

save_plot(Fig2,"Fig2_SS_analysis.pdf",10,15)
dev.off()

Fig3 <- plot_grid(subgroups_ss_ages_box_plot,
                  phenotypes_by_group_qvalues_plot,
                  ncol=1,
                  labels = c("A","B"),
                  label_size = 22)
Fig3

save_plot(Fig3,"Fig3_subgroups_phenotype_analysis.pdf",10,10)
dev.off()
