#---
#title: "Analysis of genotype-phenotype associations in ALMS"
#author: "Brais Bea Mascato"
#date: '`r format(Sys.time(), "%d %m, %Y")`'
#output:
#---
  
libraries <- c('kableExtra','ggpubr','cowplot','rstatix','readxl', 'tidyverse')
lapply(libraries,library, character.only = TRUE)

Curated_ALMS1_DB <- as.data.frame(read_xlsx("./data/Curated_ALMS1_DB.xlsx",
                                    sheet = "Genotype-Phenotype"))

ALMS1_cohort <- Curated_ALMS1_DB[complete.cases(Curated_ALMS1_DB[,c(5:19,21:37)]),c(1:37)]
colnames(ALMS1_cohort)[6:15] <- c("PV1_DNA", "EA1", "PV1_PROTEIN","TM_1","PV2_DNA", "EA2", "PV2_PROTEIN","TM_2","GS","LT")
colnames(ALMS1_cohort)[37] <- "SS"

ALMS1_cohort <- ALMS1_cohort %>% mutate(group = as.factor(ifelse(LT < 9, 1,ifelse(LT >= 9 & LT < 14,2,3))))

phenotypes <- c("VI","MT","HL","HRT","LIV","REN","PUL","SHS","REP","TYD","MEND","ABFING","INT","SCO","NER","ALO")

####Functions####

source('./src/main/utils.R')

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

df_mutations <- df_mutations %>%
                  mutate(position = str_extract(mutation_protein, "\\d{2,4}(?=\\D)"))%>%
                  arrange(as.numeric(gsub("[^0-9]", "", position)))

mutations_cohort <- ggplot(df_mutations, aes(x = mutation_protein , y = count) )+
                        geom_bar(
                          # fill= "#8EA0CB",
                          stat="identity")+
                        labs(title = "Pathogenic alleles in ALMS cohort",
                             x = "Pathogenic variant",
                             y = "Absolute frequency")+
                        # scale_fill_grey()+
                        theme_classic2() +
                        coord_flip()+
                        theme(legend.position = "none")
mutations_cohort

save_plot(mutations_cohort,"_FigS1_mutations_cohort.pdf",5,30)
dev.off()

mutations_cohort_top <- ggplot(df_mutations[which(df_mutations$count>2),], aes(x = factor(mutation_protein, level = mutation_protein) , y = count) )+
                            geom_bar(stat="identity")+
                            labs(title = "Pathogenic alleles in ALMS cohort",
                                 x = "Pathogenic variant",
                                 y = "Absolute frequency")+
                            scale_y_continuous(breaks=seq(0,40,5), limits = c(0,40))+
                            # scale_fill_grey()+
                            theme_classic2() +
                            theme(legend.position = "none",
                                  axis.text.x = element_text(angle = 45, hjust=1),
                                  plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
mutations_cohort_top

save_plot(mutations_cohort_top,"_Fig1_mutations_cohort_top.pdf",9,5)
dev.off()

## prevalence of mutations in ALMS1 exons genes

df <- data.frame(table(c(ALMS1_cohort$EA1,ALMS1_cohort$EA2)))
colnames(df) <- c("exon","count")
df <- df[df$exon!="intron 17",]
df$exon <- factor(df$exon, levels = c(1,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20))
df$percentage = round(df$count/length(c(ALMS1_cohort$EA1,ALMS1_cohort$EA2))*100,digits = 2)

Penetrance_alleles <- ggplot(df, aes(x = exon, y = count) )+
                      geom_bar(stat="identity")+
                      labs(title = "Pathogenic alleles of ALMS1 gene",
                           x = "Exon with the pathogenic variant",
                           y = "Absolute frequency")+
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
                            scale_x_discrete(labels= c("1"="G1","2"="G2","3"="G3"))+
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

#Prevalence of symptoms

pdf(file = paste0("./results/figures/", "_FigS1_prevalence_global.pdf"), width = 12, height = 5)
plot_prevalence(ALMS1_cohort, "global")
dev.off()


prevalence_top <- plot_prevalence(ALMS1_cohort, "top")
save_plot(prevalence_top,"_Fig2_prevalence_top.pdf",5,5)
dev.off()

#Syndromic score distribution complete cohort

allpatients_ss <- plot_syndromic_score(ALMS1_cohort, "All patients \n n=227")
allpatients_ss

save_plot(allpatients_ss,"_Fig2_allpatients_ss.pdf",5,7)
dev.off()

#Syndromic score distribution by sex

sexgroups_ss_box_plot <- plot_syndromic_score_box(ALMS1_cohort, ALMS1_cohort$Sex)+
                          stat_compare_means(method = "wilcox.test", label.y = 1.3 )

sexgroups_ss_box_plot

save_plot(sexgroups_ss_box_plot,"_Fig2_sexgroups_ss_boxplot.pdf",6,6)
dev.off()

#Syndromic score distribution by ages

agegroups_ss_box_plot <- plot_syndromic_score_box(ALMS1_cohort, ALMS1_cohort$Age2)+
                          stat_compare_means(method = "kruskal.test", label.y = 1.3)

agegroups_ss_box_plot

save_plot(agegroups_ss_box_plot,"_Fig2_agegroups_ss_box_plot.pdf",6,6)
dev.off()

#Syndromic score distribution in subgroups (G1, G2, G2)

#Syndromic score distribution in subgroups
stat.test <- compare_means(SS ~ group,  data = ALMS1_cohort, p.adjust.method ="BH")


subgroups_ss_box_plot <- plot_syndromic_score_box(ALMS1_cohort, ALMS1_cohort$group)+
                          scale_x_discrete(labels=c("1" = "G1", "2" = "G2", "3" = "G3"))+
                                stat_pvalue_manual(
                                  stat.test, 
                                  y.position = 1.1, step.increase = 0.1,
                                  label = "p.adj"
                                ) + # Add adj.p-value
                          stat_compare_means(method = "kruskal.test", label.y = 1.3) # Add global p-value
  
subgroups_ss_box_plot

save_plot(subgroups_ss_box_plot,"_Fig2_subgroups_ss_boxplot.pdf",6,6)
dev.off()

#Correlation plot syndromic score - ages

agegroups_ss_cor_plot <- ggplot(ALMS1_cohort, aes(AGE_ORIG, SS))+
                          ggtitle("Correlation between syndromic score and age in Alström syndrome patients")+
                          xlab("Age (Years)")+
                          ylab("Syndromic Score")+
                          geom_smooth(method=lm)+
                          stat_cor(label.x = 30, label.y = 0.3)+
                          geom_jitter()+
                          theme_bw()

agegroups_ss_cor_plot

save_plot(agegroups_ss_cor_plot,"_FigS4_agegroups_ss_cor_plot.pdf",7,7)
dev.off()

#Syndromic score distribution in subgroups by ages

stat.test <- ALMS1_cohort %>%
              group_by(Age2)%>%
              pairwise_wilcox_test(SS ~ group, p.adjust.method ="BH")

stat.test

subgroups_ss_ages_box_plot <- plot_syndromic_score_box(ALMS1_cohort, ALMS1_cohort$group)+
                                  scale_x_discrete(labels=c("1" = "G1", "2" = "G2", "3" = "G3"))+
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

#prevalence of mains phenotypes features in subgroups (G1, G2, G2)

group_exon8 <- ALMS1_cohort[ALMS1_cohort$group==1,]
group_exon10 <- ALMS1_cohort[ALMS1_cohort$group==2,]
group_exon16 <- ALMS1_cohort[ALMS1_cohort$group==3,]

Phenotypes_groups <- data.frame( "Prevalence"=c(phenotype_analysis(group_exon8$VI)$percent,
                                                phenotype_analysis(group_exon8$MT)$percent,
                                                phenotype_analysis(group_exon8$HL)$percent,
                                                phenotype_analysis(group_exon8$HRT)$percent,
                                                phenotype_analysis(group_exon8$LIV)$percent,
                                                phenotype_analysis(group_exon8$REN)$percent,
                                                phenotype_analysis(group_exon8$MEND)$percent,
                                                phenotype_analysis(group_exon8$PUL)$percent,
                                                phenotype_analysis(group_exon8$REP)$percent,
                                               phenotype_analysis(group_exon10$VI)$percent,
                                               phenotype_analysis(group_exon10$MT)$percent,
                                               phenotype_analysis(group_exon10$HL)$percent,
                                               phenotype_analysis(group_exon10$HRT)$percent,
                                               phenotype_analysis(group_exon10$LIV)$percent,
                                               phenotype_analysis(group_exon10$REN)$percent,
                                               phenotype_analysis(group_exon10$MEND)$percent,
                                               phenotype_analysis(group_exon10$PUL)$percent,
                                               phenotype_analysis(group_exon10$REP)$percent,
                                               phenotype_analysis(group_exon16$VI)$percent,
                                               phenotype_analysis(group_exon16$MT)$percent,
                                               phenotype_analysis(group_exon16$HL)$percent,
                                               phenotype_analysis(group_exon16$HRT)$percent,
                                               phenotype_analysis(group_exon16$LIV)$percent,
                                               phenotype_analysis(group_exon16$REN)$percent,
                                               phenotype_analysis(group_exon16$MEND)$percent,
                                               phenotype_analysis(group_exon16$PUL)$percent,
                                               phenotype_analysis(group_exon16$REP)$percent),
                                    
                                  "Symptom" =c(rep(c("VI","MT","HL","HRT","LIV","REN","MEND","PUL","REP"),3)),
                                    
                                   "group" =as.factor(c(rep(1,9),rep(2,9),rep(3,9)))
                                 ) 

Phenotypes_groups$Symptom <- factor(Phenotypes_groups$Symptom, levels=c("VI","MT","HL","HRT","LIV","REN","MEND","PUL","REP"))

#Plot prevalence of mains phenotypes features in subgroups (G1, G2, G2) using beta-distributions

Phenotypes_groups_qvalues <- data.frame( "Prevalence"=c(unlist(phenotype_analysis(group_exon8$VI)[,5:7]),
                                                unlist(phenotype_analysis(group_exon8$MT)[,5:7]),
                                                unlist(phenotype_analysis(group_exon8$HL)[,5:7]),
                                                unlist(phenotype_analysis(group_exon8$HRT)[,5:7]),
                                                unlist(phenotype_analysis(group_exon8$LIV)[,5:7]),
                                                unlist(phenotype_analysis(group_exon8$REN)[,5:7]),
                                                unlist(phenotype_analysis(group_exon8$MEND)[,5:7]),
                                                unlist(phenotype_analysis(group_exon8$PUL)[,5:7]),
                                                unlist(phenotype_analysis(group_exon8$REP)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$VI)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$MT)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$HL)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$HRT)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$LIV)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$REN)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$MEND)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$PUL)[,5:7]),
                                                unlist(phenotype_analysis(group_exon10$REP)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$VI)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$MT)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$HL)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$HRT)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$LIV)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$REN)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$MEND)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$PUL)[,5:7]),
                                                unlist(phenotype_analysis(group_exon16$REP)[,5:7])),
                                 
                                 "Symptom" =c(rep(c(rep("VI",3),rep("MT",3),rep("HL",3),rep("HRT",3),rep("LIV",3),rep("REN",3),rep("MEND",3),rep("PUL",3),rep("REP",3)),3)),
                                 
                                 "group" =as.factor(c(rep(1,27),rep(2,27),rep(3,27)))
) 

Phenotypes_groups_qvalues$Symptom <- factor(Phenotypes_groups_qvalues$Symptom, levels=c("VI","MT","HL","HRT","LIV","REN","MEND","PUL","REP"))

df.summary <- Phenotypes_groups_qvalues %>%
                  group_by(group, Symptom) %>%
                  summarise(
                    sd = sd(Prevalence, na.rm = TRUE),
                    Prevalence = mean(Prevalence)
                  )
df.summary

for(i in c(21:27,29,31)){
  
  if(i==21){df <- data.frame()}
  
  df <- rbind(df,data.frame(list(Symptom=phenotypes[i-20], 
                                 .y.= rep("Prevalence",3), 
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
stat.test$Symptom <- factor(stat.test$Symptom, levels=c("VI","MT","HL","HRT","LIV","REN","MEND","PUL","REP"))

kable(stat.test, format = "html") %>%
  kable_styling(full_width = F, font_size = 9,bootstrap_options = c("striped", "hover", "condensed", "responsive"))

phenotypes_by_group_qvalues_plot<- ggplot(Phenotypes_groups_qvalues, aes(x = group, y = Prevalence))+
                                geom_bar(data = Phenotypes_groups, aes(fill = group), stat = "identity")+
                                geom_jitter(position = position_jitter(0.1))+
                                geom_errorbar(data = df.summary, aes(ymin = Prevalence-sd, ymax = Prevalence+sd), width = 0.5)+
                                facet_grid(cols = vars(Symptom))+
                                scale_y_continuous(breaks=seq(0,120,20), limits = c(0,120))+
                                scale_x_discrete(labels=c("1" = "G1", "2" = "G2", "3" = "G3"))+
                                scale_fill_brewer(palette = "Set2")+
                                labs(title = "Prevalence by symptom in each group",
                                     x="All patients n = 227",
                                     y="Prevalence of the symptom(%)")+
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

#Penetrance of mains phenotypes features in sexgroups (F, M)

group_F <- ALMS1_cohort[ALMS1_cohort$Sex=="F",]
group_M <- ALMS1_cohort[ALMS1_cohort$Sex=="M",]

Phenotypes_sex_groups <- data.frame( "Prevalence"=c(phenotype_analysis(group_F$VI)$percent,
                                                phenotype_analysis(group_F$MT)$percent,
                                                phenotype_analysis(group_F$HL)$percent,
                                                phenotype_analysis(group_F$HRT)$percent,
                                                phenotype_analysis(group_F$LIV)$percent,
                                                phenotype_analysis(group_F$REN)$percent,
                                                phenotype_analysis(group_F$MEND)$percent,
                                                phenotype_analysis(group_F$PUL)$percent,
                                                phenotype_analysis(group_F$REP)$percent,
                                                phenotype_analysis(group_M$VI)$percent,
                                                phenotype_analysis(group_M$MT)$percent,
                                                phenotype_analysis(group_M$HL)$percent,
                                                phenotype_analysis(group_M$HRT)$percent,
                                                phenotype_analysis(group_M$LIV)$percent,
                                                phenotype_analysis(group_M$REN)$percent,
                                                phenotype_analysis(group_M$MEND)$percent,
                                                phenotype_analysis(group_M$PUL)$percent,
                                                phenotype_analysis(group_M$REP)$percent),
                                 
                                 "Symptom" =c(rep(c("VI","MT","HL","HRT","LIV","REN","MEND","PUL","REP"),2)),
                                 
                                 "group" =as.factor(c(rep("F",9),rep("M",9)))
) 

Phenotypes_sex_groups$Symptom <- factor(Phenotypes_groups$Symptom, levels=c("VI","MT","HL","HRT","LIV","REN","MEND","PUL","REP"))


#Plot prevalence of mains phenotypes features in sex groups (F,M) using beta-distributions

Phenotypes_sex_groups_qvalues <- data.frame( "Prevalence"=c(unlist(phenotype_analysis(group_F$VI)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$MT)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$HL)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$HRT)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$LIV)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$REN)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$MEND)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$PUL)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$REP)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$VI)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$MT)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$HL)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$HRT)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$LIV)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$REN)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$MEND)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$PUL)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$REP)[,5:7])),
                                         
                                         "Symptom" =c(rep(c(rep("VI",3),rep("MT",3),rep("HL",3),rep("HRT",3),rep("LIV",3),rep("REN",3),rep("MEND",3),rep("PUL",3),rep("REP",3)),2)),
                                         
                                         "group" =as.factor(c(rep("F",27),rep("M",27)))
) 

Phenotypes_sex_groups_qvalues$Symptom <- factor(Phenotypes_groups_qvalues$Symptom, levels=c("VI","MT","HL","HRT","LIV","REN","MEND","PUL","REP"))

df.summary <- Phenotypes_sex_groups_qvalues %>%
  group_by(group, Symptom) %>%
  summarise(
    sd = sd(Prevalence, na.rm = TRUE),
    Prevalence = mean(Prevalence)
  )
df.summary

for(i in c(21:27,29,31)){
  
  if(i==21){df <- data.frame()}
  
  df <- rbind(df,data.frame(list(Symptom=phenotypes[i-19], 
                                 .y.= rep("Prevalence",1), 
                                 pairwise_fisher_test(matrix( c(
                                   phenotype_analysis(group_F[i])$no,
                                   phenotype_analysis(group_M[i])$no,
                                   phenotype_analysis(group_F[i])$yes,
                                   phenotype_analysis(group_M[i])$yes),
                                   byrow = FALSE, 
                                   ncol = 2,
                                   dimnames = list(c("F","M"), c("No","Yes"))),
                                   p.adjust.method = "BH"))))
  df                         
}

stat.test <- as_tibble(df)
stat.test$Symptom <- factor(stat.test$Symptom, levels=c("VI","MT","HL","HRT","LIV","REN","MEND","PUL","REP"))

kable(stat.test, format = "html") %>%
  kable_styling(full_width = F, font_size = 9,bootstrap_options = c("striped", "hover", "condensed", "responsive"))

phenotypes_by_sex_group_qvalues_plot<- ggplot(Phenotypes_sex_groups_qvalues, aes(x = group, y = Prevalence))+
                                        geom_bar(data = Phenotypes_sex_groups, aes(fill = group), stat = "identity")+
                                        geom_jitter(position = position_jitter(0.1))+
                                        geom_errorbar(data = df.summary, aes(ymin = Prevalence-sd, ymax = Prevalence+sd), width = 0.5)+
                                        facet_grid(cols = vars(Symptom))+
                                        scale_y_continuous(breaks=seq(0,120,20), limits = c(0,120))+
                                        scale_x_discrete()+
                                        scale_fill_brewer(palette = "Set2")+
                                        labs(title = "Prevalence by symptom in each group",
                                             x="All patients n = 227",
                                             y="Prevalence of the symptom(%)")+
                                        theme_bw()+
                                        theme(legend.position = "none")+
                                        stat_pvalue_manual(
                                          stat.test,
                                          y.position = 75,
                                          step.increase = 0.12,
                                          hide.ns = TRUE,
                                          label = "p.adj")# Add adj.p-value

phenotypes_by_sex_group_qvalues_plot

save_plot(phenotypes_by_sex_group_qvalues_plot,"_FigS3_phenotypes_by_sex_group_qvalue.pdf",9,5)
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

Fig2 <- plot_grid(prevalence_top,
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
