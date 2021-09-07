# from the simulated transitions work out probability of transition to a
# state for cases and controls
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

suppressPackageStartupMessages({
  library(R.utils)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  library(ROCR)
  library(colorspace)
  library(tidyr)
  library(xtable)
  library(scales)
  library(caTools)
  library(caret)
})

args=commandArgs(trailingOnly=TRUE, asValues=TRUE)
filepath=as.character(args$filepath)
name=as.character(args$name)
n_params=as.numeric(args$n_params)
n_trans=as.numeric(args$n_trans)
n_iter=as.character(format(args$n_iter, scientific=F))
type=as.character(args$type)

filepath="/rds/general/user/aj1520/home/summer_project/covid_msm/"
name="iwrd"
status_filepath=paste0(filepath, "configs/", name, "/status.txt")
visual_filepath=paste0(filepath, "results/model_comparisons/")
covars_filepath=paste0(filepath, "configs/", name, "/covars.txt")
n_iter=1000000

# define transitions and states
status_names=c("I", "W", "R", "D")
to_trans_names=c("IR", "IW", "WR", "WD")
trans_names=c("II", "IW", "IR", "WW", "WR", "WD")

# get data
status=read.table(status_filepath, skip=2)
colnames(status) = status_names

# get the indicator list
trans_list = list()
for(ind in 1:nrow(status)){
  ind_trans=c("IR"=0, "IW"=0, "WR"=0,"WD"=0)
  for(j in 1:ncol(status)){
    trans_from = status[ind,j]
    if(trans_from==1){
      from=status_names[j]
      k=j+1
      while(k<=length(status_names)){
        trans_to = status[ind,k]
        if(trans_to==0){
          k=k+1
          next
        }
        else{
          to=status_names[k]
          trans_name=paste0(from, to)
          ind_trans[trans_name]=1
          k=k+1
          break
        }
      }
    }
  }
  trans_list[[ind]] = ind_trans
}

transition_indicator = data.frame()
for(i in 1:nrow(status)){transition_indicator = rbind(transition_indicator, trans_list[[i]])}
colnames(transition_indicator) = c("Recovered", "Hospitalised",
                                   "Recovered_hospitalisation",
                                   "Died_hospitalisation")
transition_indicator$Recovered = factor(transition_indicator$Recovered)
transition_indicator$Hospitalised = factor(transition_indicator$Hospitalised)
transition_indicator$Recovered_hospitalisation = factor(transition_indicator$Recovered_hospitalisation)
transition_indicator$Died_hospitalisation = factor(transition_indicator$Died_hospitalisation)

transition_indicator = cbind(individual=rownames(transition_indicator), transition_indicator)

model_list = c("time", "medical", "health", "social", "demo")
categories=c("Model 5", "Model 4", "Model 3", "Model 2", "Model 1")
# model_list = c("medical", "health", "social", "demo")
# categories=c("Model 4", "Model 3", "Model 2", "Model 1")
transitions = data.frame()

for(i in 1:length(model_list)){
  transitions_folder = list.files(paste0(filepath, "results/probabilities"),
                                  pattern=paste0(as.character(format(n_iter, scientific=F)), "_", model_list[i]),
                                  full.names = TRUE)
  model_transitions = read.table(paste0(transitions_folder, "/transitions.txt"), header=F,stringsAsFactors=F)
  colnames(model_transitions) = trans_names
  model_transitions = model_transitions[,to_trans_names]
  
  model_transitions = cbind(individual=as.integer(rownames(model_transitions)),
                            model_transitions)
  
  model_transitions = merge(transition_indicator, model_transitions, by="individual", sort=FALSE)

  model_transitions = cbind("Category" = categories[i],
                            model_transitions)
  
  transitions = rbind(transitions, model_transitions)
}

summary(transitions)

#### Case fatality

covars=read.table(covars_filepath, skip=2)
colnames(covars) = c("Age", "Sex", "Ethnicity", "BMI", "smoking_status", "alcohol_status",
                     "education", "accommodation",
                     "Cancer", "Cardiovascular", "Hypertension", "Diabetes", "Respiratory", "Autoimmune",
                     "ACE", "Angiotensin", "Steroids", "Statins")
breaks=c(45, 55, 65, 75, 85)
covars$Age = cut(covars$Age, breaks=breaks, right = FALSE)
levels(covars$Age) = c("45-55", "55-65", "65-75", "75-85")

fatality = transitions[transitions$Category == "Model 5",]
fatality$Age = covars$Age

head(fatality)

print("Case-hospitalisation")
mean(ifelse(fatality$Hospitalised == "Hospitalised", 1, 0))
mean(fatality$IW)

print("Case-fatality")
mean(ifelse(fatality$Died_hospitalisation == "Died", 1, 0))
mean(fatality$WD)
str(fatality)

total_IHR = c("Total",
          round(mean(as.numeric(as.character(fatality$Hospitalised))),3),
          round(mean(fatality$IW),3),
          round(mean(as.numeric(as.character(fatality$Died_hospitalisation))),3),
          round(mean(fatality$WD), 3))
age_IHR = fatality %>%
  group_by(Age=as.character(Age)) %>%
  summarise("Cohort IHR" = round(mean(as.numeric(as.character(Hospitalised))), 3),
            "Model IHR" = round(mean(IW), 3),
            "Cohort IFR" = round(mean(as.numeric(as.character(Died_hospitalisation))), 3),
            "Model IFR" = round(mean(WD), 3))
age_IHR=rbind(age_IHR, total_IHR)
age_IHR[,2:5] = apply(age_IHR[,2:5], 2, as.numeric)
age_IHR[,2:5] = apply(age_IHR[,2:5], 2, percent, accuracy=0.1)

print(xtable(age_IHR, type = "latex"), include.rownames=FALSE, file = "../results/model_comparisons/age_IHR.tex")

#### Case-control plots

levels(transitions$Hospitalised) = c("Not Hospitalised", "Hospitalised")
levels(transitions$Died_hospitalisation) = c("Recovered", "Died")

# normalise
transitions$IW = transitions$IW / (transitions$IW + transitions$IR)
transitions$IW = ifelse(is.na(transitions$IW), 0, transitions$IW)

transitions$WD = transitions$WD / (transitions$WD + transitions$WR)
transitions$WD = ifelse(is.na(transitions$WD), 0, transitions$WD)

mycolours=darken(c("navy", "tomato"), amount=0.2)

head(transitions)

# median
median_IW = transitions %>%
  group_by(Category, Hospitalised) %>%
  summarise(median_IW = median(IW)) %>%
  pivot_wider(names_from = Hospitalised, values_from = median_IW)

# median
median_WD = transitions %>%
  filter(Hospitalised == "Hospitalised") %>%
  group_by(Category, Died_hospitalisation) %>%
  summarise(median_WD = median(WD)) %>%
  pivot_wider(names_from = Died_hospitalisation, values_from = median_WD)

median_df = inner_join(median_IW, median_WD, by="Category")
median_df[,2:5] = round(median_df[,2:5], 3)

write.csv(median_df, paste0(visual_filepath, "case_medians.csv"))

#### Confusion matrices


conf_df = transitions[transitions$Category == "Model 5",]
table(as.numeric(conf_df$Hospitalised))


# IW
p_hospital = ifelse(conf_df$IW > 0.5, 1, 0)
p_hclass = factor(p_hospital)
levels(p_hclass) = c("Not Hospitalised", "Hospitalised")
confusionMatrix(p_hclass, conf_df$Hospitalised)

# WD
head(conf_df)
conf_iw_df = conf_df[conf_df$Hospitalised=="Hospitalised",]
p_death = ifelse(conf_iw_df$WD > 0.5, 1, 0)
p_dclass = factor(p_death)
levels(p_dclass) = c("Recovered", "Died")
confusionMatrix(p_dclass, conf_iw_df$Died_hospitalisation)

#### Get the plots

p_hospitalisations = ggplot(transitions, aes(x=IW, linetype=Hospitalised)) +
  geom_density(size=0.2) +
  xlim(0,1) +
  ylim(0,16) +
  facet_wrap(Category~., ncol=1, scales="free_x") +
  theme_minimal() +
  ggtitle("O-I") +
  labs(x="Probability") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks.x = element_line(colour = "black", size = 0.2),
        axis.text.x = element_text(size=6),
        legend.position="bottom",
        legend.title=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_linetype_manual(values=c("dotted", "solid"))
p_hospitalisations

p_death = transitions %>%
  filter(Hospitalised == "Hospitalised") %>%
  ggplot(aes(x=IW, linetype=Died_hospitalisation)) +
  geom_density(size=0.2) +
  xlim(0,1) +
  ylim(0,16) +
  facet_wrap(Category~., ncol=1, scales="free_x") +
  theme_minimal() +
  ggtitle("I-D") +
  labs(x="Probability") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks.x = element_line(colour = "black", size = 0.2),
        axis.text.x = element_text(size=6),
        legend.position="bottom",
        legend.title=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_linetype_manual(values=c("dotted", "solid"))
p_death

#### Get the plots for the presentation

p_hospitalisations_pres = transitions %>%
  filter(Category=="Model 5") %>%
  ggplot(aes(x=IW, linetype=Hospitalised)) +
  geom_density(size=0.2) +
  xlim(0,1) +
  ylim(0,16) +
  theme_minimal() +
  ggtitle("Outpatient-Inpatient") +
  labs(x="Probability") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks.x = element_line(colour = "black", size = 0.2),
        axis.text.x = element_text(size=10),
        legend.position="bottom",
        legend.title=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_linetype_manual(values=c("dotted", "solid"))
p_hospitalisations_pres

p_death_pres = transitions %>%
  filter(Category == "Model 5") %>%
  filter(Hospitalised == "Hospitalised") %>%
  ggplot(aes(x=IW, linetype=Died_hospitalisation)) +
  geom_density(size=0.2) +
  xlim(0,1) +
  ylim(0,16) +
  facet_wrap(Category~., ncol=1, scales="free_x") +
  theme_minimal() +
  ggtitle("Inpatient-Death") +
  labs(x="Probability") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks.x = element_line(colour = "black", size = 0.2),
        axis.text.x = element_text(size=10),
        legend.position="bottom",
        legend.title=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_linetype_manual(values=c("dotted", "solid"))
p_death_pres

#### Combine
case_plots = ggarrange(p_hospitalisations, p_death, ncol = 2)

case_plots_pres = ggarrange(p_hospitalisations_pres, p_death_pres, ncol = 2)

# add labels
case_plots = case_plots +
  theme(plot.margin=unit(c(0,1,0,0),"cm")) +
  annotate("text", x=1.02, y=(5/28), label= "Model 1", angle=90) +
  annotate("text", x=1.02, y=(10/28), label= "Model 2", angle=90) +
  annotate("text", x=1.02, y=(15/28), label= "Model 3", angle=90) +
  annotate("text", x=1.02, y=(20/28), label= "Model 4", angle=90) +
  annotate("text", x=1.02, y=(25/28), label= "Model 5", angle=90)

# save
ggsave(paste0(visual_filepath, "case_plots.pdf"), case_plots, device="pdf",
       width=8, height=8)
ggsave("/rds/general/user/aj1520/home/summer_project/covid_msm/results/presentation/case_plots_pres.pdf", case_plots_pres, device="pdf",
       width=12, height=4)
