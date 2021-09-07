rm(list=ls())

suppressPackageStartupMessages({
  library(R.utils)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(ggpubr)
  library(colorspace)
  library(scales)
})

args=commandArgs(trailingOnly=TRUE, asValues=TRUE)
filepath=as.character(args$filepath)
name=as.character(args$name)
n_params=as.numeric(args$n_params)
n_trans=as.numeric(args$n_trans)
n_iter=as.character(format(args$n_iter, scientific=F))
type=as.character(args$type)
start_date = '2020-09-01'

filepath="/Users/adamjackson/Documents/summer_project/covid_msm/"
filepath="/rds/general/user/aj1520/home/summer_project/covid_msm/"
name="iwrd"
seed=1
n_params=52
n_trans=6
n_iter="1000000"
type="time"

covars_filepath=paste0(filepath, "configs/", name, "/covars.txt")
prob_filepath=paste0(filepath, "results/probabilities/", name, "_", n_params, "_", n_trans, "_", n_iter, "_", type, "/")

# probability names
prob_names = c("P_I_I", "P_I_R", "P_I_W",
               "P_W_W", "P_W_R", "P_W_D")
# breaks=c(min(covars$age), -1, 0, 1, max(covars$age))
breaks=c(45, 55, 65, 75, 85)

# covars
covars=read.table(covars_filepath, skip=2)
head(covars)
colnames(covars) = c("Age", "Sex", "Ethnicity", "BMI", "smoking_status", "alcohol_status",
                     "education", "accommodation",
                     "Cancer", "Cardiovascular", "Hypertension", "Diabetes", "Respiratory", "Autoimmune",
                     "ACE", "Angiotensin", "Steroids", "Statins")

covars$Age = cut(covars$Age, breaks=breaks, right = FALSE)
levels(covars$Age) = c("45-55", "55-65", "65-75", "75-85")

# factors
levels(covars$Sex) = c("Female", "Male")
levels(covars$Ethnicity) = c("White", "Other", "Black")
levels(covars$smoking_status) = c("Never", "Previous", "Current")
levels(covars$alcohol_status) = c("Never", "Previous", "Current")
levels(covars$education) = c("High", "Intermediate", "Low")
levels(covars$accommodation) = c("House", "Flat")

covars$Comorbidity = rowSums(covars[,c("Cancer", "Cardiovascular", "Hypertension", "Diabetes", "Respiratory", "Autoimmune")])
covars$Comorbidity = ifelse(covars$Comorbidity>0,1,0)
head(covars)

# check the number of persons per group
age_bar = covars %>%
  group_by(Age) %>%
  summarise(n=n()) %>%
  ggplot(aes(x=Age, y=n)) +
  geom_bar(stat="identity")

# get each probability matrix and combine them
all_probabilities = data.frame()
for (prob_name in prob_names){
  probabilities = read.table(list.files(prob_filepath, pattern = prob_name ,full.names = TRUE),header=F,stringsAsFactors=F)
  colnames(probabilities) = seq(as.Date(start_date), by = "day", length.out = ncol(probabilities))
  
  # join with covars
  probabilities = cbind(prob_name, covars, probabilities)
  
  # add to all probabilities
  all_probabilities = rbind(all_probabilities, probabilities)
}

# adjust
levels(all_probabilities$Age) = c("45-55", "55-65", "65-75", "75-85")
mycolours=darken(c("navy", "orange", "forestgreen", "tomato"), amount=0.2)

##### Age/sex


# melt, group by and average
probs_age_sex = all_probabilities %>% 
  dplyr::select(-Ethnicity, -BMI, -smoking_status, -alcohol_status, -education, 
                -accommodation,
                -Cancer, -Cardiovascular, -Hypertension, -Diabetes, -Respiratory, -Comorbidity, -Autoimmune,
                -ACE, -Angiotensin, -Steroids, -Statins) %>%
  melt(id.vars = c('prob_name', 'Age', 'Sex'), 
       variable.name="date", value.name='probability') %>%
  group_by(prob_name, Age, Sex, date) %>%
  summarise(probability=mean(probability)) %>%
  ungroup()
head(probs_age_sex)

# convert
probs_age_sex$date = as.Date(probs_age_sex$date)
probs_age_sex$date = as.POSIXct(probs_age_sex$date)
probs_age_sex$time = 1 / probs_age_sex$probability
colnames(probs_age_sex) = c("name", "Age", "Sex", "date", "probability", "time")
probs_age_sex = probs_age_sex %>% filter(name %in% c("P_I_W", "P_W_D"))
probs_age_sex$name = factor(probs_age_sex$name)
levels(probs_age_sex$name)=c("Outpatient-Inpatient", "Inpatient-Death")

# plot
plot_age_sex = probs_age_sex %>%
  ggplot(aes(x=date, y=probability, colour=Age, linetype=Sex)) +
  geom_line() +
  theme_minimal() +
  labs(x="Date", y="Probability") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks.x = element_line(colour = "black", size = 0.2),
        axis.title.x = element_blank(),
        legend.position="bottom",
        legend.title=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14, hjust = 0)) +
  scale_linetype_manual(values=c("dotted", "solid"))+
  scale_color_manual(values=mycolours)+
  scale_x_datetime(breaks=seq.POSIXt(min(probs_age_sex$date), max(probs_age_sex$date), by = "1 month"),
                   labels = date_format("%b-%y")) +
  facet_wrap(name~., ncol=2, scales = "free")
plot_age_sex

# plot
plot_age_sex_time = probs_age_sex %>%
  ggplot(aes(x=date, y=time, colour=Age, linetype=Sex)) +
  geom_line() +
  theme_minimal() +
  labs(x="Date", y="Days") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks.x = element_line(colour = "black", size = 0.2),
        axis.title.x = element_blank(),
        legend.position="bottom",
        legend.title=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14, hjust = 0)) +
  scale_linetype_manual(values=c("dotted", "solid"))+
  scale_color_manual(values=mycolours)+
  scale_x_datetime(breaks=seq.POSIXt(min(probs_age_sex$date), max(probs_age_sex$date), by = "1 month"),
                   labels = date_format("%b-%y")) +
  facet_wrap(name~., ncol=2, scales = "free")
plot_age_sex_time

#### Age + Comorbidity

# melt, group by and average
probs_age_comorb = all_probabilities %>% 
  dplyr::select(-Sex, -Ethnicity, -BMI, -smoking_status, -alcohol_status, -education, 
                -accommodation,
                -Cancer, -Cardiovascular, -Hypertension, -Diabetes, -Respiratory, -Autoimmune,
                -ACE, -Angiotensin, -Steroids, -Statins) %>%
  melt(id.vars = c('prob_name', 'Age', 'Comorbidity'), 
       variable.name="date", value.name='probability') %>%
  group_by(prob_name, Age, Comorbidity, date) %>%
  summarise(probability=mean(probability)) %>%
  ungroup()

# edit
probs_age_comorb$date = as.Date(probs_age_comorb$date)
probs_age_comorb$date = as.POSIXct(probs_age_comorb$date)
probs_age_comorb$time = 1 / probs_age_comorb$probability
colnames(probs_age_comorb) = c("name", "Age", "Comorbidity", "date", "probability", "time")
probs_age_comorb = probs_age_comorb %>% filter(name %in% c("P_I_W", "P_W_D"))
probs_age_comorb$name = factor(probs_age_comorb$name)
probs_age_comorb$Comorbidity = factor(probs_age_comorb$Comorbidity)
levels(probs_age_comorb$Comorbidity)=c("No Comorbidites", "1+ Comorbidities")
levels(probs_age_comorb$name)=c("Outpatient-Inpatient", "Inpatient-Death")

# plot
plot_age_comorb = probs_age_comorb %>%
  ggplot(aes(x=date, y=probability, colour=Age, linetype=Comorbidity)) +
  geom_line() +
  theme_minimal() +
  labs(x="Date", y="Probability") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks.x = element_line(colour = "black", size = 0.2),
        legend.position="bottom",
        legend.title=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_linetype_manual(values=c("dotted", "solid"))+
  scale_color_manual(values=mycolours)+
  scale_x_datetime(breaks=seq.POSIXt(min(probs_age_comorb$date), max(probs_age_comorb$date), by = "1 month"),
                   labels = date_format("%b-%y")) +
  facet_wrap(name~., ncol=2, scales = "free")
plot_age_comorb
# plot time
plot_age_comorb_time = probs_age_comorb %>%
  ggplot(aes(x=date, y=time, colour=Age, linetype=Comorbidity)) +
  geom_line() +
  theme_minimal() +
  labs(x="Date", y="Days") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks.x = element_line(colour = "black", size = 0.2),
        legend.position="bottom",
        legend.title=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_linetype_manual(values=c("dotted", "solid"))+
  scale_color_manual(values=mycolours)+
  scale_x_datetime(breaks=seq.POSIXt(min(probs_age_comorb$date), max(probs_age_comorb$date), by = "1 month"),
                   labels = date_format("%b-%y")) +
  facet_wrap(name~., ncol=2, scales = "free")
plot_age_comorb_time
#### Save

prob_plot = ggarrange(plot_age_sex, plot_age_comorb, 
                      common.legend = FALSE, nrow=2)
prob_plot
time_plot = ggarrange(plot_age_sex_time, plot_age_comorb_time, 
                      common.legend = FALSE, nrow=2)
# add labels
prob_plot = prob_plot +
  theme(plot.margin=unit(c(0,1,0,0),"cm")) +
  annotate("text", x=1.02, y=(7/24), label= "Age + Comorbidity", angle=90, size=5) +
  annotate("text", x=1.02, y=(19/24), label= "Age + Sex", angle=90, size=5)

# add labels
time_plot = time_plot +
  theme(plot.margin=unit(c(0,1,0,0),"cm")) +
  annotate("text", x=1.02, y=(7/24), label= "Age + Comorbidity", angle=90, size=5) +
  annotate("text", x=1.02, y=(19/24), label= "Age + Sex", angle=90, size=5)

ggsave(paste0(filepath, "results/model_comparisons/", "prob_plot.pdf"), prob_plot, device="pdf",
       height=8, width=8)

ggsave(paste0(filepath, "results/model_comparisons/", "time_plot.pdf"), time_plot, device="pdf",
       height=8, width=8)

ggsave(paste0(filepath, "results/presentation/", "prob_age_sex.pdf"), plot_age_sex, device="pdf",
       width=12, height=4)

