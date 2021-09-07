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
  library(colorspace)
})

filepath="/rds/general/user/aj1520/home/summer_project/covid_msm/"
name="iwrd"
n_params=52
n_trans=6
n_iter="1000000"
type="time"

mycolours=lighten(c("red", "royalblue"), amount=0.2)
mycolours

covars_filepath=paste0(filepath, "configs/", name, "/covars.txt")
conditional_filepath=paste0(filepath, "results/probabilities/", name, "_", n_params, "_", n_trans, "_", n_iter, "_", type, "/conditional.txt")
visual_filepath=paste0(filepath, "results/mcmc/", name, "_", n_params, "_", n_trans, "_", n_iter, "_", type, "/figures/")
dates_filepath=paste0(filepath, "configs/", name, "/dates.txt")

# define transitions and states
trans_names=c("II", "IW", "IR", "WW", "WR", "WD")
to_trans_names=c("IR", "IW", "WR", "WD")

# get data
predicted=read.table(conditional_filepath)
head(predicted)

# reset the column names
colnames(predicted) = trans_names
predicted = predicted[,to_trans_names]

head(predicted)

# get the covariates
covars=read.table(covars_filepath, skip=2)
colnames(covars) = c("Age", "Sex", "Ethnicity", "BMI", "smoking_status", "alcohol_status",
                     "education", "accommodation",
                     "Cancer", "Cardiovascular", "Hypertension", "Diabetes", "Respiratory", "Autoimmune",
                     "ACE", "Angiotensin", "Steroids", "Statins")
breaks=c(45, 55, 65, 75, 85)
covars$Age = cut(covars$Age, breaks=breaks, right = FALSE)
levels(covars$Age) = c("45-55", "55-65", "65-75", "75-85")

# dates
dates=read.table(dates_filepath, skip=2)
colnames(dates) = c("I", "W", "R", "D")
true=dates
true$IR = ifelse((dates$W == -1) & (dates$R != -1),
                 dates$R-dates$I, NA)
true$IW = ifelse((dates$W != -1),
                 dates$W-dates$I, NA)
true$WR = ifelse((dates$W != -1) & (dates$R != -1),
                 dates$R-dates$W, NA)
true$WD = ifelse((dates$W != -1) & (dates$D != -1),
                 dates$D-dates$W, NA)
true = true[,to_trans_names]


# get the individual column and melt
predicted = cbind(individual=rownames(predicted), predicted)
predicted$Age = covars$Age
predicted_melt = melt(predicted, id.vars=c("individual", "Age"))
head(predicted_melt)
colnames(predicted_melt) = c("individual", "Age", "transition", "predicted")
predicted_melt$Type = "Simulated"

# get the individual column and melt
true = cbind(individual=rownames(true), true)
true$Age = covars$Age
true_melt = melt(true, id.vars=c("individual", "Age"))
head(true_melt)
colnames(true_melt) = c("individual", "Age", "transition", "predicted")
true_melt$Type = "True"

melt = rbind(predicted_melt, true_melt)

### Numerical tables
IW_median = melt %>%
  filter(transition == "IW") %>%
  group_by(Type, Age) %>%
  summarise(median_IW=median(predicted, na.rm=TRUE),
            lower_IW=quantile(predicted, 0.25, na.rm=TRUE),
            upper_IW=quantile(predicted, 0.75, na.rm=TRUE)) %>%
  mutate("O-I" = paste0(format(round(median_IW, 1), nsmall = 1),
                        " [", format(round(lower_IW, 1), nsmall = 1), "-",
                        format(round(upper_IW, 1), nsmall = 1), "]")) %>%
  select(Type, Age, `O-I`)

IR_median = melt %>%
  filter(transition == "IR") %>%
  group_by(Type, Age) %>%
  summarise(median_IR=median(predicted, na.rm=TRUE),
            lower_IR=quantile(predicted, 0.25, na.rm=TRUE),
            upper_IR=quantile(predicted, 0.75, na.rm=TRUE)) %>%
  mutate("O-R" = paste0(format(round(median_IR, 1), nsmall = 1),
                        " [", format(round(lower_IR, 1), nsmall = 1), "-",
                        format(round(upper_IR, 1), nsmall = 1), "]")) %>%
  select(Type, Age, `O-R`)

WD_median = melt %>%
  filter(transition == "WD") %>%
  group_by(Type, Age) %>%
  summarise(median_WD=median(predicted, na.rm=TRUE),
            lower_WD=quantile(predicted, 0.25, na.rm=TRUE),
            upper_WD=quantile(predicted, 0.75, na.rm=TRUE)) %>%
  mutate("I-D" = paste0(format(round(median_WD, 1), nsmall = 1),
                        " [", format(round(lower_WD, 1), nsmall = 1), "-",
                        format(round(upper_WD, 1), nsmall = 1), "]")) %>%
  select(Type, Age, `I-D`)

WR_median = melt %>%
  filter(transition == "WR") %>%
  group_by(Type, Age) %>%
  summarise(median_WR=median(predicted, na.rm=TRUE),
            lower_WR=quantile(predicted, 0.25, na.rm=TRUE),
            upper_WR=quantile(predicted, 0.75, na.rm=TRUE)) %>%
  mutate("I-R" = paste0(format(round(median_WR, 1), nsmall = 1),
                        " [", format(round(lower_WR, 1), nsmall = 1), "-",
                        format(round(upper_WR, 1), nsmall = 1), "]")) %>%
  select(Type, Age, `I-R`)

sojourn_median = IW_median %>%
  inner_join(IR_median, by=c("Type", "Age")) %>%
  inner_join(WD_median, by=c("Type", "Age")) %>%
  inner_join(WR_median, by=c("Type", "Age"))
sojourn_median
write.csv(sojourn_median, "../results/model_comparisons/sojourn_medians.csv")

#### Box plots

# all-predicted to transitions
IW_transitions = melt %>%
  filter(transition == "IW") %>%
  ggplot(aes(x=predicted, group=Type)) +
  geom_boxplot(aes(fill=Type)) +
  theme_minimal() +
  facet_wrap(Age~., ncol=1) +
  ggtitle("Outpatient-Inpatient") +
  labs(x="Days") +
  xlim(0,30) +
  scale_fill_manual(values=mycolours)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks.x = element_line(colour = "black", size = 0.2),
        axis.text.y=element_blank(),
        strip.text.x = element_blank())
IW_transitions

# all-predicted to transitions
IR_transitions = melt %>%
  filter(transition == "IR") %>%
  ggplot(aes(x=predicted, group=Type)) +
  geom_boxplot(aes(fill=Type)) +
  theme_minimal() +
  facet_wrap(Age~., ncol=1) +
  ggtitle("Outpatient-Recovery") +
  labs(x="Days") +
  xlim(0,15) +
  scale_fill_manual(values=mycolours)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks.x = element_line(colour = "black", size = 0.2),
        axis.text.y=element_blank(),
        strip.text.x = element_blank())
IR_transitions

# all-predicted to transitions
WD_transitions = melt %>%
  filter(transition == "WD") %>%
  ggplot(aes(x=predicted, group=Type)) +
  geom_boxplot(aes(fill=Type)) +
  theme_minimal() +
  facet_wrap(Age~., ncol=1) +
  ggtitle("Inpatient-Death") +
  labs(x="Days") +
  xlim(0,30) +
  scale_fill_manual(values=mycolours)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks.x = element_line(colour = "black", size = 0.2),
        axis.text.y=element_blank(),
        strip.text.x = element_blank())
WD_transitions

# all-predicted to transitions
WR_transitions = melt %>%
  filter(transition == "WR") %>%
  ggplot(aes(x=predicted, group=Type)) +
  geom_boxplot(aes(fill=Type)) +
  theme_minimal() +
  facet_wrap(Age~., ncol=1) +
  ggtitle("Inpatient-Recovery") +
  labs(x="Days") +
  xlim(0,30) +
  scale_fill_manual(values=mycolours)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks.x = element_line(colour = "black", size = 0.2),
        axis.text.y=element_blank(),
        strip.text.x = element_blank())
WR_transitions




# sojourn plots
sojourn_plot = ggarrange(IW_transitions, IR_transitions, WD_transitions, WR_transitions,
                      common.legend = TRUE, ncol=2, nrow=2)
# add labels
sojourn_plot = sojourn_plot +
  theme(plot.margin=unit(c(0,0,0,1),"cm"))

sojourn_plot = sojourn_plot +
  annotate("text", x=-0.5/20, y=(2.1/20), label= "75-85", angle=90, size=4) +
  annotate("text", x=-0.5/20, y=(4/20), label= "65-75", angle=90, size=4) +
  annotate("text", x=-0.5/20, y=(5.8/20), label= "55-65", angle=90, size=4) +
  annotate("text", x=-0.5/20, y=(7.6/20), label= "45-55", angle=90, size=4) +
  annotate("text", x=-0.5/20, y=(11.5/20), label= "75-85", angle=90, size=4) +
  annotate("text", x=-0.5/20, y=(13.4/20), label= "65-75", angle=90, size=4) +
  annotate("text", x=-0.5/20, y=(15.2/20), label= "55-65", angle=90, size=4) +
  annotate("text", x=-0.5/20, y=(17.0/20), label= "45-55", angle=90, size=4)

ggsave(paste0(filepath, "results/model_comparisons/", "sojourn_plot.pdf"),
       sojourn_plot, device="pdf",
       height=9, width=9)

# sojourn plots
sojourn_pres = ggarrange(IW_transitions, WD_transitions,
                         common.legend = TRUE, ncol=2, nrow=1)
sojourn_pres = sojourn_pres +
  theme(plot.margin=unit(c(0,0,0,1),"cm"))
sojourn_pres = sojourn_pres +
  annotate("text", x=-0.1/20, y=(3.5/20), label= "75-85", angle=90, size=4) +
  annotate("text", x=-0.1/20, y=(7.5/20), label= "65-75", angle=90, size=4) +
  annotate("text", x=-0.1/20, y=(11.5/20), label= "55-65", angle=90, size=4) +
  annotate("text", x=-0.1/20, y=(15.5/20), label= "45-55", angle=90, size=4)
ggsave(paste0(filepath, "results/presentation/", "sojourn_pres.pdf"),
       sojourn_pres, device="pdf",
       width=12, height=5)

sojourn_pres
