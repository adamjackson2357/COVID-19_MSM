# rm(list=ls())

suppressPackageStartupMessages({
  library(R.utils)
  library(grDevices)
  library(colorspace)
  library(MASS)
  library(scales)
  library(pheatmap)
  library(RColorBrewer)
  library(data.table)
  library(grid)
  library(gridExtra)
})

##### Arguments

args=commandArgs(trailingOnly=TRUE, asValues=TRUE)
filepath=as.character(args$filepath)
name=as.character(args$name)
seed=as.character(args$seed)
n_params=as.numeric(args$n_params)
n_trans=as.numeric(args$n_trans)
n_iter=as.numeric(args$n_iter)
BI=as.numeric(args$burn_in)
type=as.character(args$type)

# filepath="/rds/general/user/aj1520/home/summer_project/covid_msm/"
# results_filepath = paste0(filepath, "results/mcmc/")
# figure_path = paste0(filepath, "results/model_comparisons/")
# name="iwrd"
# seed=1
# n_params=52
# n_trans=6
# n_iter=400000
# BI=200000
# type = "time"

model_path=paste0(filepath, "results/mcmc/", name, "_", n_params, "_", n_trans, "_", as.character(format(n_iter, scientific = F)), "_", type, "_", seed,"/")
source(paste0(filepath, "visualisations/monitoring_functions.R"))

print(paste("File path: ", filepath))
print(paste("Model Folder:", model_path))
print(paste("Name:", name))
print(paste("seed:", seed))
print(paste("Number of params:", n_params))
print(paste("Number of trans:", n_trans))
print(paste("Number of iterations:", n_iter))
print(paste("Burn in:", BI))

transition_names = list("_0" = " O-I",
                        "_1" = " O-R",
                        "_2" = " I-D",
                        "_3" = " I-R")

variables = list("lambda_1_" = "Age (years)",
                      "lambda_2_" = "Male Sex (ref. Female)",
                      "lambda_3_" = "Other Ethnicity (ref. White)",
                      "lambda_4_" = "Black Ethnicity (ref. White)",
                      "lambda_5_" = "Intermediate Education (ref. High)",
                      "lambda_6_" = "Low Education (ref. High)",
                      "lambda_7_" = "Flat (ref. House)",
                      "lambda_8_" = "Body Mass Index (kg/m^2'*')",
                      "lambda_9_" = "Former Smoker (ref. Never)",
                      "lambda_10_" = "Current Smoker (ref. Never)",
                      "lambda_11_" = "Former Alcohol Drinker (ref. Never)",
                      "lambda_12_" = "Current Alcohol Drinker (ref. Never)",
                      "lambda_13_" = "Cancer (ref. No)",
                      "lambda_14_" = "Cardiovascular (ref. No)",
                      "lambda_15_" = "Hypertension (ref. No)",
                      "lambda_16_" = "Diabetes (ref. No)",
                      "lambda_17_" = "Respiratory (ref. No)",
                      "lambda_18_" = "Autoimmune (ref. No)",
                      "lambda_19_" = "ACE inhibitors (ref. No)",
                      "lambda_20_" = "Angiotensin II receptor blocker (ref. No)",
                      "lambda_21_" = "Oral steroid (ref. No)",
                      "lambda_22_" = "Statin (ref. No)",
                      "lambda_23_" = "Hospitalisations (t-1)",
                      "lambda_24_" = "Age Hospitalisations Interaction (t-1)")

myvariables_labels=unlist(variables)
myvariables=unlist(variables)

variable_cat=c(rep("Demographics",4),
               rep("Social", 3),
               rep("Health Risk Factors", 5),
               rep("Medical", 10), 
               rep("Hospitalisations", 2))

##### Loading QOMiC outputs

mymodels=NULL
for (k in seed){
  # for (k in 1){
  print(k)
  tmp <- read.table(paste0(model_path, "history.txt"),header=T,stringsAsFactors=F)
  assign(paste0("Res",k), tmp)
  mymodels=c(mymodels, paste0("Res",k))
}

param_cols <- names(fread(paste0(model_path, "/history.txt"), header=TRUE, nrows=0, select=seq(3, (n_params+2), 1)))
for(param in names(variables)){param_cols=ifelse(grepl(param, param_cols),
                                                      paste0(variables[[param]], substr(param_cols,nchar(param_cols)-1,nchar(param_cols))),
                                                      param_cols)}
for(param in names(transition_names)){param_cols=ifelse(grepl(param, param_cols),
                                                        paste0(substr(param_cols,1,nchar(param_cols)-2), transition_names[[param]]),
                                                        param_cols)}

##### Path to outputs

figure_path=paste0(model_path,"figures/")
dir.create(paste0(figure_path), showWarnings=FALSE)
print(paste("Output path:", figure_path))

##### Correlation of the proposal

mycov=matrix(as.numeric(Res1[nrow(Res1),grep("sigma", colnames(Res1))]),
             ncol=sqrt(length(grep("sigma", colnames(Res1)))), byrow=TRUE)
colnames(mycov) = param_cols
rownames(mycov) = param_cols
param_cols = param_cols[!grepl("mu ", param_cols)]
mycov = mycov[param_cols,param_cols]
mycor=cov2cor(mycov)

variables = c(unlist(variables), unlist(variables))
categories = c(variable_cat, variable_cat)

# Add annotation as described above, and change the name of annotation
annotation <- data.frame(Var1 = factor(categories))
levels(annotation$Var1) = c("Demographics", "Social", "Health Risk Factors", "Medical", "Hospitalisations")
rownames(annotation) <- colnames(mycor) # check out the row names of annotation

# change the color of annotation to what you want: (eg: "navy", "darkgreen")
cat_colours        <- lighten(c("grey30", "goldenrod2", "red", "royalblue", "darkolivegreen"), amount=0.2)
names(cat_colours) <- unique(variable_cat)
anno_colors <- list(Var1 = cat_colours)
mycor
breaks = seq(-1, 1, 0.02)
cor_plot = pheatmap(mycor, cluster_rows=FALSE, cluster_cols=FALSE,
                    labels_row=variables, labels_col = variables,
                    annotation_row = annotation, annotation_col = annotation, annotation_legend = FALSE,
                    annotation_colors = anno_colors, breaks=breaks,
                    cellheight=16, cellwidth = 16,
                    gaps_col=(nrow(mycor)/2), gaps_row = (nrow(mycor)/2), angle_col = "45")

pdf(paste0(figure_path, 'parameters_correlation.pdf'), height=14, width=16.5)
cor_plot
grid.text("Hospitalisation", x=(11.25/42), y=(35/36), gp=gpar(fontsize=20))
grid.text("Death following hospitalisation", x=(24/42), y=(35/36), gp=gpar(fontsize=20))

grid.text("Hospitalisation", x=(2/42), y=(26/36), gp=gpar(fontsize=20), rot=90)
grid.text("Death following hospitalisation", x=(2/42), y=(13/36), gp=gpar(fontsize=20), rot=90)
dev.off()
