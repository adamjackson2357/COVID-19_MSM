# rm(list=ls())

suppressPackageStartupMessages({
  library(R.utils)
  library(data.table)
  library(grDevices)
  library(colorspace)
  library(MASS)
  library(scales)
  library(pheatmap)
  library(RColorBrewer)
})

##### Arguments

args=commandArgs(trailingOnly=TRUE, asValues=TRUE)
execpath=as.character(args$execpath)
filepath=as.character(args$filepath)
name=as.character(args$name)
seed=as.character(args$seed)
n_params=as.numeric(args$n_params)
n_trans=as.numeric(args$n_trans)
NIter=as.character(format(args$n_iter, scientific=F))
BI=as.character(args$burn_in)
type=as.character(args$type)

# execpath="/rds/general/user/aj1520/home/summer_project/covid_msm/visualisations/"
# filepath="/rds/general/user/aj1520/home/summer_project/covid_msm/"
# name="iwrd"
# n_params=52
# n_trans=6
# NIter=as.character(format(1000000, scientific = F))
# BI=400000
# type="time"
NIter1=613800
NIter2=386200

seed = c(1,2,3,4,5)
# seed = 1

model_path=paste0(filepath, "results/mcmc/", name, "_", n_params, "_", n_trans, "_", NIter, "_", type)
model_path2=paste0(filepath, "results/mcmc/", name, "_", n_params, "_", n_trans, "_", as.character(format(NIter2, scientific = F)), "_", type)

source(paste0(execpath, "/monitoring_functions.R"))

print(paste("File path: ", filepath))
print(paste("Model Folder:", model_path))
print(paste("Name:", name))
print(paste("seed:", seed))
print(paste("Number of params:", n_params))
print(paste("Number of trans:", n_trans))
print(paste("Number of iterations:", NIter))
print(paste("Burn in:", BI))

legCol <- c('darkgreen', 'orange', 'darkred', 'skyblue', "goldenrod2")

transition_names = list("_0" = " O-I",
                        "_1" = " O-R",
                        "_2" = " I-D",
                        "_3" = " I-R")

variable_names = list("mu_0" = "Mu",
                      "mu_1" = "Mu",
                      "mu_2" = "Mu",
                      "mu_3" = "Mu",
                  "lambda_1_" = "Age (years)",
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

##### Loading QOMiC outputs

mymodels=NULL
for (k in seed){

  print(k)

  # tmp <- read.table(paste0(model_path, "_", k, "/history.txt"),header=T,stringsAsFactors=F)
  cols <- (fread(paste0(model_path, "_", k, "/history.txt"), header=TRUE, nrows=0, select=seq(1, (n_params+2), 1)))
  tmp1 <- fread(paste0(model_path, "_", k, "/history.txt"),header=FALSE, nrows=NIter1, select=seq(1, (n_params+2), 1))
  tmp2 <- fread(paste0(model_path2, "_", k, "_2/history.txt"),header=FALSE, nrows=NIter2, select=seq(1, (n_params+2), 1))
  tmp = rbind(tmp1, tmp2)
  colnames(tmp) = names(cols)
  tmp$Iter = seq(1, (NIter1+NIter2), 1)
  
  print(dim(tmp))

  # Reassign columns names
  cols = colnames(tmp)
  iter_cols = cols[1:2]
  param_cols = cols[3:(n_params+2)]
  # sigma_cols = cols[(n_params+3):length(cols)]
  
  for(param in names(variable_names)){param_cols=ifelse(grepl(param, param_cols),
                                                        paste0(variable_names[[param]], substr(param_cols,nchar(param_cols)-1,nchar(param_cols))),
                                                        param_cols)}
  for(param in names(transition_names)){param_cols=ifelse(grepl(param, param_cols),
                                                          paste0(substr(param_cols,1,nchar(param_cols)-2), transition_names[[param]]),
                                                          param_cols)}
  colnames(tmp) = c(iter_cols, param_cols)
  
  assign(paste0("Res",k), tmp)
  mymodels=c(mymodels, paste0("Res",k))
}

##### Path to outputs
figure_path=paste0(model_path,"/figures/")
dir.create(paste0(figure_path), showWarnings=FALSE)
print(paste("Output path:", figure_path))
param_names = colnames(Res1)[2:(n_params+2)]
param_names = paste0("`", param_names, "`")

##### Parameters and likelihood for the presentation

{pdf("/rds/general/user/aj1520/home/summer_project/covid_msm/results/presentation/diagnostics_pres.pdf",height=5,width=24)
  par(mfrow=c(1,3), mar=c(5,5,2,2))
  # layout(matrix(c(1:6, rep(7,3)), ncol=3, byrow=TRUE))
  
  i=3
  MonitorParams(models=mymodels, param=param_names[i], name=param_names[i], init=1, legend=FALSE, sigma=FALSE, factor=100)
  MonitorParams(models=mymodels, param=param_names[i], name=param_names[i], init=BI, legend=FALSE, sigma=FALSE, factor=100)
  GetPosterior(models=mymodels, param=param_names[i], name=param_names[i], init=BI, legend=FALSE)

  dev.off()
}

##### Parameters and likelihood over MCMC iterations

# {pdf(paste0(figure_path, 'monitor_params.pdf'),height=15,width=20)
#   par(mfrow=c(3,2), mar=c(5,5,2,2))
#   # layout(matrix(c(1:6, rep(7,3)), ncol=3, byrow=TRUE))
# 
#   for(i in 2:length(param_names)){
#     MonitorParams(models=mymodels, param=param_names[i], name=param_names[i], init=1, legend=TRUE, sigma=FALSE, factor=100)
#   }
# 
#   MonitorParams(models=mymodels, param="L", name="L", init=1, legend=FALSE, sigma=FALSE, factor=100)
#   dev.off()
# }
# 
# {pdf(paste0(figure_path, 'monitor_params_after_bi.pdf'),height=15,width=20, )
#   par(mfrow=c(3,3), mar=c(5,5,2,2))
#   # layout(matrix(c(1:6, rep(7,3)), ncol=3, byrow=TRUE))
# 
#   for(i in 2:length(param_names)){
#     MonitorParams(models=mymodels, param=param_names[i], name=param_names[i], init=BI, legend=TRUE, sigma=FALSE, factor=100)
#   }
# 
#   MonitorParams(models=mymodels, param="L", name="L", init=BI, legend=FALSE, sigma=FALSE, factor=100)
#   dev.off()
# }
# 
# ##### Posterior distribution (starting from BI)
# 
# {pdf(paste0(figure_path, 'posterior_iter.pdf'),height=15,width=20)
#   par(mfrow=c(3,3), mar=c(5,5,2,2))
#   # layout(matrix(c(1:6, rep(7,3)), ncol=3, byrow=TRUE))
#   
#   for(i in 2:length(param_names)){
#     GetPosterior(models=mymodels, param=param_names[i], name=param_names[i], init=BI, legend=TRUE)
#   }
# 
#   GetPosterior(models=mymodels, param="L", name="L", init=BI, legend=TRUE)
#   
#   dev.off()
# }

##### Likelihood as a function of parameters 

# {pdf(paste0(figure_path, 'monitor_likelihood.pdf'),
#      height=10,width=20)
#   par(mfrow=c(2,3), mar=c(5,5,2,2))
#   
#   for(i in 2:length(param_names)){
#     MonitorLikelihood(models=mymodels, param=param_names[i], name=param_names[i], init=BI, legend=TRUE)
#   }
# 
#   dev.off()
# }
# 
# 
# ##### Correlation between parameters
# 
# Res1=eval(parse(text=mymodels[1]))
# dim(Res1)
# param_names=names(Res1)[3:(2+n_params)]
# param_matrix = vector()
# for (i in 1:length(param_names)){
#   for (j in 1:length(param_names)){
#     if (i==j){next}
#     else{
#       param_matrix = rbind(param_matrix,c(param_names[i], param_names[j]))
#     }
#   }
# }
# 
# loglik_cuts=as.numeric(as.character(cut(Res1$L[BI:nrow(Res1)], breaks = 100, labels = 1:100)))
# loglik_col=colorRampPalette(c("gold","darkred"))(100)
# 
# {pdf(paste0(figure_path, 'scatter_plot_params_after_bi.pdf'),
#      height=15,width=20)
#   par(mfrow=c(3,3), mar=c(5,5,2,2))
# 
#   for (k in 1:nrow(param_matrix)){
#     mybandwidth=0.3
#     myngrid=200
#     mylwd=1
#     x=eval(parse(text=paste0("Res1$",param_matrix[k,1],"[",BI,":",nrow(Res1),"]")))
#     y=eval(parse(text=paste0("Res1$",param_matrix[k,2],"[",BI,":",nrow(Res1),"]")))
#     z <- kde2d(x,y,
#                n=myngrid)
#     contour(z, las=1, main="", cex.main=2,
#             xlab=eval(parse(text=paste0("expression(",param_matrix[k,1],")"))),
#             ylab=eval(parse(text=paste0("expression(",param_matrix[k,2],")"))),
#             cex.axis=1.5, cex.lab=1.5,
#             lwd=mylwd, cex.lab=1.5, col=darken(legCol[1],amount=0.5),
#             pch=19, cex=0.5, las=1, lty=1, drawlabels=FALSE)
#     abline(h=axTicks(2),lty=3,col="grey",lwd=1.5)
#     abline(v=axTicks(1),lty=3,col="grey",lwd=1.5)
#     points(x,y,
#            #col=alpha(lighten(legCol[1], amount=0.0001),0.4),
#            col=alpha(loglik_col[loglik_cuts],0.4),
#            pch=16, cex=0.5)
#     contour(z, las=1, main="", cex.main=2,
#             lwd=mylwd, cex.lab=1.5, col=darken(legCol[1],amount=0.5),
#             pch=19, cex=0.5, las=1, lty=1, drawlabels=FALSE, add=TRUE)
#     #mycor=formatC(cor(x,y, method="spearman"), format = "f", digits=2)
#     mycor=formatC(cor(x,y), format = "f", digits=2)
#     legend("topleft", legend=eval(parse(text=paste0("expression(rho*'=",mycor,"')"))),
#            bty="n", cex=2)
#   }
#   dev.off()
# }
# 
# 
# ##### Correlation of the proposal
# 
# mycov=matrix(as.numeric(Res1[nrow(Res1),grep("sigma", colnames(Res1))]),
#   ncol=sqrt(length(grep("sigma", colnames(Res1)))), byrow=TRUE)
# 
# pheatmap(mycov, cluster_rows=FALSE, cluster_cols=FALSE, border=NA, display_numbers=TRUE, fontsize=20, number_format = "%.2e",
#   filename=paste0(figure_path, 'heatmap_proposal_cov.png'))
# 
# pheatmap(cov2cor(mycov), cluster_rows=FALSE, cluster_cols=FALSE, border=NA, display_numbers=TRUE, fontsize=20,
#   filename=paste0(figure_path, 'heatmap_proposal_cor.png'))


