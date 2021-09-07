# rm(list=ls())
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

suppressPackageStartupMessages({
  library(R.utils)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(openxlsx)
  library(plotrix)
  library(colorspace)
  library(data.table)
})

# args=commandArgs(trailingOnly=TRUE, asValues=TRUE)
# filepath=as.character(args$filepath)
# name=as.character(args$name)
# n_trans=as.numeric(args$n_trans)
# seeds=as.numeric(args$seeds)
# n_iter=as.numeric(args$n_iter)
# BI=as.numeric(args$burn_in)

filepath="/rds/general/user/aj1520/home/summer_project/covid_msm/"
name="iwrd"
n_iter = 1000000
BI = 400000
n = n_iter-BI
n_trans = 6
seeds = 5
NIter1=613800
NIter2=386200

seed = c(1,2,3,4,5)

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

get_posterior = function(filepath, name, n_params, n_trans, n_iter, type){
  
  # get the model
  model_path=paste0(filepath, "results/mcmc/", name, "_", n_params, "_", n_trans, "_", as.character(format(n_iter, scientific = F)), "_", type)
  posterior = fread(paste0(model_path, "/history_pooled.txt"),header=TRUE, select=seq(3, (n_params+2), 1))
  
  return(posterior)
}

get_gelman_rubin = function(posterior, chains, n){
  theta_mean_mat = NULL
  sigma_mat = NULL
  for(chain in 1:chains){
    # get the start and end index of the chain
    start=1+((chain-1)*n)
    end=(chain*n)
    print(c(start,end))
    posterior_m = posterior[start:end,]
    
    # get the posterior means
    theta_m_mean = apply(posterior_m, 2, mean)
    
    # get the intra-chain variance
    sigma_m = colSums(t(apply(posterior_m, 1, function(x) (x-theta_m_mean)^2))) / (n-1)
    
    # add to the respective vectors
    theta_mean_mat = rbind(theta_mean_mat, theta_m_mean)
    sigma_mat = rbind(sigma_mat, sigma_m)
  }
  
  # calculate the mean of all chains
  theta_mean = apply(theta_mean_mat, 2, mean)
  
  # compute how the individual means vary around the join mean
  B = colSums(t(apply(theta_mean_mat, 1, function(x) (x-theta_mean)^2))) * (n/(chains-1))
  
  # compute the averaged variance of all chains
  W = colMeans(sigma_mat)
  
  # define V, an unbiased estimator of the true variance
  V = (((n-1)/n) * W) + ((chains+1)/(chains*n))*B
  
  # if the chains have converged, then W is also an unbiased estimate of the true variance
  R = sqrt(V/W)
  
  return(R)
}


n_params=52
type="time"
model_path=paste0(filepath, "results/mcmc/", name, "_", n_params, "_", n_trans, "_", as.character(format(n_iter, scientific = F)), "_", type)
model_path2=paste0(filepath, "results/mcmc/", name, "_", n_params, "_", n_trans, "_", as.character(format(NIter2, scientific = F)), "_", type)
cols <- (fread(paste0(model_path, "_", 1, "/history.txt"), header=TRUE, nrows=0, select=seq(3, (n_params+2), 1)))

mytable5=NULL
for (k in seed){
  print(k)
  tmp1 <- fread(paste0(model_path, "_", k, "/history.txt"),header=FALSE, skip=(BI+1), nrows=(NIter1-BI), select=seq(3, (n_params+2), 1))
  tmp2 <- fread(paste0(model_path2, "_", k, "_2/history.txt"),header=FALSE, nrows=NIter2, select=seq(3, (n_params+2), 1))
  tmp = rbind(tmp1, tmp2)
  print(dim(tmp))
  
  mytable5 = rbind(mytable5, tmp)
}
colnames(mytable5) = names(cols)
print(dim(mytable5))


# get the posterior
# mytable1=get_posterior(filepath, name, 12, n_trans, n_iter, "demo")
# mytable2=get_posterior(filepath, name, 18, n_trans, n_iter, "social")
# mytable3=get_posterior(filepath, name, 28, n_trans, n_iter, "health")
# mytable4=get_posterior(filepath, name, 48, n_trans, n_iter, "medical")
# mytable5=get_posterior(filepath, name, 52, n_trans, n_iter, "time")

# get the Gelman-Rubin score for each model
# R_1 = get_gelman_rubin(mytable1, seeds, n)
# R_2 = get_gelman_rubin(mytable2, seeds, n)
# R_3 = get_gelman_rubin(mytable3, seeds, n)
# R_4 = get_gelman_rubin(mytable4, seeds, n)
R_5 = get_gelman_rubin(mytable5, seeds, n)

# convert to dataframes
# R_1=data.frame(parameters=names(R_1), R_1)
# R_2=data.frame(parameters=names(R_2), R_2)
# R_3=data.frame(parameters=names(R_3), R_3)
# R_4=data.frame(parameters=names(R_4), R_4)
R_5=data.frame(parameters=names(R_5), R_5, rank=seq(1,52,1))

# join
# gelman_rubin = R_1 %>%
#   full_join(R_2, by="parameters") %>%
#   full_join(R_3, by="parameters") %>%
#   full_join(R_4, by="parameters") %>%
#   full_join(R_5, by="parameters")
gelman_rubin = R_5
# add transition and parameter names
param_cols=as.character(gelman_rubin$parameters)
for(param in names(variable_names)){param_cols=ifelse(grepl(param, param_cols),
                                                      variable_names[[param]],
                                                      param_cols)}
gelman_rubin$parameters = param_cols
gelman_rubin = gelman_rubin[order(gelman_rubin$rank),]
gelman_rubin = gelman_rubin[,-7]

trans_cols=c(rep("O-I", 25), "O-R", rep("I-D", 25), "IR")
gelman_rubin=cbind("Transition"=trans_cols, gelman_rubin)

gelman_rubin[,3] = round(gelman_rubin[,3], 3)

write.csv(gelman_rubin, paste0(filepath, "results/model_comparisons/gelman_rubin.csv"),na = "")
