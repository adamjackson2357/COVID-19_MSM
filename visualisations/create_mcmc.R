# rm(list=ls())

suppressPackageStartupMessages({
  library(R.utils)
  library(data.table)
})

write_inputs <- function(df, df_name, output_path){
  ## Save simulated dates
  write.table(matrix(c(nrow(df), ncol(df)), ncol=1), paste0(output_path, df_name, "_index.txt"), 
              quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(df, paste0(output_path, "pre_", df_name, ".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
  system(paste0("cat ",output_path, df_name, "_index.txt ", output_path, "pre_", df_name, ".txt > ", output_path, df_name, ".txt"))
  system(paste0("rm ", output_path, "pre_", df_name, ".txt"))
  system(paste0("rm ", output_path, df_name, "_index.txt"))
}

## Parameters
args=commandArgs(trailingOnly=TRUE, asValues=TRUE)
execpath=as.character(args$execpath)
filepath=as.character(args$filepath)
name=as.character(args$name)
n_params=as.numeric(args$n_params)
n_trans=as.numeric(args$n_trans)
n_iter=as.numeric(args$n_iter)
BI=as.numeric(args$burn_in)
type=as.character(args$type)

# execpath="/rds/general/user/aj1520/home/summer_project/covid_msm/visualisations/"
# filepath="/rds/general/user/aj1520/home/summer_project/covid_msm/"
# name="iwrd"
# n_params=52
# n_trans=6
# n_iter=400000
# BI=200000
# type="time"
NIter1=613800
NIter2=386200

seed = c(1)

model_path=paste0(filepath, "results/mcmc/", name, "_", n_params, "_", n_trans, "_", as.character(format(n_iter, scientific = F)), "_", type)
# model_path2=paste0(filepath, "results/mcmc/", name, "_", n_params, "_", n_trans, "_", NIter2, "_", type)

print(paste("File path: ", filepath))
print(paste("Model Folder:", model_path))
print(paste("Name:", name))
print(paste("seed:", seed))
print(paste("Number of params:", n_params))
print(paste("Number of trans:", n_trans))
print(paste("Number of iterations:", n_iter))
print(paste("Burn in:", BI))

cols <- (fread(paste0(model_path, "_", 1, "/history.txt"), header=TRUE, nrows=0, select=seq(1, (n_params+2), 1)))

mcmc=NULL
for (k in seed){
  print(k)
  # tmp <- fread(paste0(model_path, "_", k, "/history.txt"),header=FALSE, skip=(BI+1),select=seq(1, (n_params+2), 1))
  tmp <- fread(paste0(model_path, "_", k, "/history.txt"),header=FALSE, skip=(BI+1), nrows=(n_iter-BI), select=seq(1, (n_params+2), 1))
  # tmp2 <- fread(paste0(model_path2, "_", k, "_2/history.txt"),header=FALSE, nrows=NIter2, select=seq(1, (n_params+2), 1))
  # tmp = rbind(tmp1, tmp2)
  print(dim(tmp))

  mcmc = rbind(mcmc, tmp)
}

colnames(mcmc) = names(cols)
print(dim(mcmc))

# save
print(paste0("Resaving history file to: ", model_path, "/"))
write_inputs(mcmc, "mcmc", paste0(model_path, "/"))
write.table(mcmc, paste0(model_path, "/history_pooled.txt"), row.names = FALSE)
