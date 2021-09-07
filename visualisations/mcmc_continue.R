# For failed runs, pick up the values where they last were.

rm(list=ls())

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
execpath="/rds/general/user/aj1520/home/summer_project/covid_msm/visualisations/"
filepath="/rds/general/user/aj1520/home/summer_project/covid_msm/"
name="iwrd"
n_params=12
n_trans=6
n_iter=1000000
type="demo"

seed = c(1)

model_path=paste0(filepath, "results/mcmc/", name, "_", n_params, "_", n_trans, "_", as.character(format(n_iter, scientific = F)), "_", type)
input_path = '/rds/general/user/aj1520/home/summer_project/covid_msm/configs/iwrd/'

# iter_vec=NULL
# for (k in seed){
#   print(k)
#   iter <- (fread(paste0(model_path, "_", k, "/history.txt"), header=TRUE, select=1))
#   iter = as.numeric(iter[nrow(iter),])
#   print(iter)
#   iter_vec=c(iter_vec, iter)
#   print(iter_vec)
# }

iter=913900
# round down
# iter=round(min(iter_vec), digits = -2)
print(iter)

mcmc=NULL
for (k in seed){
  print(k)
  theta <- fread(paste0(model_path, "_", k, "/history.txt"),header=FALSE, nrow=1, skip=(iter+1), select=seq(3, (n_params+2), 1))
  sigma <- fread(paste0(model_path, "_", k, "/history.txt"),header=FALSE, nrow=1, skip=(iter+1), select=seq((n_params+3), (n_params+2+(n_params*n_params)), 1))
  sigma = matrix(as.vector(sigma[1,]), nrow=length(theta), ncol=length(theta), byrow=TRUE)
  
  # write theta
  write_inputs(as.data.frame(t(theta)), paste0("theta_", type, "_", k), input_path)
  write_inputs(sigma, paste0("sigma_", type, "_", k), input_path)
}

