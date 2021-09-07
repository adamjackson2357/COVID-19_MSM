rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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

args=commandArgs(trailingOnly=TRUE, asValues=TRUE)
filepath=as.character(args$filepath)
name=as.character(args$name)
n_trans=as.numeric(args$n_trans)
seed=as.numeric(args$seed)
n_iter=as.numeric(args$n_iter)
BI=as.numeric(args$burn_in)

filepath="/rds/general/user/aj1520/home/summer_project/covid_msm/"
name="iwrd"
n_iter = 1000000
BI = 400000
n_trans = 6
seed = 1

variable_names = list("lambda_1_" = "Age",
                      "lambda_2_" = "Sex",
                      "lambda_3_" = "Other",
                      "lambda_4_" = "Black",
                      "lambda_5_" = "Intermediate",
                      "lambda_6_" = "Low",
                      "lambda_7_" = "Flat",
                      "lambda_8_" = "BMI",
                      "lambda_9_" = "Former",
                      "lambda_10_" = "Current",
                      "lambda_11_" = "Former",
                      "lambda_12_" = "Current",
                      "lambda_13_" = "Yes",
                      "lambda_14_" = "Yes",
                      "lambda_15_" = "Yes",
                      "lambda_16_" = "Yes",
                      "lambda_17_" = "Yes",
                      "lambda_18_" = "Yes",
                      "lambda_19_" = "Yes",
                      "lambda_20_" = "Yes",
                      "lambda_21_" = "Yes",
                      "lambda_22_" = "Yes",
                      "lambda_23_" = "Hospitalisations",
                      "lambda_24_" = "Age Interaction")

variables = list("lambda_1_" = "Age (years)",
                 "lambda_2_" = "Sex (ref. Female)",
                 "lambda_3_" = "Ethnicity (ref. White)",
                 "lambda_4_" = "Ethnicity (ref. White)",
                 "lambda_5_" = "Education (ref. High)",
                 "lambda_6_" = "Education (ref. High)",
                 "lambda_7_" = "Accommodation (ref. House)",
                 "lambda_8_" = "Body Mass Index (kg/m2)",
                 "lambda_9_" = "Smoking status\n(ref. Never)",
                 "lambda_10_" = "Smoking status\n(ref. Never)",
                 "lambda_11_" = "Alcohol drinker status\n(ref. Never)",
                 "lambda_12_" = "Alcohol drinker status\n(ref. Never)",
                 "lambda_13_" = "Cancer (ref. No)",
                 "lambda_14_" = "Cardiovascular (ref. No)",
                 "lambda_15_" = "Hypertension (ref. No)",
                 "lambda_16_" = "Diabetes (ref. No)",
                 "lambda_17_" = "Respiratory (ref. No)",
                 "lambda_18_" = "Autoimmune (ref. No)",
                 "lambda_19_" = "ACE inhibitors (ref. No)",
                 "lambda_20_" = "ARB (ref. No)",
                 "lambda_21_" = "Oral steroid (ref. No)",
                 "lambda_22_" = "Statin (ref. No)",
                 "lambda_23_" = "Hospitalisations (t-1)",
                 "lambda_24_" = "Hospitalisations (t-1)")
myvariables_labels=unlist(variables)
myvariables=unlist(variables)
mynames=unlist(variable_names)

model_list = list("demo" = 12,
                  "social" = 18,
                  "health" = 28,
                  "medical" = 48,
                  "time" = 52)

variable_cat=c(rep("Demographics",4),
               rep("Social", 3),
               rep("Lifestyle", 5),
               rep("Medical", 10), 
               rep("Hospitalisations", 2))

myorder=unique(variable_cat)
names(variable_cat)=myvariables

outcomes=c("IW", "WD")

get_posterior = function(filepath, name, n_params, n_trans, n_iter, BI, type, seed){

  # get the model
  model_path=paste0(filepath, "results/mcmc/", name, "_", n_params, "_", n_trans, "_", as.character(format(n_iter, scientific = F)), "_", type)
  # cols = names(fread(paste0(model_path, "/history_pooled.txt"), header=TRUE, nrows=0, select=seq(3, (n_params+2), 1)))
  # print(cols)
  model = fread(paste0(model_path, "/history_pooled.txt"), header=TRUE, select=seq(3, (n_params+2), 1))
  # colnames(model) = cols
  dim(model)

  # get the posterior
  posterior = t(rbind(mean=apply(model, 2, mean),
                      apply(model, 2, quantile, probs=c(0.025, 0.975), na.rm=TRUE)))
  posterior=round(exp(posterior), 2)
  posterior = cbind(Parameters=rownames(posterior), posterior)
  colnames = colnames(posterior)
  row.names(posterior) = NULL
  posterior = data.frame(posterior)
  colnames(posterior) = colnames
  
  posterior=posterior[!grepl("mu", posterior$Parameters),]
  
  variables_col=as.character(posterior$Parameters)
  for(param in names(variables)){variables_col=ifelse(grepl(param, variables_col), variables[[param]], variables_col)}
  names_col=as.character(posterior$Parameters)
  for(param in names(variable_names)){names_col=ifelse(grepl(param, names_col), variable_names[[param]], names_col)}
  posterior = cbind("Transitions" = ifelse(grepl("_0", posterior$Parameters), "IW", "WD"),
                    "variables" = variables_col,
                    "names" = names_col,
                    posterior[,2:4])
  posterior[1:3] = apply(posterior[1:3], 2, as.character)
  posterior[4:6] = apply(posterior[4:6], 2, as.numeric)
  
  return(posterior)
}

get_BIC = function(filepath, name, n_params, n_trans, n_iter, BI, type, seed){

  # get the model
  model_path=paste0(filepath, "results/mcmc/", name, "_", n_params, "_", n_trans, "_", as.character(format(n_iter, scientific = F)), "_", type)
  L = fread(paste0(model_path, "/history_pooled.txt"),header=TRUE, select=2)$L
  n = (n_iter - BI)*5
  L_mean = mean(L)
  BIC = round(n_params*log(n) - (2*L_mean), 1)

  return(BIC)
}

# get the posterior
mytable1=get_posterior(filepath, name, 12, n_trans, n_iter, BI, "demo", seed)
mytable2=get_posterior(filepath, name, 18, n_trans, n_iter, BI, "social", seed)
mytable3=get_posterior(filepath, name, 28, n_trans, n_iter, BI, "health", seed)
mytable4=get_posterior(filepath, name, 48, n_trans, n_iter, BI, "medical", seed)
mytable5=get_posterior(filepath, name, 52, n_trans, n_iter, BI, "time", seed)

# # Get the mean BICs
# BIC_vector = c(get_BIC(filepath, name, 12, n_trans, n_iter, BI, "demo", seed),
#                get_BIC(filepath, name, 18, n_trans, n_iter, BI, "social", seed),
#                get_BIC(filepath, name, 28, n_trans, n_iter, BI, "health", seed),
#                get_BIC(filepath, name, 48, n_trans, n_iter, BI, "medical", seed),
#                get_BIC(filepath, name, 52, n_trans, n_iter, BI, "time", seed))

# join them all
mytable = mytable1 %>%
  full_join(mytable2, by = c("Transitions", "variables", "names")) %>%
  full_join(mytable3, by = c("Transitions", "variables", "names")) %>%
  full_join(mytable4, by = c("Transitions", "variables", "names")) %>%
  full_join(mytable5, by = c("Transitions", "variables", "names"))

# save
write.csv(mytable, paste0(filepath, "/results/model_comparisons/posterior_estimates.csv"))

credible_intervals = mytable[,c("Transitions", "variables", "names")]
colnames(credible_intervals) = c("Transitions", "Variables", "Names")
credible_intervals[,"Model 1"] = paste0(format(mytable[,4],2),
                              " [", format(mytable[,5],2), "-",
                              format(mytable[,6],2), "]")
credible_intervals[,"Model 2"] = paste0(format(mytable[,7],2),
                                        " [", format(mytable[,8],2), "-",
                                        format(mytable[,9],2), "]")
credible_intervals[,"Model 3"] = paste0(format(mytable[,10],2),
                                        " [", format(mytable[,11],2), "-",
                                        format(mytable[,12],2), "]")
credible_intervals[,"Model 4"] = paste0(format(mytable[,13],2),
                                        " [", format(mytable[,14],2), "-",
                                        format(mytable[,15],2), "]")
credible_intervals[,"Model 5"] = paste0(format(mytable[,16],2),
                                        " [", format(mytable[,17],2), "-",
                                        format(mytable[,18],2), "]")
credible_intervals=credible_intervals[order(credible_intervals$Transitions, as.numeric(rownames(credible_intervals))),]
credible_intervals[,4:8] = apply(credible_intervals[,4:8], 2, function(x){ifelse(grepl("NA", x), "", x)})

# save and put in supplementary materials
write.csv(credible_intervals, paste0(filepath, "/results/model_comparisons/credible_intervals.csv"))

tmptable1=mytable[mytable$Transitions=="IW",]
tmptable2=mytable[mytable$Transitions=="WD",]

tmptable1 = tmptable1[,4:ncol(tmptable1)]
tmptable2 = tmptable2[,4:ncol(tmptable2)]

tmptable1=as.matrix(tmptable1)
tmptable2=as.matrix(tmptable2)

mycolours=darken(c("navy", "orange"), amount=0.2)
delta=c(-1, 1)
myamount=c(0,0.2,0.4,0.6,0.7)
myamount=rep(0,5)
l=1

background=TRUE
myspacing=6
xseq=rep(myspacing, length(myvariables))
xseq=cumsum(xseq)
mypch=c(18,15,19,4)
mystly=c(1,2,3,4)
mycex=0.7
background=TRUE
mycolours=darken(c("navy", "tomato"), amount=0.2)

mycat=variable_cat[myvariables]
mycolours_cat=darken(lighten(c("grey30", "goldenrod2", "red", "royalblue", "darkolivegreen"), amount=0.2), amount=0.4)
names(mycolours_cat)=unique(variable_cat)
{pdf(paste0(filepath, "/results/model_comparisons/forest_plot.pdf"), useDingbats=FALSE, width=10, height=14)
  par(mfrow=c(8,1), mar=c(1,1,1,5))
  plot.new()
  for (l in 5:1){
    print(l)
    i=1
    tmptable=eval(parse(text=paste0("tmptable",i)))
    plotCI(x=xseq+delta[i], y=as.numeric(tmptable[,seq(1,ncol(tmptable),by=3)[l]]), 
           li=as.numeric(tmptable[,seq(2,ncol(tmptable),by=3)[l]]), 
           ui=as.numeric(tmptable[,seq(3,ncol(tmptable),by=3)[l]]), 
           ylim=range(c(0,as.numeric(c(max(tmptable1,na.rm=TRUE),max(tmptable2,na.rm=TRUE)))), na.rm=TRUE),
           xlim=c(min(xseq)-length(outcomes)-1, max(xseq)+length(outcomes)+1),
           xaxt="n", yaxt="n", pch=mypch, col=lighten(mycolours[i], amount=myamount[l]), lwd=1, cex=mycex, sfrac=0.001,
           ylab="",xlab="")
    xseqgreysep=c(min(xseq)-myspacing/2,apply(rbind(xseq[-1],xseq[-length(xseq)]),2,mean),max(xseq)+myspacing/2)
    if (background){
      for (k in seq(1,length(xseqgreysep),by=2)){
        # if (!is.na(tmptable[k/2,seq(1,ncol(tmptable),by=3)[l]])){
        polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]), 
                y=c(-10,-10,10,10), col=lighten(mycolours_cat[mycat[k]], amount=0.96), border=NA)
        # } else {
        #   polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]), 
        #           y=c(-10,-10,10,10), col="grey98", border=NA)
        # }
      }
      for (k in seq(2,length(xseqgreysep),by=2)){
        # if (!is.na(tmptable[k,seq(1,ncol(tmptable),by=3)[l]])){
        polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]), 
                y=c(-10,-10,10,10), col=lighten(mycolours_cat[mycat[k]], amount=0.975), border=NA)
        # } else {
        #   polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]), 
        #           y=c(-10,-10,10,10), col="grey98", border=NA)
        # }
      }
      box()
    }
    for (i in 1:2){
      tmptable=eval(parse(text=paste0("tmptable",i)))
      plotCI(x=xseq+delta[i], y=as.numeric(tmptable[,seq(1,ncol(tmptable),by=3)[l]]), 
             li=as.numeric(tmptable[,seq(2,ncol(tmptable),by=3)[l]]), 
             ui=as.numeric(tmptable[,seq(3,ncol(tmptable),by=3)[l]]), 
             ylim=range(as.numeric(c(tmptable)), na.rm=TRUE), 
             xlim=c(min(xseq)-length(outcomes)-1, max(xseq)+length(outcomes)+1),
             xaxt="n", yaxt="n", pch=mypch[i], col=lighten(mycolours[i], amount=myamount[l]), lwd=1, cex=mycex, sfrac=0.001,
             ylab="",xlab="", add=TRUE, slty=mystly[i])
    }
    axis(side=4, at=axTicks(2), cex.axis=0.7)
    mtext(side=4, text=paste("Model",l), line=3, cex.lab=0.6)
    # mtext(side=4, text=paste("BIC: ", BIC_vector[l]), line=4, cex=0.5)
    abline(h=1, lty=2)
    if (l==5){
      xseqblack=c(xseq[!duplicated(variable_cat[myvariables])]-myspacing/2, max(xseq)+myspacing/2)
      abline(v=xseqblack,lty=3,col="black")
      for (k in 1:(length(xseqblack)-1)){
        axis(side=3, at=xseqblack[c(k,k+1)]+c(1,-1), line=0.5, labels=NA, col=unique(mycolours_cat)[k])
      }
      for (k in 1:(length(xseqblack)-1)){
        axis(side=3, at=mean(xseqblack[c(k,k+1)]), line=0.2, tick=FALSE, labels=unique(variable_cat)[k], 
             cex.axis=1.5, col.axis=unique(mycolours_cat)[k])
      }
    }
    abline(v=xseqgreysep,lty=3, col="grey")
    xseqblack=c(xseq[!duplicated(variable_cat[myvariables])]-myspacing/2, max(xseq)+myspacing/2)
    abline(v=xseqblack,lty=3,col="black")
  }
  for (k in 1:length(xseq)){
    mytext=(mynames)[k]
    mytext=ifelse(mytext=="NA", yes=NA, no=mytext)
    if (is.na(mynames[k])){
      # mytext=WrapText(myvariables)[k]
      mytext=eval(parse(text=paste0("expression('",paste(strsplit(myvariables_labels[k], split = "^", fixed = TRUE)[[1]], collapse="'^'"), "')")))
    }
    axis(side=1, at=xseq[k], labels=mytext, las=2, cex.axis=1.1, col.axis=mycolours_cat[mycat[k]])
  }
  xseqgrey=c(xseq[!duplicated(myvariables)]-myspacing/2, max(xseq)+myspacing/2)
  tmpseq=xseqgrey
  for (k in 1:(length(tmpseq)-1)){
    if (!is.na(mynames)[!duplicated(myvariables)][k]){
      axis(side=1, at=tmpseq[c(k,k+1)]+c(2,-2), line=10.5, labels=NA, 
           col=mycolours_cat[variable_cat[myvariables[!duplicated(myvariables)][k]]])
    }
  }
  for (k in 1:(length(tmpseq)-1)){
    if (!is.na(mynames)[!duplicated(myvariables_labels)][k]){
      # mytmp=eval(parse(text=paste0("expression(atop('",paste(strsplit(myvariables[!duplicated(myvariables)][k], split = "^", fixed = TRUE)[[1]], collapse="'^'"), "'))")))
      mytext=paste(strsplit(myvariables_labels[!duplicated(myvariables_labels)][k], split = "^", fixed = TRUE)[[1]], collapse="'^'")
      if (grepl("\n", mytext, fixed=TRUE)){
        mytext=paste0("'",paste(strsplit(mytext, split="\n", fixed=TRUE)[[1]],collapse="','"),"'")
        mytmp=eval(parse(text=paste0("expression(atop(", mytext, "))")))
      } else {
        mytmp=mytext
      }
      axis(side=1, at=mean(tmpseq[c(k,k+1)]), line=10, tick=FALSE, cex.axis=1.1,
           labels=mytmp, las=2, col.axis=mycolours_cat[variable_cat[myvariables[!duplicated(myvariables)][k]]])
    }
  }
  dev.off()
}


{pdf(paste0(filepath, "/results/presentation/forest_plot_presentation.pdf"), useDingbats=FALSE, width=25, height=16)
  par(mfrow=c(8,1), mar=c(1,1,1,5), oma=c(10,0,0,0))
  plot.new()
  for (l in 5:1){
    print(l)
    i=1
    tmptable=eval(parse(text=paste0("tmptable",i)))
    plotCI(x=xseq+delta[i], y=as.numeric(tmptable[,seq(1,ncol(tmptable),by=3)[l]]), 
           li=as.numeric(tmptable[,seq(2,ncol(tmptable),by=3)[l]]), 
           ui=as.numeric(tmptable[,seq(3,ncol(tmptable),by=3)[l]]), 
           ylim=range(c(0,as.numeric(c(max(tmptable1,na.rm=TRUE),max(tmptable2,na.rm=TRUE)))), na.rm=TRUE),
           xlim=c(min(xseq)-length(outcomes)-1, max(xseq)+length(outcomes)+1),
           xaxt="n", yaxt="n", pch=mypch, col=lighten(mycolours[i], amount=myamount[l]), lwd=1, cex=mycex, sfrac=0.001,
           ylab="",xlab="")
    xseqgreysep=c(min(xseq)-myspacing/2,apply(rbind(xseq[-1],xseq[-length(xseq)]),2,mean),max(xseq)+myspacing/2)
    if (background){
      for (k in seq(1,length(xseqgreysep),by=2)){
        # if (!is.na(tmptable[k/2,seq(1,ncol(tmptable),by=3)[l]])){
        polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]), 
                y=c(-10,-10,10,10), col=lighten(mycolours_cat[mycat[k]], amount=0.96), border=NA)
        # } else {
        #   polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]), 
        #           y=c(-10,-10,10,10), col="grey98", border=NA)
        # }
      }
      for (k in seq(2,length(xseqgreysep),by=2)){
        # if (!is.na(tmptable[k,seq(1,ncol(tmptable),by=3)[l]])){
        polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]), 
                y=c(-10,-10,10,10), col=lighten(mycolours_cat[mycat[k]], amount=0.975), border=NA)
        # } else {
        #   polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]), 
        #           y=c(-10,-10,10,10), col="grey98", border=NA)
        # }
      }
      box()
    }
    for (i in 1:2){
      tmptable=eval(parse(text=paste0("tmptable",i)))
      plotCI(x=xseq+delta[i], y=as.numeric(tmptable[,seq(1,ncol(tmptable),by=3)[l]]), 
             li=as.numeric(tmptable[,seq(2,ncol(tmptable),by=3)[l]]), 
             ui=as.numeric(tmptable[,seq(3,ncol(tmptable),by=3)[l]]), 
             ylim=range(as.numeric(c(tmptable)), na.rm=TRUE), 
             xlim=c(min(xseq)-length(outcomes)-1, max(xseq)+length(outcomes)+1),
             xaxt="n", yaxt="n", pch=mypch[i], col=lighten(mycolours[i], amount=myamount[l]), lwd=1, cex=mycex, sfrac=0.001,
             ylab="",xlab="", add=TRUE, slty=mystly[i])
    }
    axis(side=4, at=axTicks(2), cex.axis=1.5)
    # mtext(side=4, text=paste("Model",l), line=3, cex.lab=5)
    # mtext(side=4, text=paste("BIC: ", BIC_vector[l]), line=4, cex=0.5)
    abline(h=1, lty=2)
    if (l==5){
      xseqblack=c(xseq[!duplicated(variable_cat[myvariables])]-myspacing/2, max(xseq)+myspacing/2)
      abline(v=xseqblack,lty=3,col="black")
      for (k in 1:(length(xseqblack)-1)){
        axis(side=3, at=xseqblack[c(k,k+1)]+c(1,-1), line=0.5, labels=NA, col=unique(mycolours_cat)[k])
      }
      for (k in 1:(length(xseqblack)-1)){
        axis(side=3, at=mean(xseqblack[c(k,k+1)]), line=0.2, tick=FALSE, labels=unique(variable_cat)[k], 
             cex.axis=3, col.axis=unique(mycolours_cat)[k])
      }
    }
    abline(v=xseqgreysep,lty=3, col="grey")
    xseqblack=c(xseq[!duplicated(variable_cat[myvariables])]-myspacing/2, max(xseq)+myspacing/2)
    abline(v=xseqblack,lty=3,col="black")
  }
  for (k in 1:length(xseq)){
    mytext=(mynames)[k]
    mytext=ifelse(mytext=="NA", yes=NA, no=mytext)
    if (is.na(mynames[k])){
      # mytext=WrapText(myvariables)[k]
      mytext=eval(parse(text=paste0("expression('",paste(strsplit(myvariables_labels[k], split = "^", fixed = TRUE)[[1]], collapse="'^'"), "')")))
    }
    axis(side=1, at=xseq[k], labels=mytext, las=2, cex.axis=2, col.axis=mycolours_cat[mycat[k]])
  }
  xseqgrey=c(xseq[!duplicated(myvariables)]-myspacing/2, max(xseq)+myspacing/2)
  tmpseq=xseqgrey
  for (k in 1:(length(tmpseq)-1)){
    if (!is.na(mynames)[!duplicated(myvariables)][k]){
      axis(side=1, at=tmpseq[c(k,k+1)]+c(2,-2), line=14, labels=NA, 
           col=mycolours_cat[variable_cat[myvariables[!duplicated(myvariables)][k]]])
    }
  }
  for (k in 1:(length(tmpseq)-1)){
    if (!is.na(mynames)[!duplicated(myvariables_labels)][k]){
      # mytmp=eval(parse(text=paste0("expression(atop('",paste(strsplit(myvariables[!duplicated(myvariables)][k], split = "^", fixed = TRUE)[[1]], collapse="'^'"), "'))")))
      mytext=paste(strsplit(myvariables_labels[!duplicated(myvariables_labels)][k], split = "^", fixed = TRUE)[[1]], collapse="'^'")
      if (grepl("\n", mytext, fixed=TRUE)){
        mytext=paste0("'",paste(strsplit(mytext, split="\n", fixed=TRUE)[[1]],collapse="','"),"'")
        mytmp=eval(parse(text=paste0("expression(atop(", mytext, "))")))
      } else {
        mytmp=mytext
      }
      axis(side=1, at=mean(tmpseq[c(k,k+1)]), line=14, tick=FALSE, cex.axis=2.2,
           labels=mytmp, las=2, col.axis=mycolours_cat[variable_cat[myvariables[!duplicated(myvariables)][k]]])
    }
  }
  dev.off()
}

{pdf(paste0(filepath, "/results/presentation/forest_model_5_presentation.pdf"), useDingbats=FALSE, width=25, height=16)
  par(mfrow=c(8,1), mar=c(1,1,1,5), oma=c(10,0,0,0))
  plot.new()
  for (l in 5){
    print(l)
    i=1
    tmptable=eval(parse(text=paste0("tmptable",i)))
    plotCI(x=xseq+delta[i], y=as.numeric(tmptable[,seq(1,ncol(tmptable),by=3)[l]]), 
           li=as.numeric(tmptable[,seq(2,ncol(tmptable),by=3)[l]]), 
           ui=as.numeric(tmptable[,seq(3,ncol(tmptable),by=3)[l]]), 
           ylim=range(c(0,as.numeric(c(max(tmptable1,na.rm=TRUE),max(tmptable2,na.rm=TRUE)))), na.rm=TRUE),
           xlim=c(min(xseq)-length(outcomes)-1, max(xseq)+length(outcomes)+1),
           xaxt="n", yaxt="n", pch=mypch, col=lighten(mycolours[i], amount=myamount[l]), lwd=1, cex=mycex, sfrac=0.001,
           ylab="",xlab="")
    xseqgreysep=c(min(xseq)-myspacing/2,apply(rbind(xseq[-1],xseq[-length(xseq)]),2,mean),max(xseq)+myspacing/2)
    if (background){
      for (k in seq(1,length(xseqgreysep),by=2)){
        # if (!is.na(tmptable[k/2,seq(1,ncol(tmptable),by=3)[l]])){
        polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]), 
                y=c(-10,-10,10,10), col=lighten(mycolours_cat[mycat[k]], amount=0.96), border=NA)
        # } else {
        #   polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]), 
        #           y=c(-10,-10,10,10), col="grey98", border=NA)
        # }
      }
      for (k in seq(2,length(xseqgreysep),by=2)){
        # if (!is.na(tmptable[k,seq(1,ncol(tmptable),by=3)[l]])){
        polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]), 
                y=c(-10,-10,10,10), col=lighten(mycolours_cat[mycat[k]], amount=0.975), border=NA)
        # } else {
        #   polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]), 
        #           y=c(-10,-10,10,10), col="grey98", border=NA)
        # }
      }
      box()
    }
    for (i in 1:2){
      tmptable=eval(parse(text=paste0("tmptable",i)))
      plotCI(x=xseq+delta[i], y=as.numeric(tmptable[,seq(1,ncol(tmptable),by=3)[l]]), 
             li=as.numeric(tmptable[,seq(2,ncol(tmptable),by=3)[l]]), 
             ui=as.numeric(tmptable[,seq(3,ncol(tmptable),by=3)[l]]), 
             ylim=range(as.numeric(c(tmptable)), na.rm=TRUE), 
             xlim=c(min(xseq)-length(outcomes)-1, max(xseq)+length(outcomes)+1),
             xaxt="n", yaxt="n", pch=mypch[i], col=lighten(mycolours[i], amount=myamount[l]), lwd=1, cex=mycex, sfrac=0.001,
             ylab="",xlab="", add=TRUE, slty=mystly[i])
    }
    axis(side=4, at=axTicks(2), cex.axis=1.5)
    # mtext(side=4, text=paste("Model",l), line=3, cex.lab=5)
    # mtext(side=4, text=paste("BIC: ", BIC_vector[l]), line=4, cex=0.5)
    abline(h=1, lty=2)
    if (l==5){
      xseqblack=c(xseq[!duplicated(variable_cat[myvariables])]-myspacing/2, max(xseq)+myspacing/2)
      abline(v=xseqblack,lty=3,col="black")
      for (k in 1:(length(xseqblack)-1)){
        axis(side=3, at=xseqblack[c(k,k+1)]+c(1,-1), line=0.5, labels=NA, col=unique(mycolours_cat)[k])
      }
      for (k in 1:(length(xseqblack)-1)){
        axis(side=3, at=mean(xseqblack[c(k,k+1)]), line=0.2, tick=FALSE, labels=unique(variable_cat)[k], 
             cex.axis=3, col.axis=unique(mycolours_cat)[k])
      }
    }
    abline(v=xseqgreysep,lty=3, col="grey")
    xseqblack=c(xseq[!duplicated(variable_cat[myvariables])]-myspacing/2, max(xseq)+myspacing/2)
    abline(v=xseqblack,lty=3,col="black")
  }
  for (k in 1:length(xseq)){
    mytext=(mynames)[k]
    mytext=ifelse(mytext=="NA", yes=NA, no=mytext)
    if (is.na(mynames[k])){
      # mytext=WrapText(myvariables)[k]
      mytext=eval(parse(text=paste0("expression('",paste(strsplit(myvariables_labels[k], split = "^", fixed = TRUE)[[1]], collapse="'^'"), "')")))
    }
    axis(side=1, at=xseq[k], labels=mytext, las=2, cex.axis=2, col.axis=mycolours_cat[mycat[k]])
  }
  xseqgrey=c(xseq[!duplicated(myvariables)]-myspacing/2, max(xseq)+myspacing/2)
  tmpseq=xseqgrey
  for (k in 1:(length(tmpseq)-1)){
    if (!is.na(mynames)[!duplicated(myvariables)][k]){
      axis(side=1, at=tmpseq[c(k,k+1)]+c(2,-2), line=14, labels=NA, 
           col=mycolours_cat[variable_cat[myvariables[!duplicated(myvariables)][k]]])
    }
  }
  for (k in 1:(length(tmpseq)-1)){
    if (!is.na(mynames)[!duplicated(myvariables_labels)][k]){
      # mytmp=eval(parse(text=paste0("expression(atop('",paste(strsplit(myvariables[!duplicated(myvariables)][k], split = "^", fixed = TRUE)[[1]], collapse="'^'"), "'))")))
      mytext=paste(strsplit(myvariables_labels[!duplicated(myvariables_labels)][k], split = "^", fixed = TRUE)[[1]], collapse="'^'")
      if (grepl("\n", mytext, fixed=TRUE)){
        mytext=paste0("'",paste(strsplit(mytext, split="\n", fixed=TRUE)[[1]],collapse="','"),"'")
        mytmp=eval(parse(text=paste0("expression(atop(", mytext, "))")))
      } else {
        mytmp=mytext
      }
      axis(side=1, at=mean(tmpseq[c(k,k+1)]), line=14, tick=FALSE, cex.axis=2.2,
           labels=mytmp, las=2, col.axis=mycolours_cat[variable_cat[myvariables[!duplicated(myvariables)][k]]])
    }
  }
  dev.off()
}

# crop locally

# plotname="../results/model_comparisons/forest_plot.pdf"
# system(paste("pdfcrop --margin 10",plotname,plotname))
# 
# plotname="/Users/adamjackson/Downloads/plots/forest_plot.pdf"
# system(paste("pdfcrop --margin 10",plotname,plotname))
