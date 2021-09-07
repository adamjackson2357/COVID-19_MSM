#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -N iwrd_demo

# navigate to correct directory
cd /rds/general/user/aj1520/home/summer_project/covid_msm/

# load anaconda
module load anaconda3/personal
source activate msm

# set run parameters
name="oird" #run name should correspond to the name of the folder in the configs directory
n_iter=1000000 #number of iterations for the MCMC
burn_in=400000 #MCMC burn-in iterations
adaptive=1 #set as 1 for adaptive iterations
iter_adaptive=100 #H: the adaptive memory iterations
update_frequency=100 #U: the update frequency for the adaptive
n_params=12 #number of parameters to include in the model
n_trans=6 #number of transitions in the model
simul_iter=10000 #number of iterations for the simulation
nb_iter=600000 #n_iter - burn_in
type="demo" #model name corresponds to the variables included
seed=1 #set random seed

# set file paths
FilePath="/rds/general/user/aj1520/home/summer_project/covid_msm/"
ExecPath=$FilePath"qomic_"$type"/main/"
SimulExecPath=$FilePath"qomic_"$type"/QOMiC_simul/"
VisualExecPath=$FilePath"visualisations/"

# set input names
inputDates="configs/"$name"/dates.txt"
inputStatus="configs/"$name"/status.txt"
inputCovars="configs/"$name"/covars_bin_scale.txt"
inputTrans="configs/"$name"/trans.txt"
inputTheta="configs/"$name"/theta_"$type"_"$seed".txt"
inputSigma="configs/"$name"/sigma_"$type"_"$seed".txt"

# create relevant directories
mkdir results
mkdir results/mcmc/
mkdir results/probabilities/
MCMCSeed="results/mcmc/${name}_${n_params}_${n_trans}_${n_iter}_${type}_${seed}_2/"
MCMCPooled="results/mcmc/${name}_${n_params}_${n_trans}_${n_iter}_${type}/"
ProbabilitiesFilePath="results/probabilities/${name}_${n_params}_${n_trans}_${n_iter}_${type}/"
mkdir $MCMCSeed
mkdir $MCMCPooled
mkdir $ProbabilitiesFilePath

output1=$(printf "%s%s" $MCMCSeed)
output2=$(printf "%s%s%s" $MCMCSeed "log.txt")
output3=$(printf "%s%s%s" $ProbabilitiesFilePath "log.txt")

echo $output1
echo $output2
echo $output3

# run the MCMC
mydate=$(date)
echo Start time: $mydate

$ExecPath"QOMiC" -dates $FilePath$inputDates -status $FilePath$inputStatus -covars $FilePath$inputCovars \
-trans $FilePath$inputTrans -theta $FilePath$inputTheta -sigma $FilePath$inputSigma \
-seed $seed -iter $n_iter -burn_in $burn_in \
-adaptive $adaptive -iter_adaptive $iter_adaptive -update_frequency $update_frequency \
-n_params $n_params -n_trans $n_trans \
-out $output1 >$output2

mydate=$(date)
echo End time: $mydate

# run the monitoring scripts
Rscript $VisualExecPath"R_monitor.R" -execpath $VisualExecPath -filepath $FilePath \
-name $name -seed $seed -n_params $n_params -n_trans $n_trans -n_iter $n_iter -burn_in $burn_in -type $type

# create mcmc output file for use in the simulation
Rscript $VisualExecPath"create_mcmc.R" -filepath $FilePath -name $name \
-n_params $n_params -n_trans $n_trans -n_iter $n_iter -burn_in $burn_in -type $type

# run the simulation
mydate=$(date)
echo Start time: $mydate

$SimulExecPath"QOMiC_simul" -dates $FilePath$inputDates -status $FilePath$inputStatus -covars $FilePath$inputCovars \
-trans $FilePath$inputTrans -seed $seed -iter $simul_iter -MCMC $MCMCPooled"mcmc.txt" \
-nb_iter $nb_iter -burn_in $burn_in  -n_params $n_params -n_trans $n_trans \
-scenario $ProbabilitiesFilePath > $output3

mydate=$(date)
echo End time: $mydate
