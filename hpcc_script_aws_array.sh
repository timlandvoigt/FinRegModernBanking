#$ -j y
#$ -o HPCC_results
# qsub-aws -N e20201117 -t 1-40 -pe openmp 12 hpcc_script_aws_array.sh


EXPERDEF="`echo "$JOB_NAME" | tr -d -c 0-9`"
echo "Running experiment definition $EXPERDEF with $SGE_TASK_LAST runs"	> output.log

EXPER_FILE=experdef_"$EXPERDEF".m
matlab -nodisplay -nodesktop -r "run('$EXPER_FILE');exit;" 

NECON="`cat experdef_"$EXPERDEF".txt | wc -l`"
ECONOMY="`sed -n "$SGE_TASK_ID"p experdef_"$EXPERDEF".txt`"


PPN=$NSLOTS
MAXITER=300
TOL=.00001

echo "Running experiment $ECONOMY with $NSLOTS workers"	>> output.log

EXPERNAME=$ECONOMY

LOGFILE=HPCC_results/${EXPERNAME}.log 
RES_NAME=HPCC_results/res_"$EXPERDEF"_"$EXPERNAME"

EXPER_FILE=experdef_"$EXPERDEF".m
GUESS_STATUS=no_guess
GUESS_NAME=res_20201117_base

echo "Starting ini0 grid run $EXPERNAME for $MAXITER iterations" > $LOGFILE

matlab -nodisplay -nodesktop -r "experdef_file='$EXPER_FILE'; expername='$EXPERNAME'; guess_mode='$GUESS_STATUS'; guess_path='$GUESS_NAME';run('mainSB_create_env.m');exit;" >> $LOGFILE


matlab -nodisplay -nodesktop -r "exper_path='env_$EXPERNAME';maxit=$MAXITER;no_par_processes=$PPN;tol_avg=$TOL;outname='$RES_NAME';run main_run_exper.m;exit;" >> $LOGFILE


