#!/bin/bash 

if [[ $# -ne 4 ]] ; then
    echo '*** Forgetting arguments !! Order : "READIN"  "SAVE" "NCOL" "NITER" ' 
    exit 0
fi

args=("$@")

echo "#!/bin/bash" > submit.sh
echo "#SBATCH --job-name=MM_${args[0]}_${args[1]}_N${args[2]}" >> submit.sh
echo "#SBATCH --mem=0" >> submit.sh 
echo "#SBATCH --nodes=1" >> submit.sh
#echo "#SBATCH -p longq" >> submit.sh 
echo "#SBATCH -p defq" >> submit.sh
#echo "#SBATCH -p debugq" >> submit.sh
echo "#SBATCH --output=/home/rjha1/MM/4MMC_closed/4MMC_${args[0]}_${args[1]}_${args[2]}_${args[3]}" >> submit.sh
echo "#SBATCH --time=24:00:00" >> submit.sh 
#echo "#SBATCH --time=150:00:00" >> submit.sh
echo "" >> submit.sh 

echo "set -euxo pipefail" >> submit.sh 
echo "" >> submit.sh 
echo "pwd" >> submit.sh 
echo "echo \"SLURM_JOB_ID=\$SLURM_JOB_ID\"" >> submit.sh 
echo "module load python/3.7" >> submit.sh 
#echo "env OMP_NUM_THREADS=4 python -u 3MMC.py $1 $2 $3 $4 " >> submit.sh
#echo "env OMP_NUM_THREADS=4 python -u IKKT_b.py $1 $2 $3 $4 4" >> submit.sh
echo "env OMP_NUM_THREADS=4 python -u 4MMC.py $1 $2 $3 $4" >> submit.sh
sbatch submit.sh 
rm submit.sh
