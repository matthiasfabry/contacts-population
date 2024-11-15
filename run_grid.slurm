#!/bin/bash -l
#SBATCH --ntasks=40
#SBATCH --nodes=5
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2g
#SBATCH --clusters=genius
#SBATCH --time=72:00:00
#SBATCH -A lppisn

cd "$SLURM_SUBMIT_DIR" || exit

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MESASDK_ROOT=$VSC_DATA/mesasdk_22.6.1
source "$MESASDK_ROOT"/bin/mesasdk_init.sh
export MESA_DIR=$VSC_DATA/mesa-r22.11.1

readarray -t dirs < <(find "${ARG1}" -mindepth 4 -maxdepth 4 -type d)

cd "${dirs["$PBS_ARRAYID" - 1]}" || exit
echo "$PWD"
./run_both


