#!/bin/bash -l
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=2
#SBATCH --clusters=wice
#SBATCH --time=72:00:00
#SBATCH -A lppisn

cd "$SLURM_SUBMIT_DIR" || exit

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MESASDK_ROOT=$VSC_DATA/mesasdk_22.6.1
source "$MESASDK_ROOT"/bin/mesasdk_init.sh
export MESA_DIR=$VSC_DATA/mesa-r22.11.1

readarray -t dirs < <(find "${ARG1}" -mindepth 5 -maxdepth 5 -type d)

cd "${dirs["$PBS_ARRAYID" - 1]}" || exit
echo "$PWD"

function start_or_restart {
  if [ -d "$1/photos" ] && [ "$(ls -A "$1/photos")" ]; then
    echo "has photos! restarting"
    ./../../../../../../execs/ret
  else
    echo "running from start" "$PWD"
    ./../../../../../../execs/rnrt
  fi
}

if ! [ -f "logs.tar.gz" ] ; then
  ./../../../../../../execs/clean
  ./../../../../../../execs/mk
  start_or_restart "$PWD"
  tar -cvzf "logs.tar.gz" "LOGS1/history.data" "LOGS2/history.data" "term.txt"
  rm "term.txt"
  rm "LOGS1/history.data"
  rm "LOGS2/history.data"
  rm "binary_history.data"
  ./../../../../../../execs/clean
else
  echo "$PWD already has logs.tar.gz, skipping run"
fi
echo "done here: $PWD"


