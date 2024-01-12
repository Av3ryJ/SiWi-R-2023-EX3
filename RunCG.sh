#!/bin/bash -l

#SBATCH --job-name=RunCG
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=20
#SBATCH --time=00:30:00
#SBATCH --mail-user=manuel.hommel@fau.de
#SBATCH --mail-type=END,FAIL
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
module load intelmpi
cd .
make all
mpirun -np 80 ./cg 4096 4096 100 -1