#!/bin/bash
#SBATCH --account=def-raveeg1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=128000
#SBATCH --time=0-23:59           # time (DD-HH:MM)
#SBATCH --mail-user=<raveeg1@mcmaster.ca>
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
cd ~/scratch/Appendix/larger_spread/z
mpif90 params.f90 mini.f90 grid.f90 instst.f90 fistst.f90 transition.f90 main.f90 -o ./main
mpirun -np 32 ./main
