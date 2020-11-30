#PBS -S /bin/bash
#PBS -q batch
#PBS -N testmodel
#PBS -l nodes=1:ppn=3:EDR
#PBS -l walltime=500:00:00
#PBS -l mem=3gb
#PBS -M XXXXX@XXX
#PBS -m abe
cd $PBS_O_WORKDIR
#mpif90 -O2 -o ens ens.f90
ml iomkl/2015.02
time mpirun -np $PBS_NP ./ens  1>> ./model.${PBS_JOBID}.out \
           2>> ./model.${PBS_JOBID}.err
