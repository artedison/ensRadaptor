#! /bin/bash
nreplicate=N_REPLICATE
ml iomkl/2015.02
cd ./1
mpif90 -O2 -o ens ens.f90
cd ../
for((i=1;i<=${nreplicate};i++))
do
  echo ${i}
  cd ${i}
  if [ ${i} -ne 1 ]; then
    cp ../1/ens ./
  fi
  qsub ens.sh
  cd ..
done
wait
