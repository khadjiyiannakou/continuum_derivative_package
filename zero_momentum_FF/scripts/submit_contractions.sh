#!/bin/bash

## Author: Konstantin Ottnad
## Last edited: 2014-10-06

## TODO:

## this script acts on a single gauge configuration conf.<NNNN>, where <NNNN> is passed as arg ${1}
## it then submits the job scripts for the 2pt and 3pt FT's / derivative correlators
## finally it moves the resulting correlation functions to ./<NNNN>/corrs/ and cleans up

#####################################################################################################
##                    INPUT PARAMS -- TO BE EDITED ON EVERY ENSEMBLE !!! --
#####################################################################################################
T=64
L=32
RMIN=$L   ## minimum radius for FT
RMAX=$L  ## maximum radius for FT
RSTEP=2  ## stepsize
WALLTIME2="04:00:00"                                    ## walltime for 3pt fourier trafos
WALLTIME3="03:00:00"                                    ## walltime for 2pt fourier trafos
WALLTIME4="02:00:00"                                    ## walltime for tar-ing corrs
declare -a PROJ=("4" "5" "6" "13")                      ## list of projectors
declare -a TS=("12" "14" "16")                          ## list of timeslices for 3pt functions
NODES=4                                                 ## number of NODES for MPI
CORES=16                                                ## number of CORES
declare -a CPUS=("1" "4" "4" "4")                       ## CPUs for MPI (txyz)



#####################################################################################################
##                                   END OF INPUT PARAMS
#####################################################################################################


## params managed by the script itself
WORKDIR=$PWD                ## path to workdir (dir of the script)
N=`printf %.4d ${1}`        ## store 4-digit gauge conf number
echo
echo "Running on gauge field conf.$N"
declare -a QUARKS=("u" "d") ## list of quark flavors
NP=1                        ## total number of cores
for i in ${CPUS[@]}
do
  NP=$(($NP * $i))          ## calculate the total number of cores
done

## perform some checks ...
if [ $NP != $(($NODES*$CORES)) ]; then  ## check if matching number of nodes, cores and MPI-CPUs was given
 echo "FATAL: (NODES=$NODES) * (CORES=$CORES) does not match calculated total number of processors NP=$NP ... "
 exit
fi

## switch to run dir
pushd $N

## create and submit job chain
## do the FT's for 3pt functions (they are extremely optimized and fast, they can all be run sequentially on a single core)
echo "#!/bin/bash -x" > job.2.${N}  ## overwrite !!
echo "#MSUB -l nodes=1:ppn=1" >> job.2.${N}
echo "#MSUB -l walltime=${WALLTIME2}" >> job.2.${N}
echo "#MSUB -e ${WORKDIR}/${N}/error2.${N}.txt" >> job.2.${N}
echo "#MSUB -o ${WORKDIR}/${N}/out2.${N}.txt" >> job.2.${N}
echo "#MSUB -M ottnad@ucy.ac.cy" >> job.2.${N}
echo "#MSUB -m a" >> job.2.${N}  ## send mail only on abort of job
echo "cd \$PBS_O_WORKDIR" >> job.2.${N}
echo "workdir: \$PBS_O_WORKDIR" >> job.2.${N}
echo "declare -a PROJ=(${PROJ[@]})" >> job.2.${N}
echo "declare -a TS=(${TS[@]})" >> job.2.${N}
echo "declare -a QUARKS=(${QUARKS[@]})" >> job.2.${N}
echo "date" >> job.2.${N}
echo "module load gsl" >> job.2.${N}
echo "for p in \${PROJ[@]}; do" >> job.2.${N}
echo "  for q in \${QUARKS[@]}; do" >> job.2.${N}
echo "    for t in \${TS[@]}; do" >> job.2.${N}
echo "      ../Dq_3pt -f Proj\${p}_ts\${t}_\${q}_p_loc.${N} -i block_Proj\${p}_ts\${t}_\${q}_p_loc.${N} -c ${RMIN} -C ${RMAX} -n ${RSTEP} -m ./ini/momentalist" >> job.2.${N}
echo "      ../Dq_3pt -f Proj\${p}_ts\${t}_\${q}_s_loc.${N} -i block_Proj\${p}_ts\${t}_\${q}_s_loc.${N} -c ${RMIN} -C ${RMAX} -n ${RSTEP} -m ./ini/momentalist" >> job.2.${N}
echo "      ../Dq_3pt -f Proj\${p}_ts\${t}_\${q}_v_loc.${N} -i block_Proj\${p}_ts\${t}_\${q}_v_loc.${N} -c ${RMIN} -C ${RMAX} -n ${RSTEP} -m ./ini/momentalist" >> job.2.${N}
echo "      ../Dq_3pt -f Proj\${p}_ts\${t}_\${q}_v_noe.${N} -i block_Proj\${p}_ts\${t}_\${q}_v_noe.${N} -c ${RMIN} -C ${RMAX} -n ${RSTEP} -m ./ini/momentalist" >> job.2.${N}
echo "      ../Dq_3pt -f Proj\${p}_ts\${t}_\${q}_a_loc.${N} -i block_Proj\${p}_ts\${t}_\${q}_a_loc.${N} -c ${RMIN} -C ${RMAX} -n ${RSTEP} -m ./ini/momentalist" >> job.2.${N}
echo "      ../Dq_3pt -f Proj\${p}_ts\${t}_\${q}_vD.${N} -i block_Proj\${p}_ts\${t}_\${q}_vD.${N} -c ${RMIN} -C ${RMAX} -n ${RSTEP} -m ./ini/momentalist" >> job.2.${N}
echo "      ../Dq_3pt -f Proj\${p}_ts\${t}_\${q}_aD.${N} -i block_Proj\${p}_ts\${t}_\${q}_aD.${N} -c ${RMIN} -C ${RMAX} -n ${RSTEP} -m ./ini/momentalist" >> job.2.${N}
echo "      ../Dq_3pt -f Proj\${p}_ts\${t}_\${q}_d1.${N} -i block_Proj\${p}_ts\${t}_\${q}_d1.${N} -c ${RMIN} -C ${RMAX} -n ${RSTEP} -m ./ini/momentalist" >> job.2.${N}
echo "      date" >> job.2.${N}
echo "    done" >> job.2.${N}
echo "  done" >> job.2.${N}
echo "done" >> job.2.${N}
echo "date" >> job.2.${N}
chmod +x ./job.2.${N}

## finally do the FT's for the 2pt functions (they are big and slow))
echo "#!/bin/bash -x" > job.3.${N}  ## overwrite !!
echo "#MSUB -l nodes=1:ppn=1" >> job.3.${N}
echo "#MSUB -l walltime=${WALLTIME3}" >> job.3.${N}
echo "#MSUB -e ${WORKDIR}/${N}/error3.${N}.txt" >> job.3.${N}
echo "#MSUB -o ${WORKDIR}/${N}/out3.${N}.txt" >> job.3.${N}
echo "#MSUB -M ottnad@ucy.ac.cy" >> job.3.${N}
echo "#MSUB -m a" >> job.3.${N}  ## send mail only on abort of job
echo "cd \$PBS_O_WORKDIR" >> job.3.${N}
echo "workdir: \$PBS_O_WORKDIR" >> job.3.${N}
echo "date" >> job.3.${N}
echo "module load gsl" >> job.3.${N}
echo "../Dq_2pt -c ${RMIN} -C ${RMAX} -n ${RSTEP} -f nucleon -m ./ini/momentalist -i block_nucleon.${N}.out" >> job.3.${N}
echo "date" >> job.3.${N}
echo "rm conf.${N} ${FWDPROPBASE}*.inverted ${SEQPROPBASE}*.inverted" >> job.3.${N}
echo "mv block_* ./corrs/" >> job.3.${N}
echo "mv Dq_* ./corrs/" >> job.3.${N}
echo "mv 3pt.r* ./corrs/" >> job.3.${N}
echo "mv 2pt.r* ./corrs/" >> job.3.${N}
echo "date" >> job.3.${N}
chmod +x ./job.3.${N}


## tar everything (and possibly gzip it?)
echo "#!/bin/bash -x" > job.4.${N}  ## overwrite !!
echo "#MSUB -l nodes=1:ppn=1" >> job.4.${N}
echo "#MSUB -l walltime=${WALLTIME4}" >> job.4.${N}
echo "#MSUB -e ${WORKDIR}/${N}/error4.${N}.txt" >> job.4.${N}
echo "#MSUB -o ${WORKDIR}/${N}/out4.${N}.txt" >> job.4.${N}
echo "#MSUB -M ottnad@ucy.ac.cy" >> job.4.${N}
echo "#MSUB -m a" >> job.4.${N}  ## send mail only on abort of job
echo "cd \$PBS_O_WORKDIR" >> job.4.${N}
echo "workdir: \$PBS_O_WORKDIR" >> job.4.${N}
echo "date" >> job.4.${N}
echo "../tar_corrs.sh ${N}" >> job.4.${N}
echo "date" >> job.4.${N}
chmod +x ./job.4.${N}


##### SUBMIT JOBS #####
JOBID2=$(msub ./job.2.${N} 2>&1 | grep -v -e '^$' | sed -e 's/\s*//')
JOBID3=$(msub -W depend=afterok:${JOBID2} ./job.3.${N} 2>&1 | grep -v -e '^$' | sed -e 's/\s*//')
msub -W depend=afterok:${JOBID3} ./job.4.${N}

popd
echo
