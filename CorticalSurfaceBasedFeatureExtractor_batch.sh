#!/bin/sh

source CorticalSurfaceBasedFeatureExtractor_init-sge.sh

GROUP=$1		# CON | FCD | TLE | HET | PMG
PREFIX=$2	# TLE or mcd
CASES=$3		# case_list.txt  
BASEDIR=$4	# "/data/noel/noel6/CTControls-1.2.0_64/NoelCIVET/"
OUTDIR=$5 	# "/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/"
QBATCH=$6	# noel64.q

LOGDIR="/local_raid/seokjun/01_project/04_IntracorticalAnalysis/04_logs/"

NUMINSURF=3					# number of intracortical surfaces
NUMSUBSURF=3				# number of subcortical surfaces
NUMMESH=81920           # 327680
SAMPLING_SPACE="native" # tal
UseLightVersion=1       # 1: no deformation of surfaces, 0: deform surfaces
KERNEL=5
PARAMTERIC='quadratic'
ITERATION=10
print=""

# DO NOT EDIT parameters below unless you figure out the process of surface generation.
THRES_ROUGH=0.2
WJt=1.0
WJd=0.9

for CASE in `cat ${CASES}`; do
	#InputFilesPreparation.sh ${GROUP} ${BASEDIR} ${CASE} ${PREFIX} ${OUTDIR} ${NUMMESH}
   #${print} qbatch -q ${QBATCH} -N CSF_${CASE} -logfile ${LOGDIR}/CSF_${CASE}.log CorticalSurfaceBasedFeatureExtractor.sh ${BASEDIR} ${OUTDIR} ${LOGDIR} ${PREFIX} ${CASE} ${NUMMESH} ${NUMINSURF} ${NUMSUBSURF} ${SAMPLING_SPACE} ${UseLightVersion}
	CorticalSurfaceBasedFeatureExtractor.sh ${BASEDIR} ${OUTDIR} ${LOGDIR} ${PREFIX} ${CASE} ${NUMMESH} ${NUMINSURF} ${NUMSUBSURF} ${SAMPLING_SPACE} ${UseLightVersion}
done
