#!/bin/sh

source CorticalSurfaceBasedFeatureExtractor_init-sge.sh

GROUP=$1		# CON | FCD | TLE | HET | PMG
PREFIX=$2	# TLE or mcd
CASES=$3		# case_list.txt  
OUTDIR=$4 	# "/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/"
DTIDIR=$5	# "/data/noel/noel6/CTControls-1.2.0_64/dti/"
QBATCH=$6	# noel64.q
  
NUMINSURF=3
NUMSUBSURF=7
NUMMESH=81920           # 327680
SAMPLING_SPACE="native" # tal
UseLightVersion=1       # 1: no deformation of surfaces, 0: deform surfaces
KERNEL=5
PARAMTERIC='quadratic'
ITERATION=10
print=""

LOGDIR="/local_raid/seokjun/01_project/IntracorticalAnalysis/04_logs/"

for CASE in `cat ${CASES}`; do
	#${print} qbatch -q ${QBATCH} -N CSF_${CASE} -logfile ${LOGDIR}/dti_sba_${CASE}.log DTISurfaceBasedFeatureExtractor.sh ${OUTDIR} ${DTIDIR} ${PREFIX} ${CASE} ${NUMMESH} ${NUMINSURF} ${NUMSUBSURF} ${SAMPLING_SPACE}
	DTISurfaceBasedFeatureExtractor.sh ${OUTDIR} ${DTIDIR} ${PREFIX} ${CASE} ${NUMMESH} ${NUMINSURF} ${NUMSUBSURF} ${SAMPLING_SPACE}
done
