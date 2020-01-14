#!/bin/sh

source CorticalSurfaceBasedFeatureExtractor_init-sge.sh

OUTDIR=${1}
PREFIX=${2}
CASES=${3}
NUMINSURF=${4}
NUMSUBSURF=${5}
NUMMESH=${6}
SAMPLING_SPACE=${7}
KERNEL=${8}
PARAMETRIC=${9}
USELIGHT=${10}
QBATCH='noel64.q'
RESULT_DIR=${OUTDIR}/${ID}/measurement/  

print=''
LOGDIR="/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/misc/logs/"

for CASE in `cat ${CASES}`; do

	${print} qbatch -q ${QBATCH} -N SRL1_${CASE} -logfile ${LOGDIR}/SRL_${CASE}_${KERNEL}.log IndividualProcessing_smoothing.sh ${OUTDIR} ${PREFIX} ${CASE} left ${NUMINSURF} ${NUMSUBSURF} ${NUMMESH} ${SAMPLING_SPACE} 5 ${PARAMETRIC} ${USELIGHT}	
	${print} qbatch -q ${QBATCH} -N SRR1_${CASE} -logfile ${LOGDIR}/SRR_${CASE}_${KERNEL}.log IndividualProcessing_smoothing.sh ${OUTDIR} ${PREFIX} ${CASE} right ${NUMINSURF} ${NUMSUBSURF} ${NUMMESH} ${SAMPLING_SPACE} 5 ${PARAMETRIC} ${USELIGHT}		
done
