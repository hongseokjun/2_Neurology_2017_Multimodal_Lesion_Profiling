#!/bin/sh

source FeatureExtractor_init-sge.sh

CASES=$1
QBATCH=$2

OUTDIR="/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/"
LOGDIR="/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/04_logs/"

PREFIX="TLE"
NUMINSURF=3
NUMSUBSURF=3
NUMMESH=81920           # 327680
SAMPLING_SPACE="native" # tal
UseLightVersion=1       # 1: no deformation of surfaces, 0: deform surfaces
KERNEL=10
PARAMTERIC='hk'
ITERATION=10
print=""

for CASE in `cat ${CASES}`; do 
    ${print} qbatch -q ${QBATCH} -N LFeE_${CASE} -logfile ${LOGDIR}/FeatureExtractor_${CASE}_left_processing.log FeatureExtractor.sh ${OUTDIR} ${PREFIX} ${CASE} left ${NUMINSURF} ${NUMSUBSURF} ${NUMMESH} ${SAMPLING_SPACE} ${KERNEL} ${PARAMTERIC} ${ITERATION} ${UseLightVersion}
    ${print} qbatch -q ${QBATCH} -N RFeE_${CASE} -logfile ${LOGDIR}/FeatureExtractor_${CASE}_right_processing.log FeatureExtractor.sh ${OUTDIR} ${PREFIX} ${CASE} right ${NUMINSURF} ${NUMSUBSURF} ${NUMMESH} ${SAMPLING_SPACE} ${KERNEL} ${PARAMTERIC} ${ITERATION} ${UseLightVersion}
done
