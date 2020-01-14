#!/bin/sh

source CorticalSurfaceBasedFeatureExtractor_init-sge.sh
QBATCH='noel64.q'
print=''
LOGDIR="/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/misc/logs/"
CASES=${1}

for CASE in `cat ${CASES}`; do
	${print} qbatch -q ${QBATCH} -N CSF_${CASE} -logfile ${LOGDIR}/CSF_${CASE}.log /data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/01_analysis/pipeline/IndividualProcessing.sh ${CASE}
done
