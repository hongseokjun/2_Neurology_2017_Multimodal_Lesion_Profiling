#!/bin/bash
source dtiProcess_transmantle_init-sge.sh

# ------------------------------------------------------------------------------
# Are we on 64 bit?
OS=`uname -m`
if [ "$OS" != "x86_64" ]
then
	echo "Not on a 64-bit architecture. Abort"
	exit
fi 

CASES=${1}
QBATCH=${2}
print=''

PREFIX=mcd
BASEDIR=/data/noel/noel7/CTNLES-3T_CIVET-1.2.1/
OUTDIR=/local_raid/seokjun/01_project/07_DTI_process/NLES/
KEEP=0

for CASE in `cat ${CASES}`; do
	${print} qbatch -q ${QBATCH} -N dreg_${CASE} -logfile dti_nlin_${CASE}_processing.log DTI_T1_nonlinear_coregistration.sh ${CASE} ${PREFIX} ${BASEDIR} ${OUTDIR} ${KEEP}
done
