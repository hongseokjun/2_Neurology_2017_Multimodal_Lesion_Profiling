#!/bin/sh

PREFIX=${1}
CASES=${2}
BASEDIR=${3}

for CASE in `cat ${CASES}`; do
	FSL_BBR.sh ${PREFIX} ${CASE} ${BASEDIR}
done

