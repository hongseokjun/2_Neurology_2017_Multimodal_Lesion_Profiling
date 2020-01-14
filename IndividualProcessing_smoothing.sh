#!/bin/bash

echo "*************************************************************************"
echo "run stage:$0 $*"
echo "on: `hostname`" 
echo "at: `date`"
echo "The Path is: ${PATH}"
echo "*************************************************************************"
echo ""
echo ""
source ${NOELSOFT_DIR}/BashTools/noel_do_cmd
source ${NOELSOFT_DIR}/BashTools/noel_do_echo
source ${NOELSOFT_DIR}/BashTools/noel_do_cmd_new

OUTDIR=${1}
PREFIX=${2}
ID=${3}
LEFT_RIGHT=${4}
NUMINSURF=${5}
NUMSUBSURF=${6}
NUMMESH=${7}
SAMPLING_SPACE=${8}
KERNEL=${9}
PARAMETRIC=${10}
USELIGHT=${11}
RESULT_DIR=${OUTDIR}/${ID}/measurement/

#---------------------------------------------------------------------
noel_do_echo "7) Smooth all the features computed using heat kernel smoothing"
noel_do_cmd DiffusionSmoother.sh ${OUTDIR} ${PREFIX} ${ID} ${LEFT_RIGHT} ${NUMINSURF} ${NUMSUBSURF} ${NUMMESH} ${SAMPLING_SPACE} ${KERNEL} ${PARAMETRIC} ${USELIGHT}

#---------------------------------------------------------------------
noel_do_echo "8) Resample smoothed features using surface-based registration for correspondence across subjects"
noel_do_cmd SurfaceBasedResampler.sh ${OUTDIR} ${PREFIX} ${ID} ${LEFT_RIGHT} ${NUMINSURF} ${NUMSUBSURF} ${NUMMESH} ${SAMPLING_SPACE} ${KERNEL} ${PARAMETRIC}


#if [ ${KERNEL} == 20 ]; then
#	rm -rf ${OUTDIR}/${ID}/surfaces/*${LEFT_RIGHT}*.vtk
#	gunzip -vf ${OUTDIR}/${ID}/surfaces/*${LEFT_RIGHT}*.obj.gz
#	noel_do_cmd VTKFileGenerator.sh ${OUTDIR} ${PREFIX} ${ID} ${SAMPLING_SPACE} ${LEFT_RIGHT} ${NUMMESH} ${NUMINSURF} ${NUMSUBSURF} ${USELIGHT} ${KERNEL} ${PARAMETRIC} 'nosurfreg' 'sm'
#	gzip -vf ${OUTDIR}/${ID}/surfaces/*${LEFT_RIGHT}*.obj
#fi
