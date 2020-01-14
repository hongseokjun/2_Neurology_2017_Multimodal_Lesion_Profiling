#!/bin/bash

print_help()
{
echo "
-----------------------------------------------
`basename $0` <OUT_DIR> <PREFIX> <ID> <LEFT_RIGHT> <NUM_INSURF> <NUM_SUBSURF> <NUMMESH> <SAMPLING_SPACE> <KERNEL> <PARAMETRIC> <USELIGHT>
-----------------------------------------------

This script extracts intra-cortical and subcortical surface-based morphometric and textural features:
1) relative intensity
2) perpendicular-/tangential-gradient
3) sulcal depth
4) curvature 
5) cortical thickness
and will generate a text file for each feature.

example:
$0 '/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/result/', ...
   'TLE', ...
	'0365_1', ...
	'left', ...
	3, ...
	3, ...	
	81920, ...
	'native', ...
	20, ...
	'hk', ... [linear|qaudratic|hk]	
   1;

Tested for 64bit platforms!

hong.seok.jun@gmail.com
November 2012
"
}

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

#---------------------------------------------------------------------
noel_do_echo "Argument check and assignment" 
if [ $# -lt 11 ]
then
    print_help
    exit 1
fi

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
ITERATION=${11}
USELIGHT=${12}
RESULT_DIR=${OUTDIR}/${ID}/measurement/

#---------------------------------------------------------------------
noel_do_echo "1) Extract relative intensities and map them on surfaces"
noel_do_cmd FeatureExtractor_RelativeIntensity.sh ${OUTDIR} ${PREFIX} ${ID} ${LEFT_RIGHT} ${NUMINSURF} ${NUMSUBSURF} ${NUMMESH} ${SAMPLING_SPACE} ${USELIGHT}

#---------------------------------------------------------------------
noel_do_echo "2) Computes perpendicular gradient using RI values sampled from each surface"
noel_do_cmd FeatureExtractor_PerpendicularGradient.sh ${OUTDIR} ${PREFIX} ${ID} ${LEFT_RIGHT} ${NUMMESH} ${NUMINSURF} ${NUMSUBSURF} ${SAMPLING_SPACE}

#---------------------------------------------------------------------
noel_do_echo "3) Computes tangential gradient using RI values sampled from each surface"
noel_do_cmd FeatureExtractor_TangentialGradient.sh ${OUTDIR} ${PREFIX} ${ID} ${LEFT_RIGHT} ${NUMMESH} ${NUMINSURF} ${NUMSUBSURF} ${SAMPLING_SPACE}

#---------------------------------------------------------------------
noel_do_echo "4) Computes cortical thickness"
noel_do_cmd FeatureExtractor_CorticalThickness.sh ${OUTDIR} ${PREFIX} ${ID} ${LEFT_RIGHT} ${NUMMESH} ${SAMPLING_SPACE}

#---------------------------------------------------------------------
noel_do_echo "5) Computes sulcal depth of each surface"
noel_do_cmd FeatureExtractor_SulcalDepth.sh ${OUTDIR} ${PREFIX} ${ID} ${LEFT_RIGHT} ${NUMMESH} ${NUMINSURF} ${SAMPLING_SPACE} ${ITERATION} 1

#---------------------------------------------------------------------
noel_do_echo "6) Computes sulcal depth of each surface using a different strategy"
noel_do_cmd FeatureExtractor_SulcalDepth.sh ${OUTDIR} ${PREFIX} ${ID} ${LEFT_RIGHT} ${NUMMESH} ${NUMINSURF} ${SAMPLING_SPACE} ${ITERATION} 2

#---------------------------------------------------------------------
noel_do_echo "6) Computes curvature of each surface"
noel_do_cmd FeatureExtractor_Curvature.sh ${OUTDIR} ${PREFIX} ${ID} ${LEFT_RIGHT} ${NUMMESH} ${NUMINSURF} ${SAMPLING_SPACE} ${ITERATION}

#---------------------------------------------------------------------
noel_do_echo "7) Smooth all the features computed using heat kernel smoothing"
noel_do_cmd DiffusionSmoother.sh ${OUTDIR} ${PREFIX} ${ID} ${LEFT_RIGHT} ${NUMINSURF} ${NUMSUBSURF} ${NUMMESH} ${SAMPLING_SPACE} ${KERNEL} ${PARAMETRIC} ${USELIGHT}

#---------------------------------------------------------------------
noel_do_echo "8) Resample smoothed features using surface-based registration for correspondence across subjects"
noel_do_cmd SurfaceBasedResampler.sh ${OUTDIR} ${PREFIX} ${ID} ${LEFT_RIGHT} ${NUMINSURF} ${NUMSUBSURF} ${NUMMESH} ${SAMPLING_SPACE} ${KERNEL} ${PARAMETRIC}

#---------------------------------------------------------------------
noel_do_echo "9) Prepare vtk file formats for visualization"
noel_do_cmd VTKFileGenerator.sh ${OUTDIR} ${PREFIX} ${ID} ${SAMPLING_SPACE} ${LEFT_RIGHT} ${NUMMESH} ${NUMINSURF} ${NUMSUBSURF} ${USELIGHT} ${KERNEL} ${PARAMETRIC} 'nosurfreg' 'sm'

echo "***************************************************************************"
echo "end stage:$0 $*" 
echo "at: `date`"
echo "***************************************************************************"
echo ""
echo ""  
