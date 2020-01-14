#!/bin/bash

print_help()
{
echo "
-----------------------------------------------
`basename $0` <OUT_DIR> <PREFIX> <ID> <LEFT_RIGHT> <NUM_INSURF> <NUM_SUBSURF> <NUMMESH> <SAMPLING_SPACE> <USELIGHT>
-----------------------------------------------

This script does:
1) computes relative intensity (RI) in the volume and 
2) sample these values by intersecting the volume using surfaces. 
3) Subsequently, it maps them onto each surface. 
4) Furthermore, using partial volume effects (pve) map, 
   it adjusts absolute intensity on the given voxel which pve of csf is over 0.5
	so as to make sure that there is no csf pve on the intensity at the end.

We compute the RI in line with the equation of our following model,
   If the tissue is GM
		100 * (AI - GM_peak) / (BG - GM_peak)
	or if the tissue is WM, then
		100 * (AI - WM_peak) / (WM_peak - BG)
	
	where AI: absolute intensity on a given voxel
	      GM_peak: the most high frequently appearing intensity value in the GM tissue
			WM_peak: for the WM tissue
			BG: the intensity of decision boundary between GM and WM
	  
example:
$0 '/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/result/', ...
   'TLE', ...
	'0365_1', ...
	'left', ...
	3, ...
	3, ...	
	81920, ...
	'native', ...
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
if [ $# -lt 9 ]
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
USELIGHT=${9}
VISUALIZATION=1

EMMADIR=${NOELSOFT_DIR}/MatlabTools/hosung/emma/
SURSTATDIR=/data/noel/noel2/local/matlab/surfstat_chicago/

#-----------------------------------------------------------------------
noel_do_echo "1) Compute RI in the volume using our new model above"
noel_do_cmd RelativeIntensityCalculator_T1_FLAIR_ratio.sh -prefix ${PREFIX} -ID ${ID} -dir ${OUTDIR} -space ${SAMPLING_SPACE}

#-----------------------------------------------------------------------
noel_do_echo "2) Sample and map RI values onto each surface"
noel_do_cmd RelativeIntensitySampler_T1_FLAIR_ratio.sh ${OUTDIR} ${PREFIX} ${ID} ${LEFT_RIGHT} ${NUMMESH} ${NUMINSURF} ${NUMSUBSURF} ${SAMPLING_SPACE}
