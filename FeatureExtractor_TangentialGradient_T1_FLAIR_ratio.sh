#!/bin/bash

print_help()
{
echo "
-----------------------------------------------
`basename $0`<OUT_DIR> <PREFIX> <ID> <LEFT_RIGHT> <NUMMESH> <NUM_INSURF> <NUM_SUBSURF> <SAMPLING_SPACE>
-----------------------------------------------

This script computes tangential gradient using RI values sampled from each surface.
The gradient is calculated as following: Grad = (sum(abs(RI(i)-RI(x))/d))/# of neighbourings;
	  
example:
$0 ...
	'/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/result/', ...
	'TLE', ...
	'0365_1', ...
	'left', ...
	 81920, ...
	 3, ...
	 3, ...	
	'native';

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

#---------------------------------------------------------------------
noel_do_echo "Argument check and assignment" 
if [ $# -lt 8 ]
then
    print_help
    exit 1
fi
 
OUTDIR=${1}
PREFIX=${2}
ID=${3}
LEFT_RIGHT=${4}
NUMMESH=${5}
NUMINSURF=${6}
NUMSUBSURF=${7}
SAMPLING_SPACE=${8}
PVE_CORR=0

EMMADIR=${NOELSOFT_DIR}/MatlabTools/hosung/emma/
SURSTATDIR=/data/noel/noel2/local/matlab/surfstat_chicago/

#---------------------------------------------------------------------
noel_do_echo "Compute perpendicular gradient"
echo TangentialGradientCalculator_T1_FLAIR_ratio '(' "'"${OUTDIR}"'", "'"${PREFIX}"'", "'"${ID}"'", "'"${LEFT_RIGHT}"'", ${NUMMESH}, ${NUMINSURF}, ${NUMSUBSURF}, "'"${SAMPLING_SPACE}"'", ${PVE_CORR}, {"'"t1_flair_ratio"'"} ')'
matlab -nodisplay << EOF > ${OUTDIR}/${ID}/logs/${PREFIX}_${ID}_${LEFT_RIGHT}_${NUMMESH}_tangential_gradient_ratio.log


  addpath( '${EMMADIR}' );
  addpath( '${SURFSTATDIR}' ); 
  
  TangentialGradientCalculator_T1_FLAIR_ratio( '${OUTDIR}', ...
		    					 							  '${PREFIX}', ...
		 													  '${ID}', ...
								      					  '${LEFT_RIGHT}', ...
  								         					${NUMMESH}, ...
								         					${NUMINSURF}, ...
								         					${NUMSUBSURF}, ...
								      					  '${SAMPLING_SPACE}', ...
										   					${PVE_CORR}, ...
																{'t1_flair_ratio'} );
EOF
