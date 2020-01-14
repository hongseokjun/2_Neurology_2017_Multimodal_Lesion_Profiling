#!/bin/bash
print_help()
{
echo "
-----------------------------------------------
`basename $0` <DATADIR> <PREFIX> <ID> <LEFT_RIGHT> <NUMMESH> <NUM_INSURF> <NUM_SUBSURF> <SAMPLING_SPACE> [EMMADIR SURFSTATDIR]
-----------------------------------------------
Sample relative intensities from the volume and map them onto the surface
Uses MATLAB11a!

If only 7 arguments are provided, the script will use:
	EMMADIR=<NOELSOFT_DIR>/MatlabTools/hosung/emma/
	SURSTATDIR=/data/noel/noel2/local/matlab/surfstat_chicago/

Tested for 64bit platforms!	
hong.seok.jun@gmail.com November 2012
"
}


echo "*************************************************************************"
echo "run stage:$*"
echo "on: `hostname`" 
echo "at: `date`"
echo "The Path is: ${PATH}"
echo "*************************************************************************"
echo ""
echo ""

#---------------------------------------------------------------------
noel_do_echo "1) Argment check and assignment" 
if [ $# -lt 6 ]
then
    print_help
    exit 1
fi 

# input handling
DATADIR=${1}
DTIDIR=${2}
PREFIX=${3}
ID=${4}
LEFT_RIGHT=${5}
NUMMESH=${6}
NUMINSURF=${7}
NUMSUBSURF=${8}
SAMPLING_SPACE=${9}
EMMADIR=${NOELSOFT_DIR}/MatlabTools/hosung/emma/
SURSTATDIR=/data/noel/noel2/local/matlab/surfstat_chicago/
overwrite=0

if [ $# -gt 9 ]
then
	EMMADIR=$10
fi
if [ $# -gt 10 ]
then
	SURFSTATDIR=$11
fi

#---------------------------------------------------------------------
noel_do_echo "Sample relative intensities from the volume and map them onto surfaces"
echo DTIIndicesIntensitySampler'(' "'"${DATADIR}"'",  "'"${DTIDIR}"'", "'"${PREFIX}"'", "'"${ID}"'", "'"${LEFT_RIGHT}"'", ${NUMMESH}, ${NUMINSURF}, ${NUMSUBSURF}, "'"${SAMPLING_SPACE}"'", {"'"dti"'"}, ${overwrite}, 2 ')'
matlab -nodisplay << EOF > ${DATADIR}/${ID}/logs/${PREFIX}_${ID}_${LEFT_RIGHT}_sample_dti.log


  addpath( '${EMMADIR}' );
  addpath( '${SURFSTATDIR}' ); 
  
  DTIIndicesSampler( '${DATADIR}', ...
  							'${DTIDIR}', ...
		    	    		'${PREFIX}', ...
							'${ID}', ...
							'${LEFT_RIGHT}', ...
							 ${NUMMESH}, ...
							 ${NUMINSURF}, ...
							 ${NUMSUBSURF}, ...
						   '${SAMPLING_SPACE}', ...
							  {'dti'}, ...
							 ${overwrite}, ... 
							 2 ); % 1: SurfStatReadVol2Surf, 2: volume_object_evaluate
EOF
