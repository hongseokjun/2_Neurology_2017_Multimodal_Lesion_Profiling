#!/bin/bash
print_help()
{
echo "
-----------------------------------------------
`basename $0` <OUTDIR> <PREFIX> <ID> <LEFT_RIGHT> <NUMMESH> <NUMINSURF> <SAMPLING_SPACE> <ITERATION>
-----------------------------------------------
This script computes sulcal depth for each cortical surface
Tested for 64bit platforms!	

hong.seok.jun@gmail.com November 2012
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
SAMPLING_SPACE=${7}
ITERATION=${8}
MODE=${9}
THRESH=5 # at which depth do we have sulci, e.g. 5 (mm)

if [ ${SAMPLING_SPACE} == 'native' ]; then
	POSTFIX='_native'
elif [ ${SAMPLING_SPACE} == 'tal' ]; then
	POSTFIX=''
fi


if [ ${MODE} == 1 ]; then ## Mode = 1: Use a depth potential function (M Boucher et al 2009)
	#---------------------------------------------------------------------
	noel_do_echo "Calculation of the sulcal depth in the defined space"
	noel_do_echo "1) Perform barycentric smoothing"
	noel_do_cmd_new ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj smooth_surface_barycentric_new.sh ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1.obj ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj ${ITERATION} 

	for (( iter=1; iter <= ${NUMINSURF}; iter++ ))
	do
		noel_do_cmd_new ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj smooth_surface_barycentric_new.sh ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1.obj ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj ${ITERATION} 
	done

	noel_do_cmd_new ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj smooth_surface_barycentric_new.sh ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1.obj ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj ${ITERATION}

	#---------------------------------------------------------------------
	noel_do_echo "2) Compute sulcal depth"
	noel_do_cmd_new ${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_sd.txt /data/noel/noel6/Oct-2010/bin/depth_potential -alpha 0.5 -depth_potential ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj  ${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_sd.txt

	for (( iter=1; iter <= ${NUMINSURF}; iter++ ))
	do
		noel_do_cmd_new ${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_sd.txt /data/noel/noel6/Oct-2010/bin/depth_potential -alpha 0.5 -depth_potential ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj ${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_sd.txt
	done

	noel_do_cmd_new ${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_sd.txt /data/noel/noel6/Oct-2010/bin/depth_potential -alpha 0.5 -depth_potential ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj ${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_sd.txt
elif [ ${MODE} == 2 ]; then ## Mode = 2: Use the conventional method (Overlaying brain mask and calculate the distance from gyral crown) (K Im et al 2008) 
	echo "gunzip"
	gunzip -f ${OUTDIR}/${ID}/final/${PREFIX}_${ID}_t1_final.mnc.gz
	gunzip -f ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_t1_nuc.mnc.gz
	
	if [ ${SAMPLING_SPACE} == 'native' ]; then
		surface_mask2 -binary_mask ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_t1_nuc.mnc ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}_native_t1.obj ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_native_mask.mnc
		
	elif [ ${SAMPLING_SPACE} == 'tal' ]; then
		surface_mask2 -binary_mask ${OUTDIR}/${ID}/final/${PREFIX}_${ID}_t1_final.mnc ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}_t1.obj ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_mask.mnc
	fi
	
	echo "********************************************************"
	echo "Distance transform in brain mask"
	echo		
	mincmorph ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}${POSTFIX}_mask.mnc -successive D ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}${POSTFIX}_mask_D.mnc -clobber
	mincmorph -distance ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}${POSTFIX}_mask_D.mnc ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_distance_${LEFT_RIGHT}${POSTFIX}.mnc -clobber

	echo "********************************************************"
	echo "Intersect with volume"
	echo		
	echo "volume_object_evaluate ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_distance_${LEFT_RIGHT}${POSTFIX}.mnc ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1.obj  ${OUTDIR}/${ID}/temp/${PREFIX}_${ID}_distance_mask_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}.txt"
			volume_object_evaluate ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_distance_${LEFT_RIGHT}${POSTFIX}.mnc ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1.obj  ${OUTDIR}/${ID}/temp/${PREFIX}_${ID}_distance_mask_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}.txt
			
	
	echo "********************************************************"
	echo "Generate Geodesic depth map using Matlab"
	echo 

	matlab -nodisplay <<EOF > ${OUTDIR}/${ID}/logs/${PREFIX}_${ID}_${LEFT_RIGHT}_${NUMMESH}_sulcaldepth.log
		addpath( '${EMMADIR}' );
		addpath( '${SURFSTATDIR}' ); 
		
		%addpath( genpath('/data/noel/noel6/seokjun/noelsoft/SurfaceAnalysis/Sulcation_Depth/') );
		mask       = SurfStatReadData({'${OUTDIR}/${ID}/temp/${PREFIX}_${ID}_distance_mask_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}.txt'})
		files      = SurfStatListDir(['${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${postfix}_t1.obj']);
		surf       = SurfStatReadSurf( files );
		surf.tri = transpose(surf.tri); 
		mask2 = mask < ${THRESH}; 
		which surfGeoDist.m
		map = surfGeoDist(surf, mask2)
		fid = fopen('${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_sd2.txt','wt');
		fprintf(fid,'%f\n',map);
		fclose(fid);
EOF

fi
