#!/bin/bash

# InputFilesPreparation.sh CON /data/noel/noel6/CTControls-1.2.0_64/NoelCIVET/ 301_1 TLE /data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/
# 1) Make directories
# 2) Copy native files and reshape all files for making all the header-file-information same between native and final images
# 3) Copy the linear transformation matrix and invert it
# 4) Copy surface files from the original CIVET directory
# 5) Copy classification files from the original CIVET directory
# 6) Copy a mask file from the original CIVET directory

#################################################################
# Files
# native:
#		${PREFIX}_${ID}_t1_nuc.mnc.gz
#		${PREFIX}_${ID}_[flair|ir]_nuc.mnc.gz						
# final:
#		${PREFIX}_${ID}_t1_final.mnc.gz	
#		${PREFIX}_${ID}_[flair|ir]_final.mnc.gz					
# xfm:
#		${PREFIX}_${ID}_t1_tal.xfm
#		${PREFIX}_${ID}_t1_tal_inverted.xfm
#		${PREFIX}_${ID}_t1_to_tal.xfm
#  	${PREFIX}_${ID}_tal_to_t1.xfm
#		${PREFIX}_${ID}_left_surfmap.sm
#		${PREFIX}_${ID}_left_surfmap_327680.sm
#		${PREFIX}_${ID}_right_surfmap.sm
#		${PREFIX}_${ID}_right_surfmap_327680.sm
#  	${PREFIX}_${ID}_[flair|ir]_to_t1.xfm						
#  	${PREFIX}_${ID}_[flair|ir]_to_tal_concat.xfm				
#  	${PREFIX}_${ID}_tal_to_[flair|ir].xfm						
# surfaces:
#		${PREFIX}_${ID}_gray_surface_left_81920_t1.obj
#		${PREFIX}_${ID}_gray_surface_left_327680_t1.obj
#		${PREFIX}_${ID}_white_surface_left_81920_t1.obj
#		${PREFIX}_${ID}_white_surface_left_327680_t1.obj		
#		${PREFIX}_${ID}_gray_surface_right_81920_t1.obj
#		${PREFIX}_${ID}_gray_surface_right_327680_t1.obj
#		${PREFIX}_${ID}_white_surface_right_81920_t1.obj
#		${PREFIX}_${ID}_white_surface_right_327680_t1.obj
# classify:
#		${PREFIX}_${ID}_ClassFinal.mnc.gz
#		${PREFIX}_${ID}_final_classify.mnc.gz
#		${PREFIX}_${ID}_ClassFinal_native.mnc.gz
#		${PREFIX}_${ID}_final_classify_native.mnc.gz
# mask:
#  	${PREFIX}_${ID}_brain_mask.mnc.gz
#  	${PREFIX}_${ID}_brain_mask_native.mnc.gz
#		${PREFIX}_${ID}_brain_mask_native_[flair|ir].mnc.gz	
# temp:
#		${PREFIX}_${ID}_[flair|ir]_reshaped_nuc.imp     		
#		${PREFIX}_${ID}_[flair|ir]_reshaped_nuc.mnc.gz  		
#  Total: 216.5Mb

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

GROUP=${1}
BASEDIR=${2}
ID=${3}
PREFIX=${4}
OUTDIR=${5}
NUMMESH=${6}

case "${GROUP}" in
'CON') NATDIR=Controls
		;;
'FCD') NATDIR=FCD
		;;
'TLE') NATDIR=TLE
		;;
esac

noel_do_echo "2) Copy native files and make the vtk format of file for visualization later"
MULT_MOD=( ir ) # t1-0.6 dti rs

for MOD_INDEX in "${MULT_MOD[@]}"
do
	#BASEDIR_temp=/data/noel/noel6/native3T/native_${NATDIR}
   BASEDIR_temp=/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/02_org_data/
	
	noel_do_echo "   - ${MOD_INDEX} is reshaped ..."
		cp ${BASEDIR_temp}/${PREFIX}_${ID}_${MOD_INDEX}.mnc* ${OUTDIR}/${ID}/native/
		mv ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_nuc.mnc.gz ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_nuc_org.mnc.gz
		gzip -vf ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}.mnc
		xdim=`mincinfo -dimlength xspace ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}.mnc.gz`
		ydim=`mincinfo -dimlength yspace ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}.mnc.gz`
		zdim=`mincinfo -dimlength zspace ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}.mnc.gz`
		xstart=`mincinfo -attvalue xspace:start ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}.mnc.gz`
		ystart=`mincinfo -attvalue yspace:start ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}.mnc.gz`
		zstart=`mincinfo -attvalue zspace:start ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}.mnc.gz`
		xstep=`mincinfo -attvalue xspace:step ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}.mnc.gz`
		ystep=`mincinfo -attvalue yspace:step ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}.mnc.gz`
		zstep=`mincinfo -attvalue zspace:step ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}.mnc.gz`

		if [ $(echo "${xstep} < 0" | bc) -eq 1 -a $(echo "${xstart} > 0" | bc) -eq 1 ]; then
				xtemp2=`echo $xstep*-1 | bc`
				xtemp1=`echo $xstart-$xdim*$xtemp2+$xtemp2 | bc`
		else
				xtemp2=${xstep}
				xtemp1=${xstart}
		fi
		
		
		if [ $(echo "${ystep} < 0" | bc) -eq 1 -a $(echo "${ystart} > 0" | bc) -eq 1 ]; then		
				ytemp2=`echo $ystep*-1 | bc`
				ytemp1=`echo $ystart-$ydim*$ytemp2+$ytemp2 | bc`
		else
				ytemp2=${ystep}
				ytemp1=${ystart}		
		fi		
		
		if [ $(echo "${zstep} < 0" | bc) -eq 1 -a $(echo "${zstart} > 0" | bc) -eq 1 ]; then		
				ztemp2=`echo $zstep*-1 | bc`
				ztemp1=`echo $zstart-$zdim*$ztemp2+$ztemp2 | bc`
		else
				ztemp2=${zstep}
				ztemp1=${zstart}		
		fi		

		noel_do_cmd_new ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_reshaped.mnc.gz mincresample -xstart $xtemp1 -ystart $ytemp1 -zstart $ztemp1 -xstep $xtemp2 -ystep $ytemp2 -zstep $ztemp2 ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}.mnc.gz ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_reshaped.mnc -keep_real_range -clobber
		noel_do_cmd_new ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_reshaped_temp.mnc mincreshape -dimorder zspace,yspace,xspace ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_reshaped.mnc ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_reshaped_temp.mnc			
		mv -f ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_reshaped_temp.mnc ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_reshaped.mnc
		gzip -vf ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_reshaped.mnc

		
		min=`mincstats -quiet -pctT 0.01  ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_reshaped.mnc.gz -mask ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_brain_mask_native_${MOD_INDEX}.mnc.gz -mask_binvalue 1` 
		max=`mincstats -quiet -pctT 99.99 ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_reshaped.mnc.gz -mask ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_brain_mask_native_${MOD_INDEX}.mnc.gz -mask_binvalue 1` 

		minccalc -clobber -expression "((A[0] < $min) ? 0.0 : (A[0] > $max) ? 1.0 : (A[0] - $min)/($max - $min)) * 100;" ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_reshaped.mnc.gz ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_nuc.mnc
	
		gzip -f ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_nuc.mnc	
done

