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

#---------------------------------------------------------------------
noel_do_echo "1) Make directories"
if [ ! -d ${OUTDIR}/${ID}/logs ];
then	
	mkdir ${OUTDIR}/${ID}
	mkdir ${OUTDIR}/${ID}/native
	mkdir ${OUTDIR}/${ID}/final
	mkdir ${OUTDIR}/${ID}/surfaces
	mkdir ${OUTDIR}/${ID}/xfm
	mkdir ${OUTDIR}/${ID}/measurement
	mkdir ${OUTDIR}/${ID}/classify
	mkdir ${OUTDIR}/${ID}/logs
	mkdir ${OUTDIR}/${ID}/temp
	mkdir ${OUTDIR}/${ID}/mask
fi

noel_do_echo "2) Copy native files and make the vtk format of file for visualization later"
MULT_MOD=( flair ) # t1-0.6 dti rs

for MOD_INDEX in "${MULT_MOD[@]}"
do
	#BASEDIR_temp=/data/noel/noel6/native3T/native_${NATDIR}
   BASEDIR_temp=/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/02_org_data/
	
	noel_do_echo "   - ${MOD_INDEX} is reshaped ..."
	cp ${BASEDIR_temp}/${PREFIX}_${ID}_${MOD_INDEX}.mnc* ${OUTDIR}/${ID}/native
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

	
	noel_do_echo "   - NUC and Clamp, and Coregister ${MOD_INDEX} to T1_tal ..."		
	#### 1. first nuc and clamp without mask and coregistration
	noel_do_cmd_new ${OUTDIR}/${ID}/temp/${PREFIX}_${ID}_${MOD_INDEX}_reshaped_clamp.mnc.gz spip3T-nuc_and_clamp.sh ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_reshaped.mnc.gz -o ${OUTDIR}/${ID}/temp/
	#noel_do_cmd_new ${OUTDIR}/${ID}/xfm/${PREFIX}_${ID}_${MOD_INDEX}_to_tal_concat.xfm spip3T-coreg_amri_new.sh ${OUTDIR}/${ID}/temp/${PREFIX}_${ID}_${MOD_INDEX}_reshaped_clamp.mnc.gz -nuc ${OUTDIR}/${ID}/native/ -tal ${OUTDIR}/${ID}/final/ -xfm ${OUTDIR}/${ID}/xfm/ -reg ${OUTDIR}/${ID}/final/
	
	noel_do_cmd_new ${OUTDIR}/${ID}/xfm/${PREFIX}_${ID}_${MOD_INDEX}_to_tal_concat.xfm mincresample -clobber ${OUTDIR}/${ID}/temp/${PREFIX}_${ID}_${MOD_INDEX}_reshaped_clamp.mnc.gz /tmp/tmp_${MOD_INDEX}_$$.mnc -like /data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/${ID}/native/${PREFIX}_${ID}_t1_nuc.mnc.gz
	noel_do_cmd_new ${OUTDIR}/${ID}/xfm/${PREFIX}_${ID}_${MOD_INDEX}_to_t1.xfm mritoself -clobber -lsq6 -mi -far -model mni_icbm152_t1_tal_nlin_sym_09a -modeldir /data/noel/noel6/Mar-07-2007-x64//install/data/mni_autoreg/ -target_talxfm /data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/${ID}/xfm/${PREFIX}_${ID}_t1_to_tal.xfm /tmp/tmp_${MOD_INDEX}_$$.mnc /data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/${ID}/native/${PREFIX}_${ID}_t1_nuc.mnc.gz ${OUTDIR}/${ID}/xfm/${PREFIX}_${ID}_${MOD_INDEX}_to_t1.xfm
	noel_do_cmd_new ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_to_t1.mnc mincresample -clobber -like /data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/${ID}/native/${PREFIX}_${ID}_t1_nuc.mnc.gz -transform ${OUTDIR}/${ID}/xfm/${PREFIX}_${ID}_${MOD_INDEX}_to_t1.xfm  ${OUTDIR}/${ID}/temp/${PREFIX}_${ID}_${MOD_INDEX}_reshaped_clamp.mnc.gz ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_to_t1.mnc
	gzip -vf ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_to_t1.mnc
	
	noel_do_cmd_new ${OUTDIR}/${ID}/xfm/${PREFIX}_${ID}_${MOD_INDEX}_to_tal_concat.xfm xfmconcat ${OUTDIR}/${ID}/xfm/${PREFIX}_${ID}_${MOD_INDEX}_to_t1.xfm /data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/${ID}/xfm/${PREFIX}_${ID}_t1_to_tal.xfm ${OUTDIR}/${ID}/xfm/${PREFIX}_${ID}_${MOD_INDEX}_to_tal_concat.xfm
	noel_do_cmd_new ${OUTDIR}/${ID}/final/${PREFIX}_${ID}_${MOD_INDEX}_to_t1_tal.mnc mincresample -clobber -like /data/noel/noel6/Mar-07-2007-x64/install/data/mni_autoreg/mni_icbm152_t1_tal_nlin_sym_09a.mnc -transform ${OUTDIR}/${ID}/xfm/${PREFIX}_${ID}_${MOD_INDEX}_to_tal_concat.xfm  ${OUTDIR}/${ID}/temp/${PREFIX}_${ID}_${MOD_INDEX}_reshaped_clamp.mnc.gz ${OUTDIR}/${ID}/final/${PREFIX}_${ID}_${MOD_INDEX}_to_t1_tal.mnc
	gzip -vf ${OUTDIR}/${ID}/final/${PREFIX}_${ID}_${MOD_INDEX}_to_t1_tal.mnc

	noel_do_cmd_new ${OUTDIR}/${ID}/xfm/${PREFIX}_${ID}_tal_to_${MOD_INDEX}.xfm xfminvert ${OUTDIR}/${ID}/xfm/${PREFIX}_${ID}_${MOD_INDEX}_to_tal_concat.xfm ${OUTDIR}/${ID}/xfm/${PREFIX}_${ID}_tal_to_${MOD_INDEX}.xfm

	#### 2. transfrom T1 final mask to the native space of multimodal image
	noel_do_cmd_new ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_brain_mask_native_${MOD_INDEX}.mnc.gz mincresample -like ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_reshaped.mnc.gz -nearest_neighbour -transformation ${OUTDIR}/${ID}/xfm/${PREFIX}_${ID}_tal_to_${MOD_INDEX}.xfm /data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/${ID}/mask/${PREFIX}_${ID}_brain_mask.mnc.gz ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_brain_mask_native_${MOD_INDEX}.mnc
	gzip -vf ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_brain_mask_native_${MOD_INDEX}.mnc
   noel_do_cmd_new ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_brain_mask_native_${MOD_INDEX}_dilated.mnc.gz mincmorph -clobber -successive DDD ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_brain_mask_native_${MOD_INDEX}.mnc.gz ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_brain_mask_native_${MOD_INDEX}_dilated.mnc
	gzip -vf ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_brain_mask_native_${MOD_INDEX}_dilated.mnc
	rm -rf ${OUTDIR}/${ID}/temp/${PREFIX}_${ID}_${MOD_INDEX}_reshaped_clamp.mnc.gz

	#### 3. second nuc and clamp using mask
	noel_do_cmd_new ${OUTDIR}/${ID}/temp/${PREFIX}_${ID}_${MOD_INDEX}_reshaped_clamp.mnc.gz spip3T-nuc_and_clamp.sh ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_reshaped.mnc.gz -o ${OUTDIR}/${ID}/temp/ -mask ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_brain_mask_native_${MOD_INDEX}_dilated.mnc.gz
	noel_do_cmd_new ${OUTDIR}/${ID}/final/${PREFIX}_${ID}_${MOD_INDEX}_final.mnc.gz mincresample -like ${OUTDIR}/${ID}/final/${PREFIX}_${ID}_${MOD_INDEX}_to_t1_tal.mnc -transformation ${OUTDIR}/${ID}/xfm/${PREFIX}_${ID}_${MOD_INDEX}_to_tal_concat.xfm ${OUTDIR}/${ID}/temp/${PREFIX}_${ID}_${MOD_INDEX}_reshaped_clamp.mnc.gz ${OUTDIR}/${ID}/final/${PREFIX}_${ID}_${MOD_INDEX}_final.mnc
	gzip -vf ${OUTDIR}/${ID}/final/${PREFIX}_${ID}_${MOD_INDEX}_final.mnc

	#### 4. clear files
#	noel_do_cmd_new ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_nuc.mnc.gz mv ${OUTDIR}/${ID}/temp/${PREFIX}_${ID}_${MOD_INDEX}_reshaped_clamp.mnc.gz ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_nuc.mnc.gz
#	rm -rf ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}.mnc.gz
#	rm -rf ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_reshaped.mnc.gz
#	rm -rf ${OUTDIR}/${ID}/native/${PREFIX}_${ID}_${MOD_INDEX}_to_t1.mnc.gz		
#	rm -rf ${OUTDIR}/${ID}/final/${PREFIX}_${ID}_${MOD_INDEX}_stx_norm.mnc.gz
#	rm -rf ${OUTDIR}/${ID}/mask/${PREFIX}_${ID}_brain_mask_native_${MOD_INDEX}_dilated.mnc.gz
done

