#!/bin/bash

#################################################################
# Files
# native:
#		${PREFIX}_${ID}_t1_nuc.mnc.gz
#		${PREFIX}_${ID}_[flair|ir]_nuc.mnc.gz
#		${PREFIX}_${ID}_[t1|flair|ir]_[GM|WM]_native.mnc.gz							*
#		${PREFIX}_${ID}_ri_[t1|flair|ir]_[gm|wm]_native.mnc.gz						*
# final:
#		${PREFIX}_${ID}_t1_final.mnc.gz	
#		${PREFIX}_${ID}_[flair|ir]_final.mnc.gz
#		${PREFIX}_${ID}_[t1|flair|ir]_[GM}WM]_final.mnc.gz								*
#		${PREFIX}_${ID}_ri_[t1|flair|ir]_[gm|wm].mnc.gz									*
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
#		${PREFIX}_${ID}_gray_surface_left_[81920|327680]_[native| ]_[t1|flair|ir].obj			
#		${PREFIX}_${ID}_white_surface_left_[81920|327680]_[native| ]_[t1|flair|ir].obj		
#		${PREFIX}_${ID}_gray_surface_right_[81920|327680]_[native| ]_[t1|flair|ir].obj		
#		${PREFIX}_${ID}_white_surface_right_[81920|327680]_[native| ]_[t1|flair|ir].obj		
# classify:
#		${PREFIX}_${ID}_ClassFinal.mnc.gz
#		${PREFIX}_${ID}_final_classify.mnc.gz
#		${PREFIX}_${ID}_ClassFinal_native.mnc.gz
#		${PREFIX}_${ID}_final_classify_native.mnc.gz
#		${PREFIX}_${ID}_final_classify_minc1.mnc.gz							
#		${PREFIX}_${ID}_final_classify_gm.mnc.gz								
#		${PREFIX}_${ID}_final_classify_wm.mnc.gz							
#		${PREFIX}_${ID}_final_classify_gm_dil.mnc.gz							
#		${PREFIX}_${ID}_final_classify_wm_dil.mnc.gz							
#		${PREFIX}_${ID}_final_classify_dil.mnc.gz								
#		${PREFIX}_${ID}_ClassFinal_csf_masked.mnc.gz							
#		${PREFIX}_${ID}_final_classify_dil_remap.mnc.gz						
#		${PREFIX}_${ID}_[GM|WM]_final_classify.mnc.gz									*
#		${PREFIX}_${ID}_[t1|flair|ir]_[GM}WM]_native_classify.mnc.gz				*
# mask:
#  	${PREFIX}_${ID}_brain_mask.mnc.gz
#  	${PREFIX}_${ID}_brain_mask_native.mnc.gz
#		${PREFIX}_${ID}_brain_mask_native_[flair|ir].mnc.gz	
# temp:
#		${PREFIX}_${ID}_[flair|ir]_reshaped_nuc.imp     		
#		${PREFIX}_${ID}_[flair|ir]_reshaped_nuc.mnc.gz
# 		${PREFIX}_${ID}_laplace_Grad[X|Y|Z].mnc.gz 							
# 		${PREFIX}_${ID}_laplace_laplace.mnc.gz								
# 		${PREFIX}_${ID}_laplace_reassigned.mnc.gz								
# 		${PREFIX}_${ID}_laplace_thickness.mnc.gz								
# 		${PREFIX}_${ID}_laplace_yezzi.mnc.gz									
# 		${PREFIX}_${ID}_[left|right]_[81920|326780]_stream_line.obj		
# 		${PREFIX}_${ID}_[left|right]_[81920|326780]_stream_line.vtk		
# log:
#		${PREFIX}_${ID}_[left|right]_innerwhite_remapclassify.log
#		${PREFIX}_${ID}_[left|right]_SubcorticalLaminarGeneration.log
#		${PREFIX}_${ID}_[left|right]_[81920|326780]_IntraCorticalGeneration.log
# measurement:
#		${PREFIX}_${ID}_left_right_WM_intervertex_dist_avg_[81920|327680]_t1.txt
#		${PREFIX}_${ID}_left_right_GM_intervertex_dist_avg_[81920|327680]_t1.txt
#		${PREFIX}_${ID}_[gray|intracortical|white]_surface_[|1|2|3]_[left|right]_[81920|327680][_native| ]_[t1|flair|ir]_RI.txt		*

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
noel_do_cmd RelativeIntensityCalculator.sh -prefix ${PREFIX} -ID ${ID} -dir ${OUTDIR} -space ${SAMPLING_SPACE}

#-----------------------------------------------------------------------
noel_do_echo "2) Sample and map RI values onto each surface"
noel_do_cmd RelativeIntensitySampler.sh ${OUTDIR} ${PREFIX} ${ID} ${LEFT_RIGHT} ${NUMMESH} ${NUMINSURF} ${NUMSUBSURF} ${SAMPLING_SPACE}

#-----------------------------------------------------------------------
noel_do_echo "3) Regress out the pve of csf from AI and recalculate RI"
noel_do_echo "   over the only vertices which pve is over 0.5"
noel_do_cmd PartialVolumeEffectRemover.sh ${OUTDIR} ${PREFIX} ${ID} ${SAMPLING_SPACE} ${LEFT_RIGHT} ${NUMMESH} ${NUMINSURF}
