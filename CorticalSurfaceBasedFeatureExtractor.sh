#!/bin/sh

# CorticalSurfaceBasedFeatureExtractor.sh /data/noel/noel6/CTControls-1.2.0_64/NoelCIVET/ /data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/ /data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/04_logs/ TLE 1 3 3 81920 native 1
#
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
#		${PREFIX}_${ID}_gray_surface_left_[81920|327680]_[native| ]_[t1|flair|ir].obj			*
#		${PREFIX}_${ID}_white_surface_left_[81920|327680]_[native| ]_[t1|flair|ir].obj		*
#		${PREFIX}_${ID}_gray_surface_right_[81920|327680]_[native| ]_[t1|flair|ir].obj		*
#		${PREFIX}_${ID}_white_surface_right_[81920|327680]_[native| ]_[t1|flair|ir].obj		*
# classify:
#		${PREFIX}_${ID}_ClassFinal.mnc.gz
#		${PREFIX}_${ID}_final_classify.mnc.gz
#		${PREFIX}_${ID}_ClassFinal_native.mnc.gz
#		${PREFIX}_${ID}_final_classify_native.mnc.gz
#		${PREFIX}_${ID}_final_classify_minc1.mnc.gz							*
#		${PREFIX}_${ID}_final_classify_gm.mnc.gz								*
#		${PREFIX}_${ID}_final_classify_wm.mnc.gz								*
#		${PREFIX}_${ID}_final_classify_gm_dil.mnc.gz							*
#		${PREFIX}_${ID}_final_classify_wm_dil.mnc.gz							*
#		${PREFIX}_${ID}_final_classify_dil.mnc.gz								*
#		${PREFIX}_${ID}_ClassFinal_csf_masked.mnc.gz							*
#		${PREFIX}_${ID}_final_classify_dil_remap.mnc.gz						*
# mask:
#  	${PREFIX}_${ID}_brain_mask.mnc.gz
#  	${PREFIX}_${ID}_brain_mask_native.mnc.gz
#		${PREFIX}_${ID}_brain_mask_native_[flair|ir].mnc.gz	
# temp:
#		${PREFIX}_${ID}_[flair|ir]_reshaped_nuc.imp     		
#		${PREFIX}_${ID}_[flair|ir]_reshaped_nuc.mnc.gz
# 		${PREFIX}_${ID}_laplace_Grad[X|Y|Z].mnc.gz 							*
# 		${PREFIX}_${ID}_laplace_laplace.mnc.gz									*
# 		${PREFIX}_${ID}_laplace_reassigned.mnc.gz								*
# 		${PREFIX}_${ID}_laplace_thickness.mnc.gz								*
# 		${PREFIX}_${ID}_laplace_yezzi.mnc.gz									*
# 		${PREFIX}_${ID}_[left|right]_[81920|326780]_stream_line.obj		*
# 		${PREFIX}_${ID}_[left|right]_[81920|326780]_stream_line.vtk		*
# log:
#		${PREFIX}_${ID}_[left|right]_innerwhite_remapclassify.log
#		${PREFIX}_${ID}_[left|right]_SubcorticalLaminarGeneration.log
#		${PREFIX}_${ID}_[left|right]_[81920|326780]_IntraCorticalGeneration.log
# measurement:
#		${PREFIX}_${ID}_left_right_WM_intervertex_dist_avg_[81920|327680]_t1.txt
#		${PREFIX}_${ID}_left_right_GM_intervertex_dist_avg_[81920|327680]_t1.txt

print_help()
{
echo "
-----------------------------------------------
`basename $0` <BASEDIR> <OUTDIR> <LOGDIR> <PREFIX> <CASE> <NUMINSURF> <NUMSUBSURF> <SAMPLING_SPACE> <USELIGHT>
-----------------------------------------------

This script generates intra/sub-cortical surfaces.
1) For intracortical surfaces, this script equi-distantly positions
   vertices of each surface between GM-CSF and GM-WM boundary while constraining 
   intervertex uniform spacing and deformation magnitude.
2) For subcortical surfaces, this script firstly generates a laplacian field 
   between GM-WM and ventricle boundary and deforms inward the GM-WM surface
	to get subcortical surfaces according to the Laplacian field. 
	Each iteration brings the surface 1mm below.
	  
example:
$0 '/data/noel/noel6/CTControls-1.2.0_64/', ...
   'TLE', ...
	'0365_1', ...
	'left', ...
	3, ...
	3, ...
	'/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/result/', ...
	81920, ...
	0.2, ...
	1.0, ...
	0.9, ...
	'native', ...
   1;

Tested for 64bit platforms!

hong.seok.jun@gmail.com
November 2012
March    2012
"
}

source ${NOELSOFT_DIR}/BashTools/noel_do_cmd
source ${NOELSOFT_DIR}/BashTools/noel_do_echo
source ${NOELSOFT_DIR}/BashTools/noel_do_cmd_new

################## Example for input arguments ##################
# BASEDIR="/data/noel/noel6/CTControls-1.2.0_64/NoelCIVET/"
# OUTDIR="/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/"
# LOGDIR="/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/04_logs/"
# PREFIX=TLE
# CASE=301_1
# NUMINSURF=3
# NUMSUBSURF=3
# NUMMESH=81920           # 81920 | 327680
# SAMPLING_SPACE="native" # native | tal
# UseLightVersion=1       # 1: no deformation of surfaces | 0: deform surfaces

BASEDIR=${1}
OUTDIR=${2}
LOGDIR=${3}
PREFIX=${4}
CASE=${5}
NUMMESH=${6}
NUMINSURF=${7}
NUMSUBSURF=${8}
SAMPLING_SPACE=${9}     
UseLightVersion=${10}  
     
KERNEL=5                  # FWHM for diffusion smoothing 
PARAMETERIC='quadratic'   # DIffusion smooting mode: linear | quadratic | hk 
ITERATION=10              # Barycentric smoothing iteration
print=""

# DO NOT EDIT parameters below unless you figure out the process of surface generation.
THRES_ROUGH=0.2
WJt=1.0
WJd=0.9

gunzip -vf ${OUTDIR}/${CASE}/*/*.obj.gz

if [ ${NUMMESH} == 81920 ];
then
	# The number of Mesh: 81920 | ${THRES_ROUGH} ${WJt} ${WJd} are unused...
	IntraCorticalSurfaceGenerator.sh ${OUTDIR} ${PREFIX} ${CASE} left ${NUMINSURF} ${NUMSUBSURF} 81920 ${THRES_ROUGH} ${WJt} ${WJd} ${SAMPLING_SPACE} ${UseLightVersion}
	IntraCorticalSurfaceGenerator.sh ${OUTDIR} ${PREFIX} ${CASE} right ${NUMINSURF} ${NUMSUBSURF} 81920 ${THRES_ROUGH} ${WJt} ${WJd} ${SAMPLING_SPACE} ${UseLightVersion}
	FeatureExtractor.sh ${OUTDIR} ${PREFIX} ${CASE} left ${NUMINSURF} ${NUMSUBSURF} 81920 ${SAMPLING_SPACE} ${KERNEL} ${PARAMETERIC} ${ITERATION} ${UseLightVersion}
	FeatureExtractor.sh ${OUTDIR} ${PREFIX} ${CASE} right ${NUMINSURF} ${NUMSUBSURF} 81920 ${SAMPLING_SPACE} ${KERNEL} ${PARAMETERIC} ${ITERATION} ${UseLightVersion}	
elif [ ${NUMMESH} == 327680 ];
then
	# The number of Mesh: 327680
	IntraCorticalSurfaceGenerator.sh ${OUTDIR} ${PREFIX} ${CASE} left ${NUMINSURF} ${NUMSUBSURF} 327680 ${THRES_ROUGH} ${WJt} ${WJd} ${SAMPLING_SPACE} ${UseLightVersion}
	IntraCorticalSurfaceGenerator.sh ${OUTDIR} ${PREFIX} ${CASE} right ${NUMINSURF} ${NUMSUBSURF} 327680 ${THRES_ROUGH} ${WJt} ${WJd} ${SAMPLING_SPACE} ${UseLightVersion}
	FeatureExtractor.sh ${OUTDIR} ${PREFIX} ${CASE} left ${NUMINSURF} ${NUMSUBSURF} 327680 ${SAMPLING_SPACE} ${KERNEL} ${PARAMETERIC} ${ITERATION} ${UseLightVersion}
	FeatureExtractor.sh ${OUTDIR} ${PREFIX} ${CASE} right ${NUMINSURF} ${NUMSUBSURF} 327680 ${SAMPLING_SPACE} ${KERNEL} ${PARAMETERIC} ${ITERATION} ${UseLightVersion}
fi

noel_do_cmd VTKFileGenerator.sh ${OUTDIR} ${PREFIX} ${CASE} ${SAMPLING_SPACE} left ${NUMMESH} ${NUMINSURF} ${NUMSUBSURF} ${UseLightVersion} ${KERNEL} ${PARAMETERIC} 'nosurfreg' 'sm'
noel_do_cmd VTKFileGenerator.sh ${OUTDIR} ${PREFIX} ${CASE} ${SAMPLING_SPACE} right ${NUMMESH} ${NUMINSURF} ${NUMSUBSURF} ${UseLightVersion} ${KERNEL} ${PARAMETERIC} 'nosurfreg' 'sm'

gzip -vf ${OUTDIR}/${CASE}/*/*.mnc
#gzip -vf ${OUTDIR}/${CASE}/*/*.obj
