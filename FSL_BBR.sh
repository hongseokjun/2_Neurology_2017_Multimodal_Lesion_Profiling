#!/bin/sh
source ${NOELSOFT_DIR}/BashTools/noel_do_cmd
source ${NOELSOFT_DIR}/BashTools/noel_do_echo
source ${NOELSOFT_DIR}/BashTools/noel_do_cmd_new

PREFIX=${1}
CASE=${2}
BASEDIR=${3}

mkdir ${BASEDIR}/${CASE}/bbr/
noel_do_cmd_new ${BASEDIR}/${CASE}/masks/${PREFIX}_${CASE}_brain_mask_native.mnc mincresample -transformation ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_tal_to_t1.xfm -like ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_nuc.mnc -nearest ${BASEDIR}/${CASE}/masks/${PREFIX}_${CASE}_brain_mask.mnc ${BASEDIR}/${CASE}/masks/${PREFIX}_${CASE}_brain_mask_native.mnc

AimsFileConvert -i ${BASEDIR}/${CASE}/masks/${PREFIX}_${CASE}_brain_mask_native.mnc -o ${BASEDIR}/${CASE}/masks/${PREFIX}_${CASE}_brain_mask_native_nii.nii --orient "-1 -1 -1"
AimsFileConvert -i ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_nuc.mnc -o ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_nuc_nii.nii --orient "-1 -1 -1"
AimsFileConvert -i ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_nuc_masked.mnc -o ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_nuc_masked_nii.nii --orient "-1 -1 -1"
AimsFileConvert -i ${BASEDIR}/${CASE}/masks/${PREFIX}_${CASE}_brain_mask.mnc -o ${BASEDIR}/${CASE}/masks/${PREFIX}_${CASE}_brain_mask_nii.nii
AimsFileConvert -i ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_final.mnc -o ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_final_nii.nii

noel_do_cmd_new ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_to_tal_fsl_mat.mat fsl5.0-flirt -inweight ${BASEDIR}/${CASE}/masks/${PREFIX}_${CASE}_brain_mask_native_nii.nii -refweight ${BASEDIR}/${CASE}/masks/${PREFIX}_${CASE}_brain_mask_nii.nii -in ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_nuc_nii.nii -ref ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_final_nii.nii -out ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_to_tal_fsl.nii -omat ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_to_tal_fsl_mat.mat
noel_do_cmd_new ${BASEDIR}/${CASE}/temp/${PREFIX}_${CASE}_final_class_WM.mnc minccalc -expr 'if(A[0]==3){out=1}' ${BASEDIR}/${CASE}/temp/${PREFIX}_${CASE}_final_classify_insulaFix.mnc ${BASEDIR}/${CASE}/temp/${PREFIX}_${CASE}_final_class_WM.mnc -clobber
noel_do_cmd_new ${BASEDIR}/${CASE}/temp/${PREFIX}_${CASE}_final_class_WM_native.mnc mincresample -like ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_nuc.mnc -nearest -transformation ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_tal_to_t1.xfm ${BASEDIR}/${CASE}/temp/${PREFIX}_${CASE}_final_class_WM.mnc ${BASEDIR}/${CASE}/temp/${PREFIX}_${CASE}_final_class_WM_native.mnc
AimsFileConvert -i ${BASEDIR}/${CASE}/temp/${PREFIX}_${CASE}_final_class_WM_native.mnc -o ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_nuc_masked_nii_wmseg.nii --orient "-1 -1 -1"

noel_do_cmd_new ${BASEDIR}/${CASE}/bbr/${PREFIX}_${CASE}_dti_t1_bbr.mat fsl5.0-epi_reg --epi=${BASEDIR}/${CASE}/fsldti/nodif.nii.gz --t1=${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_nuc_nii.nii --t1brain=${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_nuc_masked_nii.nii --out=${BASEDIR}/${CASE}/bbr/${PREFIX}_${CASE}_dti_t1_bbr --echospacing=0.00085 --pedir=-y -v
noel_do_cmd_new ${BASEDIR}/${CASE}/bbr/${PREFIX}_${CASE}_t1_dti_bbr.mat fsl5.0-convert_xfm -omat ${BASEDIR}/${CASE}/bbr/${PREFIX}_${CASE}_t1_dti_bbr.mat -inverse ${BASEDIR}/${CASE}/bbr/${PREFIX}_${CASE}_dti_t1_bbr.mat
noel_do_cmd_new ${BASEDIR}/${CASE}/bbr/${PREFIX}_${CASE}_tal_t1_fsl_mat.mat fsl5.0-convert_xfm -omat ${BASEDIR}/${CASE}/bbr/${PREFIX}_${CASE}_tal_t1_fsl_mat.mat -inverse ${BASEDIR}/${CASE}/spip/${PREFIX}_${CASE}_t1_to_tal_fsl_mat.mat
noel_do_cmd_new ${BASEDIR}/${CASE}/bbr/${PREFIX}_${CASE}_tal_t1_dti_bbr.mat fsl5.0-convert_xfm -omat ${BASEDIR}/${CASE}/bbr/${PREFIX}_${CASE}_tal_t1_dti_bbr.mat -concat ${BASEDIR}/${CASE}/bbr/${PREFIX}_${CASE}_t1_dti_bbr.mat ${BASEDIR}/${CASE}/bbr/${PREFIX}_${CASE}_tal_t1_fsl_mat.mat
