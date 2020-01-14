#!/bin/sh

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

# ID=027
# PREFIX=mcd
# BASEDIR=/data/noel/noel6/CTFCD-1.2.0_64/
# OUTDIR=/local_raid/seokjun/99_temp/DTI_processing/FCD/

ID=${1}
PREFIX=${2}
BASEDIR=${3}
OUTDIR=${4}

if [ ${5} == 1 ]; then
	KEEP=-keep_tmp
else
	KEEP=''
fi

mkdir ${OUTDIR}/${ID}/Noel_DTI_nlin
gunzip -vf ${BASEDIR}/NoelCIVET/${ID}/classify/${PREFIX}_${ID}_pve_wm.mnc.gz
gunzip -vf ${BASEDIR}/NoelCIVET/${ID}/classify/${PREFIX}_${ID}_pve_gm.mnc.gz
gunzip -vf ${BASEDIR}/NoelCIVET/${ID}/temp/${PREFIX}_${ID}_final_classify_insulaFix.mnc.gz

mincresample -like ${OUTDIR}/${ID}/mncdti/${PREFIX}_${ID}_dtifsl_FA.mnc -transformation ${OUTDIR}/${ID}/spip/${PREFIX}_${ID}_tal_to_t1_to_dti_lin.xfm ${OUTDIR}/${ID}/masks/${PREFIX}_${ID}_brain_mask.mnc ${OUTDIR}/${ID}/Noel_DTI_nlin/${PREFIX}_${ID}_brain_mask_dti.mnc -nearest -clobber
minccalc -expr 'if(A[0]+A[1]*0.4<1) out=A[0]+A[1]*0.4 else out=1' ${BASEDIR}/NoelCIVET/${ID}/classify/${PREFIX}_${ID}_pve_wm.mnc ${BASEDIR}/NoelCIVET/${ID}/classify/${PREFIX}_${ID}_pve_gm.mnc ${OUTDIR}/${ID}/Noel_DTI_nlin/${PREFIX}_${ID}_pve_wm_0.4gm.mnc -clobber
minccalc -expr 'if(A[0]>1) out=A[1]' -signed -short ${BASEDIR}/NoelCIVET/${ID}/temp/${PREFIX}_${ID}_final_classify_insulaFix.mnc ${OUTDIR}/${ID}/Noel_DTI_nlin/${PREFIX}_${ID}_pve_wm_0.4gm.mnc ${OUTDIR}/${ID}/Noel_DTI_nlin/${PREFIX}_${ID}_pve_wm_0.4gm_masked.mnc -clobber


gunzip -vf ${OUTDIR}/${ID}/mncdti/${PREFIX}_${ID}_dtifsl_*.mnc.gz
noel_dti_to_tal.sh ${OUTDIR}/${ID}/mncdti/${PREFIX}_${ID}_dtifsl_FA.mnc ${OUTDIR}/${ID}/mncdti/${PREFIX}_${ID}_dtifsl_B0.mnc ${OUTDIR}/${ID}/mncdti/${PREFIX}_${ID}_dtifsl_MD.mnc ${OUTDIR}/${ID}/Noel_DTI_nlin/${PREFIX}_${ID}_pve_wm_0.4gm_masked.mnc ${OUTDIR}/${ID}/Noel_DTI_nlin/${PREFIX}_${ID}_dti2tal.xfm -mask ${OUTDIR}/${ID}/Noel_DTI_nlin/${PREFIX}_${ID}_brain_mask_dti.mnc ${KEEP}
xfminvert ${OUTDIR}/${ID}/Noel_DTI_nlin/${PREFIX}_${ID}_dti2tal_nlin.xfm ${OUTDIR}/${ID}/Noel_DTI_nlin/${PREFIX}_${ID}_tal2dti_nlin.xfm
mincresample -like ${OUTDIR}/${ID}/mncdti/${PREFIX}_${ID}_dtifsl_FA.mnc -transformation ${OUTDIR}/${ID}/Noel_DTI_nlin/${PREFIX}_${ID}_tal2dti_nlin.xfm  ${OUTDIR}/${ID}/spip/${PREFIX}_${ID}_t1_final.mnc ${OUTDIR}/${ID}/Noel_DTI_nlin/${PREFIX}_${ID}_t1_final_to_dti.mnc -clobber 
mincresample -like ${OUTDIR}/${ID}/mncdti/${PREFIX}_${ID}_dtifsl_FA.mnc -transformation ${OUTDIR}/${ID}/Noel_DTI_nlin/${PREFIX}_${ID}_tal2dti_nlin.xfm ${BASEDIR}/Lesion/lesion_vol/${PREFIX}_${ID}_seg.mnc ${OUTDIR}/${ID}/Noel_DTI_nlin/${PREFIX}_${ID}_seg_dti.mnc -nearest -clobber 
