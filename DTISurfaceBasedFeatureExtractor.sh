#!/bin/sh

source ${NOELSOFT_DIR}/BashTools/noel_do_cmd
source ${NOELSOFT_DIR}/BashTools/noel_do_echo
source ${NOELSOFT_DIR}/BashTools/noel_do_cmd_new

################## Example for input arguments ##################
# OUTDIR="/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/"
# DTIDIR="/data/noel/noel6/CTFCD-1.2.0_64/dti/ | /local_raid/seokjun/99_temp/DTI_processing/FCD/"
# PREFIX=mcd
# CASE=002
# NUMINSURF=3
# NUMSUBSURF=7
# NUMMESH=81920           # 81920 | 327680
# SAMPLING_SPACE='native'

OUTDIR=${1}
DTIDIR=${2}
PREFIX=${3}
CASE=${4}
NUMMESH=${5}
NUMINSURF=${6}
NUMSUBSURF=${7}
SAMPLING_SPACE=${8}
     
KERNEL=5                  # FWHM for diffusion smoothing 
PARAMETERIC='quadratic'   # DIffusion smooting mode: linear | quadratic | hk 
ITERATION=10              # Barycentric smoothing iteration
print=""

${print} gunzip -vf ${OUTDIR}/${CASE}/surfaces/*81920_t1.obj.gz

#---------------------------------------------------------------------
noel_do_echo "1) Generate subcortical surface generation"
${print} noel_do_cmd SubcorticalLaminarGenerator.sh ${OUTDIR} ${PREFIX} ${CASE} left  ${NUMMESH} ${NUMSUBSURF} 1
${print} noel_do_cmd SubcorticalLaminarGenerator.sh ${OUTDIR} ${PREFIX} ${CASE} right ${NUMMESH} ${NUMSUBSURF} 1

#---------------------------------------------------------------------
noel_do_echo "2) Project cortical/subcortical surfaces into the native space"
${print} noel_do_cmd TransformTALSurfaceIntoNativeDTI.sh ${PREFIX} ${CASE} ${OUTDIR} ${DTIDIR}/${CASE}/Noel_DTI_nlin/ ${NUMMESH} left  ${NUMINSURF} ${NUMSUBSURF}
${print} noel_do_cmd TransformTALSurfaceIntoNativeDTI.sh ${PREFIX} ${CASE} ${OUTDIR} ${DTIDIR}/${CASE}/Noel_DTI_nlin/ ${NUMMESH} right ${NUMINSURF} ${NUMSUBSURF}

#-----------------------------------------------------------------------
noel_do_echo "3) Sample and map DTI indices onto each surface"
${print} noel_do_cmd DTIIndicesSampler.sh ${OUTDIR} ${DTIDIR} ${PREFIX} ${CASE} left  ${NUMMESH} ${NUMINSURF} ${NUMSUBSURF} ${SAMPLING_SPACE}
${print} noel_do_cmd DTIIndicesSampler.sh ${OUTDIR} ${DTIDIR} ${PREFIX} ${CASE} right ${NUMMESH} ${NUMINSURF} ${NUMSUBSURF} ${SAMPLING_SPACE}

#---------------------------------------------------------------------
noel_do_echo "4) Smooth all the features computed using heat kernel smoothing"
${print} noel_do_cmd DiffusionSmoother.sh ${OUTDIR} ${PREFIX} ${CASE} left  ${NUMINSURF} ${NUMSUBSURF} ${NUMMESH} ${SAMPLING_SPACE} ${KERNEL} ${PARAMETERIC} 1
${print} noel_do_cmd DiffusionSmoother.sh ${OUTDIR} ${PREFIX} ${CASE} right ${NUMINSURF} ${NUMSUBSURF} ${NUMMESH} ${SAMPLING_SPACE} ${KERNEL} ${PARAMETERIC} 1

#---------------------------------------------------------------------
noel_do_echo "5) Resample smoothed features using surface-based registration for correspondence across subjects"
${print} noel_do_cmd SurfaceBasedResampler.sh ${OUTDIR} ${PREFIX} ${CASE} left ${NUMINSURF} ${NUMSUBSURF} ${NUMMESH} ${SAMPLING_SPACE} ${KERNEL} ${PARAMETERIC}
${print} noel_do_cmd SurfaceBasedResampler.sh ${OUTDIR} ${PREFIX} ${CASE} right ${NUMINSURF} ${NUMSUBSURF} ${NUMMESH} ${SAMPLING_SPACE} ${KERNEL} ${PARAMETERIC}

#---------------------------------------------------------------------
#noel_do_echo "6) Make vtk file formats for visualization"
#${print} noel_do_cmd VTKFileGenerator.sh ${OUTDIR} ${PREFIX} ${CASE} ${SAMPLING_SPACE} left ${NUMMESH} ${NUMINSURF} ${NUMSUBSURF} 1 ${KERNEL} ${PARAMETERIC} 'nosurfreg' 'sm'
#${print} noel_do_cmd VTKFileGenerator.sh ${OUTDIR} ${PREFIX} ${CASE} ${SAMPLING_SPACE} right ${NUMMESH} ${NUMINSURF} ${NUMSUBSURF} 1 ${KERNEL} ${PARAMETERIC} 'nosurfreg' 'sm'

${print} gzip -vf ${OUTDIR}/${CASE}/surfaces/*.obj
