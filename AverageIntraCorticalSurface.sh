#!/bin/bash

#################################### 
## GM-CSF && GM-WM surfaces 
OUTDIR=${1}
PREFIX=${2}
CASES=${3}
LEFT_RIGHT=${4}
NUMINSURF=${5}
NUMSUBSURF=${6}
NUMMESH=${7}
SAMPLING_SPACE=${8}
POSTFIX=''
POSTFIX_NUMMESH=''
PRINT=''

#################################### T1 volume & Laplacian maps
if [ ${SAMPLING_SPACE} == 'native' ]; then
	POSTFIX='_native'

fi

MULT_MOD=( t1 flair ir ) #t2cor-0.4
for MOD_INDEX in "${MULT_MOD[@]}"
do
	SURF_SET=( gray white )
	for SURF_INDEX in "${SURF_SET[@]}"
	do
		SURF_DATA=''
		for CASE in `cat ${CASES}`;
		do
			${PRINT} gunzip -vf ${OUTDIR}/${CASE}/surfaces/${PREFIX}_${CASE}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}.obj.gz
			${PRINT} sphere_resample_obj ${OUTDIR}/${CASE}/surfaces/${PREFIX}_${CASE}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}.obj ${OUTDIR}/${CASE}/xfm/${PREFIX}_${CASE}_${LEFT_RIGHT}_surfmap.sm ${OUTDIR}/${CASE}/surfaces/${PREFIX}_${CASE}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_res.obj
			SURF_DATA=${SURF_DATA}' '${OUTDIR}/${CASE}/surfaces/${PREFIX}_${CASE}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_res.obj			
			${PRINT} gzip -vf ${OUTDIR}/${CASE}/surfaces/${PREFIX}_${CASE}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}.obj
		done
		
		${PRINT} average_surfaces ${OUTDIR}/average_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}_${MOD_INDEX}.obj none none 1 ${SURF_DATA}
		${PRINT} rm -rf ${OUTDIR}/*/surfaces/*${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_res.obj
		obj2vtk.sh ${OUTDIR}/average_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}_${MOD_INDEX}.obj ${OUTDIR}/average_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}_${MOD_INDEX}.vtk 
	done

	## Intracortical surfaces
	for (( iter=1; iter <= ${NUMINSURF}; iter++ ))
	do
		SURF_DATA=''
		for CASE in `cat ${CASES}`;
		do
			${PRINT} gunzip -vf ${OUTDIR}/${CASE}/surfaces/${PREFIX}_${CASE}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}.obj.gz
			${PRINT} sphere_resample_obj ${OUTDIR}/${CASE}/surfaces/${PREFIX}_${CASE}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}.obj ${OUTDIR}/${CASE}/xfm/${PREFIX}_${CASE}_${LEFT_RIGHT}_surfmap.sm ${OUTDIR}/${CASE}/surfaces/${PREFIX}_${CASE}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_res.obj 
			SURF_DATA=${SURF_DATA}' '${OUTDIR}/${CASE}/surfaces/${PREFIX}_${CASE}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_res.obj 		
			${PRINT} gzip -vf ${OUTDIR}/${CASE}/surfaces/${PREFIX}_${CASE}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}.obj
		done	
		
		${PRINT} average_surfaces ${OUTDIR}/average_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}_${MOD_INDEX}.obj none none 1 ${SURF_DATA}
		${PRINT} rm -rf ${OUTDIR}/*/surfaces/*${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_res.obj
		obj2vtk.sh ${OUTDIR}/average_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}_${MOD_INDEX}.obj ${OUTDIR}/average_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}_${MOD_INDEX}.vtk
	done

	## immediate subcortical surfaces
	for (( iter=1; iter <= ${NUMSUBSURF}; iter++ ))
	do	
		SURF_DATA=''
		for CASE in `cat ${CASES}`;
		do
			${PRINT} gunzip -vf ${OUTDIR}/${CASE}/surfaces/${PREFIX}_${CASE}_white_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}.obj.gz
			${PRINT} sphere_resample_obj ${OUTDIR}/${CASE}/surfaces/${PREFIX}_${CASE}_white_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}.obj ${OUTDIR}/${CASE}/xfm/${PREFIX}_${CASE}_${LEFT_RIGHT}_surfmap.sm ${OUTDIR}/${CASE}/surfaces/${PREFIX}_${CASE}_white_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_res.obj 
			SURF_DATA=${SURF_DATA}' '${OUTDIR}/${CASE}/surfaces/${PREFIX}_${CASE}_white_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_res.obj
			${PRINT} gzip -vf ${OUTDIR}/${CASE}/surfaces/${PREFIX}_${CASE}_white_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}.obj
		done	
		
		${PRINT} average_surfaces ${OUTDIR}/average_white_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}_${MOD_INDEX}.obj none none 1 ${SURF_DATA}
		${PRINT} rm -rf ${OUTDIR}/*/surfaces/*${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_res.obj		
		obj2vtk.sh ${OUTDIR}/average_white_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}_${MOD_INDEX}.obj ${OUTDIR}/average_white_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}_${MOD_INDEX}.vtk
	done		
done
