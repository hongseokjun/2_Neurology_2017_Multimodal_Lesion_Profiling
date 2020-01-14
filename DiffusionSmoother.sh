#!/bin/bash
print_help()
{
echo "
-----------------------------------------------
`basename $0` <OUT_DIR> <PREFIX> <ID> <LEFT_RIGHT> <NUM_INSURF> <NUM_SUBSURF> <NUMMESH> <SAMPLING_SPACE> <KERNEL> <PARAMETRIC> <USELIGHT>
-----------------------------------------------
This script smooths the sampled data on each surface using surface-based kernel to increase SNR.
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
if [ $# -lt 10 ]
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
KERNEL=${9}
PARAMETRIC=${10}
USELIGHT=${11}
RESULT_DIR=${OUTDIR}/${ID}/measurement/
POSTFIX=''

#################################### T1 volume & Laplacian maps
if [ ${SAMPLING_SPACE} == 'native' ]; then
	POSTFIX='_native'

fi

if [ ${PARAMETRIC} == 'linear' ]; then
	DIFFUSION=0
elif [ ${PARAMETRIC} == 'quadratic' ]; then	
	DIFFUSION=1
elif [ ${PARAMETRIC} == 'hk' ]; then
	DIFFUSION=2
fi

#################################### intra-/sub-cortical surfaces
## GM-CSF && GM-WM surfaces 
SURF_SET=( gray white )
FEATURE_SET=( RI RI_corrected pg pg_GM_WM tg ct FA MD pg_abs_dist )
MULT_MOD=( t1 flair ir dti t1_flair_ratio ) #t2cor-0.4

for MOD_INDEX in "${MULT_MOD[@]}"
do
	MOD_INDEX2=${MOD_INDEX}
	if [ ${MOD_INDEX} == 't1_flair_ratio' ]
	then
		MOD_INDEX2=t1
	fi
	
	for SURF_INDEX in "${SURF_SET[@]}"
	do
		DATASET=''
		for FEATURE_INDEX in "${FEATURE_SET[@]}"
		do
			if [ -e ${RESULT_DIR}/${PREFIX}_${ID}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_${FEATURE_INDEX}.txt ] ; then
				noel_do_cmd_new ${RESULT_DIR}/${PREFIX}_${ID}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt diffuse -kernel ${KERNEL} -iterations 1000 -parametric ${DIFFUSION} ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX2}.obj ${RESULT_DIR}/${PREFIX}_${ID}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_${FEATURE_INDEX}.txt ${RESULT_DIR}/${PREFIX}_${ID}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt
			fi
		done	
	done

	## Intracortical surfaces
	for (( iter=1; iter <= ${NUMINSURF}; iter++ ))
	do
		DATASET=''
		for FEATURE_INDEX in "${FEATURE_SET[@]}"
		do
			if [ -e ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_${FEATURE_INDEX}.txt ] ; then
				noel_do_cmd_new ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt diffuse -kernel ${KERNEL} -iterations 1000 -parametric ${DIFFUSION} ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX2}.obj ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_${FEATURE_INDEX}.txt ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt			
			fi
		done	

		## if there are deformed surfaces, then ...
		if [ ${USELIGHT} -eq 0 ]
		then
			DATASET=''
			for FEATURE_INDEX in "${FEATURE_SET[@]}"
			do			
				if [ -e ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}_deformed${POSTFIX}_${MOD_INDEX}_${FEATURE_INDEX}.txt ] ; then
					noel_do_cmd_new ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}_deformed${POSTFIX}_${MOD_INDEX}_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt diffuse -kernel ${KERNEL} -iterations 1000 -parametric ${DIFFUSION} ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}_deformed${POSTFIX}_${MOD_INDEX2}.obj ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}_deformed${POSTFIX}_${MOD_INDEX}_${FEATURE_INDEX}.txt ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}_deformed${POSTFIX}_${MOD_INDEX}_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt							
				fi
			done	
		fi							  
	done

	## immediate subcortical surfaces
	for (( iter=1; iter <= ${NUMSUBSURF}; iter++ ))
	do	
		DATASET=''
		for FEATURE_INDEX in "${FEATURE_SET[@]}"
		do	
			if [ -e ${RESULT_DIR}/${PREFIX}_${ID}_white_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_${FEATURE_INDEX}.txt ] ; then
				noel_do_cmd_new ${RESULT_DIR}/${PREFIX}_${ID}_white_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt diffuse -kernel ${KERNEL} -iterations 1000 -parametric ${DIFFUSION} ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_white_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX2}.obj ${RESULT_DIR}/${PREFIX}_${ID}_white_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_${FEATURE_INDEX}.txt ${RESULT_DIR}/${PREFIX}_${ID}_white_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_${MOD_INDEX}_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt	
			fi	
		done		  
	done		
done

############################################################################################
FEATURE_SET=( sd sd2 mc )

for SURF_INDEX in "${SURF_SET[@]}"
do
	DATASET=''
	for FEATURE_INDEX in "${FEATURE_SET[@]}"
	do
		if [ -e ${RESULT_DIR}/${PREFIX}_${ID}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}.txt -a ${FEATURE_INDEX} == 'sd' ] ; then
			cp ${RESULT_DIR}/${PREFIX}_${ID}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}.txt ${RESULT_DIR}/${PREFIX}_${ID}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt
		fi

		if [ -e ${RESULT_DIR}/${PREFIX}_${ID}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}.txt -a ${SURF_INDEX} == 'white' -a ${FEATURE_INDEX} == 'sd2' ] ; then
			noel_do_cmd_new ${RESULT_DIR}/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt diffuse -kernel ${KERNEL} -iterations 1000 -parametric ${DIFFUSION} ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1.obj ${RESULT_DIR}/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}.txt ${RESULT_DIR}/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt			
		fi				
		
		if [ -e ${RESULT_DIR}/${PREFIX}_${ID}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}.txt -a ${FEATURE_INDEX} == 'mc' ] ; then
			noel_do_cmd_new ${RESULT_DIR}/${PREFIX}_${ID}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt diffuse -kernel ${KERNEL} -iterations 1000 -parametric ${DIFFUSION} ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1.obj ${RESULT_DIR}/${PREFIX}_${ID}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}.txt ${RESULT_DIR}/${PREFIX}_${ID}_${SURF_INDEX}_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt			
		fi	
	done	
done

## Intracortical surfaces
for (( iter=1; iter <= ${NUMINSURF}; iter++ ))
do
	DATASET=''
	for FEATURE_INDEX in "${FEATURE_SET[@]}"
	do
		if [ -e ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}.txt -a ${FEATURE_INDEX} == 'sd' ] ; then
			cp ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}.txt ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt			
		fi
		
		if [ -e ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}.txt -a ${FEATURE_INDEX} == 'mc' ] ; then
			noel_do_cmd_new ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt diffuse -kernel ${KERNEL} -iterations 1000 -parametric ${DIFFUSION} ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1.obj ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}.txt ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt	
		fi	
	done	
	
	## if there are deformed surfaces, then ...
	if [ ${USELIGHT} -eq 0 ]
	then
		DATASET=''
		for FEATURE_INDEX in "${FEATURE_SET[@]}"
		do			
			if [ -e ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}_deformed${POSTFIX}_t1_${FEATURE_INDEX}.txt ] ; then
				cp ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}_deformed${POSTFIX}_t1_${FEATURE_INDEX}.txt ${RESULT_DIR}/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}_deformed${POSTFIX}_t1_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt							
			fi
		done	
	fi							  
done

## immediate subcortical surfaces
for (( iter=1; iter <= ${NUMSUBSURF}; iter++ ))
do	
	DATASET=''
	for FEATURE_INDEX in "${FEATURE_SET[@]}"
	do	
		if [ -e ${RESULT_DIR}/${PREFIX}_${ID}_white_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}.txt ] ; then
			cp ${RESULT_DIR}/${PREFIX}_${ID}_white_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}.txt ${RESULT_DIR}/${PREFIX}_${ID}_white_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_${FEATURE_INDEX}_${PARAMETRIC}_sm_${KERNEL}.txt	
		fi	
	done		  
done		
