#!/bin/bash
print_help()
{
echo "
-----------------------------------------------
`basename $0` <OUTDIR> <PREFIX> <ID> <LEFT_RIGHT> <NUMMESH> <NUMINSURF> <SAMPLING_SPACE> <ITERATION>
-----------------------------------------------
This script computes curvature for each cortical surface
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

if [ ${SAMPLING_SPACE} == 'native' ]; then
	POSTFIX='_native'
elif [ ${SAMPLING_SPACE} == 'tal' ]; then
	POSTFIX=''
fi

#---------------------------------------------------------------------
noel_do_echo "Calculation of the curvature in the defined space"
noel_do_echo "1) Perform barycentric smoothing"
noel_do_cmd_new ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj smooth_surface_barycentric_new.sh ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1.obj ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj ${ITERATION} 

for (( iter=1; iter <= ${NUMINSURF}; iter++ ))
do
	noel_do_cmd_new ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj smooth_surface_barycentric_new.sh ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1.obj ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj ${ITERATION} 
done

noel_do_cmd_new ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj smooth_surface_barycentric_new.sh ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1.obj ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj ${ITERATION}

#---------------------------------------------------------------------
noel_do_echo "2) Compute mean curvature"
noel_do_cmd_new ${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_mc.txt /data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/mean_curvature ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj ${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_mc.txt

for (( iter=1; iter <= ${NUMINSURF}; iter++ ))
do
	noel_do_cmd_new ${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_mc.txt /data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/mean_curvature ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj ${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_intracortical_surface_${iter}_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_mc.txt
done

noel_do_cmd_new ${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_mc.txt /data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/mean_curvature ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_barycentric_sm_${ITERATION}.obj ${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_mc.txt

matlab -nodisplay <<EOF > ${OUTDIR}/${ID}/logs/${PREFIX}_${ID}_${LEFT_RIGHT}_${NUMMESH}_absolute_curvature.log
		addpath( '${EMMADIR}' );
		addpath( '${SURFSTATDIR}' ); 
		intracortical_surface_number = 3;
		
		for i = 1 : intracortical_surface_number + 2
        if(i == 1)
            basename  = [ '${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_' num2str(${NUMMESH}) '${POSTFIX}_t1_mc.txt' ];
        elseif(i == 5)
            basename  = [ '${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_' num2str(${NUMMESH}) '${POSTFIX}_t1_mc.txt' ];
        else
            basename  = [ '${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_intracortical_surface_' num2str(i-1) '_${LEFT_RIGHT}_' num2str(${NUMMESH}) '${POSTFIX}_t1_mc.txt' ];
        end
		  
		  data = SurfStatReadData(basename);
		  data = abs(data);
		  SurfStatWriteData(basename, data);
    end
EOF
