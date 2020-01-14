#!/bin/bash
print_help()
{
echo "
-----------------------------------------------
`basename $0`  <OUTDIR> <PREFIX> <ID> <LEFT_RIGHT> <NUMMESH> <SAMPLING_SPACE>
-----------------------------------------------
This script computes cortical thickness in [native|tal] without smoothing process.
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
if [ $# -lt 6 ]
then
    print_help
    exit 1
fi

# ------------------------------------------------------------------------------
do_cmd() 
{
    local l_command=""
    local l_sep=""
    local l_index=1
    while [ ${l_index} -le $# ]; do
	eval arg=\${$l_index}
	l_command="${l_command}${l_sep}${arg}"
	l_sep=" "
	l_index=$[${l_index}+1]
    done
    echo
    echo "-------------------------------------------------------------------------------"
    echo "${log_header} ${l_command}"
    echo "-------------------------------------------------------------------------------"
    $l_command
}

OUTDIR=${1}
PREFIX=${2}
ID=${3}
LEFT_RIGHT=${4}
NUMMESH=${5}
SAMPLING_SPACE=${6}

#---------------------------------------------------------------------
noel_do_echo "Calculation of the cortical thickness in the defined space"
if [ ${SAMPLING_SPACE} == 'native' ]; then
	POSTFIX='_native'
elif [ ${SAMPLING_SPACE} == 'tal' ]; then
	POSTFIX=''
fi

noel_do_cmd_new ${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_intracortical_surface_2_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_ct.txt cortical_thickness -tlink ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_white_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1.obj ${OUTDIR}/${ID}/surfaces/${PREFIX}_${ID}_gray_surface_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1.obj ${OUTDIR}/${ID}/measurement/${PREFIX}_${ID}_intracortical_surface_2_${LEFT_RIGHT}_${NUMMESH}${POSTFIX}_t1_ct.txt
