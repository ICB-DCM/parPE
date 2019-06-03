{% extends "loadleveler_multistep.tpl" %}
{% block before_any_step %}
{{ super() }}
set -x # show commands in log files
set -e # exit on error

# Get job file directory
# NOTE: $BASH_SOURCE does not work here, because script seems to be copied by LoadLeveler
# NOTE: $LOADL_STEP_INITDIR is the submission dir, not necessarily the script dir
SCRIPT_PATH="`dirname \"${LOADL_STEP_COMMAND}\"`"
SCRIPT_PATH="`( cd \"${SCRIPT_PATH}\" && pwd )`"
cd $SCRIPT_PATH
{% endblock %}

{% block after_any_step %}
{{ super() }}
# send custom notification email
#if [ -n "${EMAIL_RECIPIENTS}" ]; then
#   echo -e "`pwd`\n\n `ls -l` \n\n `head -n 200 ${LOADL_STEP_ERR}`" | mail -s "${LOADL_JOB_NAME} finshed" ${EMAIL_RECIPIENTS}
#fi
{% endblock %}
{#


# Input file for optimization
INPUT_DATA_FILE="./datafiles/CV1_hierarchical_withProtein_withSigmaPerProtein_proteinOffset_MS1.h5"
# Prefix for optimization result HDF5 files
OPTIMIZATION_RESULT_PREFIX="results-${LOADL_JOB_NAME}-${LOADL_STEP_ID}/${LOADL_JOB_NAME}-${LOADL_STEP_ID}"
# replace step id by step 0
OPTIMIZATION_RESULT_PREFIX=$(echo "$OPTIMIZATION_RESULT_PREFIX" | sed -re 's/\.[0-9]+$/\.0/')
# File name of result file from master
OPTIMIZATION_RESULT_FILENAME="${OPTIMIZATION_RESULT_PREFIX}_rank00000.h5"
# Blank-separated list of recipients for custom notification mail
# (LoadLeveler only supports a single email address)
EMAIL_RECIPIENTS="daniel.weindl@helmholtz-muenchen.de leonard.schmiester@helmholtz-muenchen.de"
# Time stamp command run between individual steps
TIMESTAMP="echo `date` $LOADL_STEP_NAME >> time-${LOADL_JOB_NAME}.out"

$TIMESTAMP

pwd


echo "### step: ${LOADL_STEP_NAME} ###"

case ${LOADL_STEP_NAME} in
    optimize )
        printenv

        # run optimization
        #MP_INFOLEVEL=6
        poe  ./estimate_ERBB_RAS_AKT_Drugs_r389936_syms_withSigmaPerProtein -o "${OPTIMIZATION_RESULT_PREFIX}" "${INPUT_DATA_FILE}" || dmesg | tail -n 300
        #gdb ./estimate_ERBB_RAS_AKT_Drugs_r389936_syms_withSigmaPerProtein .
        # copy input data to rank 0 result file
        # Need hdf5 version >= 1.10
        # module load hdf5/serial/1.10.0
    ;;

    simulate )
        # OpenMP settings, use 1 thread per core
        export OMP_NUM_THREADS=28
        export KMP_AFFINITY="granularity=core,compact,1"

        SIMULATION_INFILE="${OPTIMIZATION_RESULT_FILENAME}"
        SIMULATION_OUTFILE=simulationsAtOptima.h5
        # TODO read starts from file
        NUMSTARTS=10
        for MS in `seq 0 $((NUMSTARTS - 1))`
        do
            ./simulateERBB_RAS_AKT_Drugs_r389936_withSigma ${SIMULATION_INFILE} ${SIMULATION_OUTFILE} $MS || dmesg | tail -n 300
        done
    ;;
esac



$TIMESTAMP
echo `date` $LOADL_STEP_NAME >> time-${LOADL_JOB_NAME}.out
    #}
