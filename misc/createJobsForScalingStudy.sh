#!/bin/bash
# This script creates job files for scaling studies
# using a template and creating new files for the specified
# number of nodes
#
# Daniel Weindl 2017

OUTPREFIX="scaling"
TPL_CORES=TPL_CORES
TPL_NODES=TPL_NODES

# invalid options:
if [ $# -lt 3 ] || [ $# -gt 4 ] ; then
    echo "USAGE: createJobsForScalingStudy.sh TEMPLATE_FILE CORES_PER_NODE NUM_NODES1,NUM_NODES2,... OUTFILE_PREFIX"
    echo 
    echo "Will replace ${TPL_CORES} and ${TPL_NODES} in the respective copy of TEMPLATE_FILE accordingly."
    exit 1;
fi

# parse arguments
if [ $# == 4 ]; then
    OUTPREFIX=$4
fi

INFILE=$1
CORES_PER_NODE=$2
NODELIST=$3

# generate files
IFS=','
for NODES in $NODELIST
do
    CORES=$[$NODES*$CORES_PER_NODE]
    echo "$NODES nodes -> $CORES cores"
    OUTFILE="${OUTPREFIX}_${NODES}n_${CORES}c.job"
    if [ -f $OUTFILE ]
    then
        echo "  $OUTFILE exists... skipping"
    else
        sed -r s/${TPL_CORES}/${CORES}/ ${INFILE} | \
        sed -r s/${TPL_NODES}/${NODES}/ > ${OUTFILE}
    fi
done