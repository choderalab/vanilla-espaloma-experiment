#!/bin/bash

cyc1=1
cyc2=10
basename="tyk2-esp"
DIR=${PWD}

for (( i=${cyc1}; i<=${cyc2}; i++ ))
do
    j=$(($i - 1))

    # create soft link    
    if [ ${i} == 1 ]
    then
        ln -s ../../prep md0
    fi

    mkdir md${i}
    sed -e 's/@@@JOBNAME@@@/'${basename}''${i}'/g' \
        -e 's/@@@RESTART_PREFIX@@@/..\/md'${j}'/g' \
        lsf-template.sh > ./md${i}/run.sh
    chmod u+x ./md${i}/run.sh

    echo "submit job ${i}"
    cd md${i}

    # make input and submit job
    if [ ${i} == ${cyc1} ]
    then
        bsub < run.sh
        cd ${DIR}
    else
        jobname=${basename}${j}
        bsub -w 'done('${jobname}')' < run.sh
        cd ${DIR}
    fi

    sleep 1
done