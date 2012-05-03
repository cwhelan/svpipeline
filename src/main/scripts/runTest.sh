#!/bin/bash

set -e
set -u

source test_params

hadoop jar ../lib/svpipeline-1.0-SNAPSHOT-exe.jar -Dmapred.reduce.tasks=200 incrementalUpdateSingleEndDeletionScores \
    --inputHDFSDir /user/whelanch/svpipeline/$SAMPLE/$ALIGNMENTS \
    --outputHDFSDir /user/whelanch/svpipeline/$SAMPLE/$NAME \
    --targetIsize $ISIZE --targetIsizeSD $ISIZE_SD --faidx $HDFS_FAI \
    --mapabilityWeighting $MAPABILITY \
    --excludePairsMappingIn $EXCLUSIONS --resolution $RESOLUTION

hadoop jar ../lib/svpipeline-1.0-SNAPSHOT-exe.jar exportWigAndBedFiles \
    --inputHDFSDir /user/whelanch/svpipeline/$SAMPLE/$NAME \
    --faidx $LOCAL_FAI \
    --resolution $RESOLUTION --medianFilterWindow $MEDIAN_FILTER_WINDOW --outputPrefix $NAME

gzip *.wig *.bed

python ../src/main/scripts/evalWigFile.py ${NAME}_piledup_deletion_scores.wig.gz \
    $TRUTH \
    $LOCAL_FAI $MEDIAN_FILTER_WINDOW $THRESHOLD_MIN > ${NAME}_perf.txt