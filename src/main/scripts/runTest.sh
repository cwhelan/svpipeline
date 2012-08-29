#!/bin/bash

set -e
set -u

if ! [ -a ./cloudbreak_params ]
then
    echo "Could not find cloudbreak_params file"
fi

source ./cloudbreak_params

if ! ( [ -n $SAMPLE_NAME ] && [ -n $HDFS_SAMPLE_DIR ] && [ -n $BUILD_DIR ] && [ -n $LOCAL_FAI ] && [ -n $HDFS_FAI ] )
then
	echo cloudbreak_params needs to have:
	echo
	echo SAMPLE_NAME: overall name of the experiment
	echo HDFS_SAMPLE_DIR: HDFS dir to store the score pileup in
	echo BUILD_DIR: location on the local filesystem to build cloudbreak in
	echo LOCAL_FAI: location of the genome .fai on the local filesystem
	echo HDFS_FAI: location of the genome .fai on the HDFS filesystem    
fi

READ_GROUP_FILE=$1
MAPABILITY=$2
FILTER=$3
RESOLUTION=$4
MEDIAN_FILTER_WINDOW=$5
ALIGNER=$6
MAX_INSERT=$7
TRUTH=$8
THRESHOLD_MIN=$9
TAG=$10

pushd $BUILD_DIR
git pull
mvn clean
mvn assembly:assembly
SHORT_GIT_TAG=`git rev-parse --short HEAD`
popd

NAME=`basename $READ_GROUP_FILE`_`basename $MAPABILITY`_`basename $FILTER`_${MAX_INSERT}_${RESOLUTION}_${MEDIAN_FILTER_WINDOW}_${ALIGNER}_${SHORT_GIT_TAG}_${TAG}

echo Experiment name: $NAME

mkdir $NAME
pushd $NAME

hadoop jar $BUILD_DIR/target/cloudbreak-1.0-SNAPSHOT-exe.jar -Dmapred.reduce.tasks=200 incrementalUpdateSingleEndDeletionScores \
    --inputFileDescriptor $READ_GROUP_FILE \
    --outputHDFSDir $HDFS_DIRNAME/$NAME \
    --faidx $HDFS_FAI \
    --maxInsertSize $MAX_INSERT_SIZE \
    --mapabilityWeighting $MAPABILITY \
    --excludePairsMappingIn $EXCLUSIONS \
    --resolution $RESOLUTION \
    --aligner $ALIGNER

hadoop jar ../lib/cloudbreak-1.0-SNAPSHOT-exe.jar exportWigAndBedFiles \
    --inputHDFSDir /user/whelanch/svpipeline/$SAMPLE/$NAME \
    --faidx $LOCAL_FAI \
    --resolution $RESOLUTION --medianFilterWindow $MEDIAN_FILTER_WINDOW --outputPrefix $NAME

gzip *.wig *.bed

python ../src/main/scripts/evalWigFile.py ${NAME}_piledup_deletion_scores.wig.gz \
    $TRUTH \
    $LOCAL_FAI $MEDIAN_FILTER_WINDOW $THRESHOLD_MIN > ${NAME}_perf.txt