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
SHORT_NAME=${10}

pushd $BUILD_DIR
git pull
mvn clean
mvn assembly:assembly
SHORT_GIT_TAG=`git rev-parse --short HEAD`
popd

NAME=${SHORT_NAME}_`basename $READ_GROUP_FILE | awk -F'.' '{print $1}'`_`basename $MAPABILITY | awk -F'.' '{print $1}'`_`basename $FILTER | awk -F'.' '{print $1}'`_${MAX_INSERT}_${RESOLUTION}_${MEDIAN_FILTER_WINDOW}_${ALIGNER}_${SHORT_GIT_TAG}

echo Experiment name: $NAME

mkdir $NAME
pushd $NAME

exec > >(tee -a test.log)

cat <<EOF
Parameters:

EXPERIMENT=$NAME

SOFTWARE_VERSION=$SHORT_GIT_TAG

SAMPLE_NAME=$SAMPLE_NAME
HDFS_SAMPLE_DIR=$HDFS_SAMPLE_DIR
BUILD_DIR=$BUILD_DIR
LOCAL_FAI=$LOCAL_FAI
HDFS_FAI=$HDFS_FAI
READ_GROUP_FILE=$READ_GROUP_FILE
MAPABILITY=$MAPABILITY
FILTER=$FILTER
RESOLUTION=$RESOLUTION
MEDIAN_FILTER_WINDOW=$MEDIAN_FILTER_WINDOW
ALIGNER=$ALIGNER
MAX_INSERT=$MAX_INSERT
TRUTH=$TRUTH
THRESHOLD_MIN=$THRESHOLD_MIN
SHORT_NAME=$SHORT_NAME

EOF

cat <<EOF

hadoop jar $BUILD_DIR/target/cloudbreak-1.0-SNAPSHOT-exe.jar -Dmapred.reduce.tasks=200 GMMFitSingleEndInsertSizes
   --inputFileDescriptor $READ_GROUP_FILE 
   --outputHDFSDir $HDFS_SAMPLE_DIR/$NAME 
   --faidx $HDFS_FAI 
   --maxInsertSize $MAX_INSERT 
   --mapabilityWeighting $MAPABILITY 
   --excludePairsMappingIn $FILTER 
   --resolution $RESOLUTION 
   --aligner $ALIGNER

EOF

hadoop jar $BUILD_DIR/target/cloudbreak-1.0-SNAPSHOT-exe.jar -Dmapred.reduce.tasks=200 GMMFitSingleEndInsertSizes \
    --inputFileDescriptor $READ_GROUP_FILE \
    --outputHDFSDir $HDFS_SAMPLE_DIR/$NAME \
    --faidx $HDFS_FAI \
    --maxInsertSize $MAX_INSERT \
    --mapabilityWeighting $MAPABILITY \
    --excludePairsMappingIn $FILTER \
    --resolution $RESOLUTION \
    --aligner $ALIGNER

cat <<EOF

hadoop jar $BUILD_DIR/target/cloudbreak-1.0-SNAPSHOT-exe.jar exportWigAndBedFiles 
   --inputHDFSDir $HDFS_SAMPLE_DIR/$NAME 
   --faidx $LOCAL_FAI 
   --resolution $RESOLUTION --medianFilterWindow $MEDIAN_FILTER_WINDOW --outputPrefix $NAME

EOF

hadoop jar $BUILD_DIR/target/cloudbreak-1.0-SNAPSHOT-exe.jar exportGMMResults \
    --inputHDFSDir $HDFS_SAMPLE_DIR/$NAME \
    --faidx $LOCAL_FAI \
    --resolution $RESOLUTION --medianFilterWindow $MEDIAN_FILTER_WINDOW --outputPrefix $NAME

echo gzip *.wig
gzip *.wig

~/software/IGVTools/igvtools tile $NAME_w0.wig.gz $NAME_w0.wig.gz.tdf hg18
~/software/IGVTools/igvtools tile $NAME_mu1.wig.gz $NAME_mu1.wig.gz.tdf hg18
~/software/IGVTools/igvtools tile $NAME_l1.wig.gz $NAME_l1.wig.gz.tdf hg18
~/software/IGVTools/igvtools tile $NAME_l2.wig.gz $NAME_l2.wig.gz.tdf hg18
~/software/IGVTools/igvtools tile $NAME_l1f.wig.gz $NAME_l1f.wig.gz.tdf hg18
~/software/IGVTools/igvtools tile $NAME_lrHet.wig.gz $NAME_lrHet.wig.gz.tdf hg18
~/software/IGVTools/igvtools tile $NAME_lrHom.wig.gz $NAME_lrHom.wig.gz.tdf hg18


