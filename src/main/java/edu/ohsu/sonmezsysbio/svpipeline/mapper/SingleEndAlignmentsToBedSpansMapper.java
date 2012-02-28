package edu.ohsu.sonmezsysbio.svpipeline.mapper;

import edu.ohsu.sonmezsysbio.svpipeline.NovoalignNativeRecord;
import edu.ohsu.sonmezsysbio.svpipeline.SVPipeline;
import edu.ohsu.sonmezsysbio.svpipeline.SingleEndAlignmentScorer;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;

import java.io.IOException;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/23/11
 * Time: 10:12 AM
 */
public class SingleEndAlignmentsToBedSpansMapper extends MapReduceBase implements Mapper<LongWritable, Text, Text, Text> {

    private boolean matePairs;
    private Integer maxInsertSize = 500000;

    public Integer getMaxInsertSize() {
        return maxInsertSize;
    }

    public void setMaxInsertSize(Integer maxInsertSize) {
        this.maxInsertSize = maxInsertSize;
    }

    public boolean isMatePairs() {
        return matePairs;
    }

    public void setMatePairs(boolean matePairs) {
        this.matePairs = matePairs;
    }


    public Double getTargetIsize() {
        return targetIsize;
    }

    public void setTargetIsize(Double targetIsize) {
        this.targetIsize = targetIsize;
    }

    public Double getTargetIsizeSD() {
        return targetIsizeSD;
    }

    public void setTargetIsizeSD(Double targetIsizeSD) {
        this.targetIsizeSD = targetIsizeSD;
    }

    private Double targetIsize;
    private Double targetIsizeSD;

    @Override
    public void configure(JobConf job) {
        super.configure(job);
        targetIsize = Double.parseDouble(job.get("pileupDeletionScore.targetIsize"));
        targetIsizeSD = Double.parseDouble(job.get("pileupDeletionScore.targetIsizeSD"));

        maxInsertSize = Integer.parseInt(job.get("pileupDeletionScore.maxInsertSize"));
        matePairs = Boolean.parseBoolean(job.get("pileupDeletionScore.isMatePairs"));
    }

    public void map(LongWritable key, Text value, OutputCollector<Text, Text> output, Reporter reporter)
            throws IOException {

        String line = value.toString();
        //System.err.println("LINE: " + line);

        int firstTabIndex = line.indexOf('\t');
        String readPairId = line.substring(0,firstTabIndex);
        String lineValues = line.substring(firstTabIndex + 1);

        String[] readAligments = lineValues.split(SVPipeline.READ_SEPARATOR);
        String read1AlignmentsString = readAligments[0];
        String[] read1Alignments = read1AlignmentsString.split(SVPipeline.ALIGMENT_SEPARATOR);
        List<NovoalignNativeRecord> read1AlignmentRecords = NovoalignSingleEndMapperHelper.parseAlignmentsIntoRecords(read1Alignments);

        String read2AlignmentsString = readAligments[1];
        String[] read2Alignments = read2AlignmentsString.split(SVPipeline.ALIGMENT_SEPARATOR);
        List<NovoalignNativeRecord> read2AlignmentRecords = NovoalignSingleEndMapperHelper.parseAlignmentsIntoRecords(read2Alignments);

        emitDeletionScoresForAllPairs(read1AlignmentRecords, read2AlignmentRecords, output);
    }

    private void emitDeletionScoresForAllPairs(List<NovoalignNativeRecord> read1AlignmentRecords, List<NovoalignNativeRecord> read2AlignmentRecords, OutputCollector<Text, Text> output) throws IOException {
        for (NovoalignNativeRecord record1 : read1AlignmentRecords) {
            for (NovoalignNativeRecord record2 : read2AlignmentRecords) {
                emitBedSpanForPair(record1, record2, output);
            }
        }
    }

    public void emitBedSpanForPair(NovoalignNativeRecord record1, NovoalignNativeRecord record2, OutputCollector<Text, Text> output) throws IOException {

        // todo: not handling translocations for now
        if (! record1.getChromosomeName().equals(record2.getChromosomeName())) return;

        int endPosterior1 = record1.getPosteriorProb();
        int endPosterior2 = record2.getPosteriorProb();

        int insertSize = -1;
        Double isizeMean;
        Double isizeSD;

        NovoalignNativeRecord leftRead = record1.getPosition() < record2.getPosition() ?
                record1 : record2;
        NovoalignNativeRecord rightRead = record1.getPosition() < record2.getPosition() ?
                record2 : record1;

        // todo: not handling inversions for now
        if (!SingleEndAlignmentScorer.validateMappingOrientations(record1, record2, isMatePairs())) return;

        isizeMean = targetIsize;
        isizeSD = targetIsizeSD;
        if (matePairs) {
            if (!SingleEndAlignmentScorer.isMatePairNotSmallFragment(record1, record2)) {
                isizeMean = 150.0;
                isizeSD = 15.0;
            }
        }

        insertSize = rightRead.getPosition() - leftRead.getPosition();

        if (SingleEndAlignmentScorer.validateInsertSize(insertSize, record1.getReadId(), maxInsertSize)) return;

        double deletionScore = SingleEndAlignmentScorer.computeDeletionScore(
                endPosterior1,
                endPosterior2,
                insertSize,
                isizeMean,
                isizeSD
        );

        output.collect(new Text(leftRead.getReadId()),
                new Text(leftRead.getChromosomeName() + "\t" + leftRead.getPosition() + "\t" + rightRead.getPosition() + "\t" + leftRead.getReadId() + "\t" + deletionScore));

    }
}
