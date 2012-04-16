package edu.ohsu.sonmezsysbio.svpipeline.mapper;

import edu.ohsu.sonmezsysbio.svpipeline.*;
import org.apache.hadoop.io.DoubleWritable;
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
public class SingleEndAlignmentsToDeletionScoreMapper extends MapReduceBase implements Mapper<LongWritable, Text, Text, DoubleWritable> {

    private boolean matePairs;
    private Integer maxInsertSize = 500000;
    private PairedAlignmentScorer scorer;

    public Integer getMaxInsertSize() {
        return maxInsertSize;
    }

    public void setMaxInsertSize(Integer maxInsertSize) {
        this.maxInsertSize = maxInsertSize;
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

    public boolean isMatePairs() {
        return matePairs;
    }

    public void setMatePairs(boolean matePairs) {
        this.matePairs = matePairs;
    }

    public PairedAlignmentScorer getScorer() {
        return scorer;
    }

    public void setScorer(PairedAlignmentScorer scorer) {
        this.scorer = scorer;
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

        scorer = new ProbabilisticPairedAlignmentScorer();
    }

    public void map(LongWritable key, Text value, OutputCollector<Text, DoubleWritable> output, Reporter reporter)
            throws IOException {
        String line = value.toString();
        int firstTabIndex = line.indexOf('\t');
        String readPairId = line.substring(0,firstTabIndex);
        String lineValues = line.substring(firstTabIndex + 1);
        
        String[] readAligments = lineValues.split(SVPipeline.READ_SEPARATOR);
        String read1AlignmentsString = readAligments[0];
        String[] read1Alignments = read1AlignmentsString.split(SVPipeline.ALIGNMENT_SEPARATOR);
        List<NovoalignNativeRecord> read1AlignmentRecords = NovoalignSingleEndMapperHelper.parseAlignmentsIntoRecords(read1Alignments);

        String read2AlignmentsString = readAligments[1];
        String[] read2Alignments = read2AlignmentsString.split(SVPipeline.ALIGNMENT_SEPARATOR);
        List<NovoalignNativeRecord> read2AlignmentRecords = NovoalignSingleEndMapperHelper.parseAlignmentsIntoRecords(read2Alignments);

        emitDeletionScoresForAllPairs(read1AlignmentRecords, read2AlignmentRecords, output);
    }

    private void emitDeletionScoresForAllPairs(List<NovoalignNativeRecord> read1AlignmentRecords, List<NovoalignNativeRecord> read2AlignmentRecords, OutputCollector<Text, DoubleWritable> output) throws IOException {
        for (NovoalignNativeRecord record1 : read1AlignmentRecords) {
            for (NovoalignNativeRecord record2 : read2AlignmentRecords) {
                emitDeletionScoresForPair(record1, record2, output);
            }
        }
    }

    private void emitDeletionScoresForPair(NovoalignNativeRecord record1, NovoalignNativeRecord record2, OutputCollector<Text, DoubleWritable> output) throws IOException {

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
        if (!scorer.validateMappingOrientations(record1, record2, matePairs)) return;

        isizeMean = targetIsize;
        isizeSD = targetIsizeSD;
        if (matePairs) {
            //System.err.println("insert size: " + insertSize);
            if (!scorer.isMatePairNotSmallFragment(record1, record2)) {
                isizeMean = 150.0;
                isizeSD = 15.0;
            }
        }

        insertSize = rightRead.getPosition() + rightRead.getSequence().length() - leftRead.getPosition();

        if (! scorer.validateInsertSize(insertSize, record1.getReadId(), maxInsertSize)) return;

        Double pMappingCorrect = scorer.probabilityMappingIsCorrect(endPosterior1, endPosterior2);

        double deletionScore = scorer.computeDeletionScore(
                insertSize,
                isizeMean,
                isizeSD,
                pMappingCorrect
        );
        //System.err.println("computed deletion score : " + deletionScore);

        int genomeOffset = leftRead.getPosition() - leftRead.getPosition() % SVPipeline.RESOLUTION;

        insertSize = insertSize + leftRead.getPosition() % SVPipeline.RESOLUTION + SVPipeline.RESOLUTION - rightRead.getPosition() % SVPipeline.RESOLUTION;

        for (int i = 0; i <= insertSize; i = i + SVPipeline.RESOLUTION) {
            Text outKey = new Text(record1.getChromosomeName() + "\t" + (genomeOffset + i));
            DoubleWritable outVal = new DoubleWritable(deletionScore);
            output.collect(outKey, outVal);
        }

    }

}
