package edu.ohsu.sonmezsysbio.svpipeline.mapper;

import edu.ohsu.sonmezsysbio.svpipeline.NovoalignNativeRecord;
import edu.ohsu.sonmezsysbio.svpipeline.SVPipeline;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/23/11
 * Time: 10:12 AM
 */
public class SingleEndAlignmentsToDeletionScoreMapper extends MapReduceBase implements Mapper<LongWritable, Text, Text, DoubleWritable> {

    public static final int WINDOW_SIZE = 100;
    private Integer maxInsertSize = 500000;
    private boolean matePairs;

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

    public void map(LongWritable key, Text value, OutputCollector<Text, DoubleWritable> output, Reporter reporter)
            throws IOException {
        String line = value.toString();
        int firstTabIndex = line.indexOf('\t');
        String readPairId = line.substring(0,firstTabIndex);
        String lineValues = line.substring(firstTabIndex + 1);
        
        String[] readAligments = lineValues.split(SVPipeline.READ_SEPARATOR);
        String read1AlignmentsString = readAligments[0];
        String[] read1Alignments = read1AlignmentsString.split(SVPipeline.ALIGMENT_SEPARATOR);
        List<NovoalignNativeRecord> read1AlignmentRecords = parseAlignmentsIntoRecords(read1Alignments);

        String read2AlignmentsString = readAligments[1];
        String[] read2Alignments = read2AlignmentsString.split(SVPipeline.ALIGMENT_SEPARATOR);
        List<NovoalignNativeRecord> read2AlignmentRecords = parseAlignmentsIntoRecords(read2Alignments);

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

        if (matePairs) {
            boolean matePair = false;

            // validate mapping orientations
            if (record1.isForward() && ! record2.isForward()) {
                if (record1.getPosition() - record2.getPosition() > 0) matePair = true;
                if (record1.getPosition() - record2.getPosition() < 0 &&
                        record1.getPosition() - record2.getPosition() > -500) matePair = false;
            } else if (!record1.isForward() && record2.isForward()) {
                if (record1.getPosition() - record2.getPosition() < 0) matePair = true;
                if (record1.getPosition() - record2.getPosition() > 0 &&
                        record1.getPosition() - record2.getPosition() < 500) matePair = false;
            } else {
                return;
            }

            insertSize = rightRead.getPosition() - leftRead.getPosition();
            //System.err.println("insert size: " + insertSize);

            isizeMean = matePair ? targetIsize : 150;
            isizeSD = matePair ? targetIsizeSD : 15;
        } else {

            // validate mapping orientations
            if (record1.isForward() && ! record2.isForward()) {
                if (record1.getPosition() > record2.getPosition()) return;
            } else if (!record1.isForward() && record2.isForward()) {
                if (record1.getPosition() < record2.getPosition()) return;
            } else {
                return;
            }

            insertSize = rightRead.getPosition() - leftRead.getPosition();
            isizeMean = targetIsize;
            isizeSD = targetIsizeSD;
        }

        if (insertSize > maxInsertSize) {
            System.err.println("Pair " + record1.getReadId()  + ": Insert size would be greater than 100,000 - skipping");
            return;
        }


        double deletionScore = computeDeletionScore(
                endPosterior1,
                endPosterior2,
                insertSize,
                isizeMean,
                isizeSD
        );
        //System.err.println("computed deletion score : " + deletionScore);

        int genomeOffset = leftRead.getPosition() - leftRead.getPosition() % WINDOW_SIZE;

        insertSize = insertSize + leftRead.getPosition() % WINDOW_SIZE + WINDOW_SIZE - rightRead.getPosition() % WINDOW_SIZE;

        for (int i = 0; i <= insertSize; i = i + WINDOW_SIZE) {
            Text outKey = new Text(record1.getChromosomeName() + "\t" + (genomeOffset + i));
            DoubleWritable outVal = new DoubleWritable(deletionScore);
            output.collect(outKey, outVal);
        }

    }

    private List<NovoalignNativeRecord> parseAlignmentsIntoRecords(String[] read1Alignments) {
        List<NovoalignNativeRecord> read1AlignmentList = new ArrayList<NovoalignNativeRecord>();
        for (String alignment : read1Alignments) {
            String[] fields1 = alignment.split("\t");
            NovoalignNativeRecord record1 = NovoalignNativeRecord.parseRecord(fields1);
            read1AlignmentList.add(record1);
        }
        return read1AlignmentList;
    }

    public static double computeDeletionScore(int codedEndPosterior1, int codedEndPosterior2, int insertSize, Double targetIsize, Double targetIsizeSD) {
        //System.err.println("target isize: " + targetIsize + ", sd " + targetIsizeSD);
        NormalDistribution insertSizeDist = new NormalDistributionImpl(targetIsize, targetIsizeSD);
        // deletion score = codedEndPosterior1 * codedEndPosterior2 * P(X < insertSize - 2 * targetIsizeSD)

        double deletionProb;
        try {
            deletionProb = insertSizeDist.cumulativeProbability(Math.max(0, insertSize - 1.5 * targetIsizeSD));
            //System.err.println("Deletion prob: " + deletionProb);
        } catch (MathException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

        double vote = deletionProb > 0.5 ? 1 : -1;


        double endPosterior1 = codedEndPosterior1 == 0 ? 0.0001 : 1 - Math.pow(10.0, codedEndPosterior1 / -10.0);
        double endPosterior2 = codedEndPosterior2 == 0 ? 0.0001 : 1 - Math.pow(10.0, codedEndPosterior2 / -10.0);
        //System.err.println("posteriors: " + endPosterior1 + "," + endPosterior2);

        //return deletionProb + endPosterior1 + endPosterior2;
        return vote * endPosterior1 * endPosterior2;
    }

}
