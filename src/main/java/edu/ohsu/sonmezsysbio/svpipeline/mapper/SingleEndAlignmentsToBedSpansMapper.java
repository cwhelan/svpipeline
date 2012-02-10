package edu.ohsu.sonmezsysbio.svpipeline.mapper;

import edu.ohsu.sonmezsysbio.svpipeline.NovoalignNativeRecord;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;

import java.io.IOException;
import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/23/11
 * Time: 10:12 AM
 */
public class SingleEndAlignmentsToBedSpansMapper extends MapReduceBase implements Mapper<LongWritable, Text, Text, Text> {

    private boolean matePairs;

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

        matePairs = Boolean.parseBoolean(job.get("pileupDeletionScore.isMatePairs"));

    }

    public void map(LongWritable key, Text value, OutputCollector<Text, Text> output, Reporter reporter)
            throws IOException {
        String line = value.toString();
        //System.err.println("LINE: " + line);

        String[] lineParts = line.split("\tSEP\t");
        String[] fields1 = lineParts[0].split("\t");
        String[] alnFields1 = Arrays.copyOfRange(fields1, 1, fields1.length);

        NovoalignNativeRecord record1 = NovoalignNativeRecord.parseRecord(alnFields1);
        String[] alnFields2 = lineParts[1].split("\t");
        NovoalignNativeRecord record2 = NovoalignNativeRecord.parseRecord(alnFields2);

        // todo: not handling translocations for now
        if (! record1.getChromosomeName().equals(record2.getChromosomeName())) return;

        int endPosterior1 = 0;
        int endPosterior2 = 0;
        try {
            endPosterior1 = record1.getPosteriorProb();
            endPosterior2 = record2.getPosteriorProb();
            //System.err.println("posteriors: " + endPosterior1 + "," + endPosterior2);
        } catch (NumberFormatException e) {
            System.err.println("Problem with line: " + line);
            e.printStackTrace();
            throw new RuntimeException("Could not parse input record2 " + line);
        }

//        boolean matePair = false;
//        if (record1.isForward() && ! record2.isForward()) {
//            if (record1.getPosition() - record2.getPosition() > 0) matePair = true;
//            if (record1.getPosition() - record2.getPosition() < 0 &&
//                    record1.getPosition() - record2.getPosition() > -500) matePair = false;
//        } else if (!record1.isForward() && record2.isForward()) {
//            if (record1.getPosition() - record2.getPosition() < 0) matePair = true;
//            if (record1.getPosition() - record2.getPosition() > 0 &&
//                    record1.getPosition() - record2.getPosition() < 500) matePair = false;
//        } else {
//            return;
//        }

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

        if (insertSize > 1000000) {
            System.err.println("Pair " + fields1[0]  + ": Insert size would be greater than 100,000 - skipping");
            return;
        }


        double deletionScore = computeDeletionScore(
                endPosterior1,
                endPosterior2,
                insertSize,
                isizeMean,
                isizeSD
        );

        output.collect(new Text(leftRead.getReadId()),
                new Text(leftRead.getChromosomeName() + "\t" + leftRead.getPosition() + "\t" + rightRead.getPosition() + "\t" + leftRead.getReadId() + "\t" + deletionScore));

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
