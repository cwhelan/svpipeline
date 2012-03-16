package edu.ohsu.sonmezsysbio.svpipeline.mapper;

import edu.ohsu.sonmezsysbio.svpipeline.SAMRecord;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.hadoop.io.DoubleWritable;
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
public class MappedPairToDeletionScoreMapper extends MapReduceBase implements Mapper<LongWritable, Text, Text, DoubleWritable> {

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

    }

    public void map(LongWritable key, Text value, OutputCollector<Text, DoubleWritable> output, Reporter reporter)
            throws IOException {
        String line = value.toString();

        String[] lineParts = line.split("\tSEP\t");
        String[] fields1 = lineParts[0].split("\t");
        String[] samFields1 = Arrays.copyOfRange(fields1, 1, fields1.length);

        SAMRecord samRecord1 = SAMRecord.parseSamRecord(samFields1);
        String[] samFields2 = lineParts[1].split("\t");
        SAMRecord samRecord2 = SAMRecord.parseSamRecord(samFields2);

        if (! (samRecord1.isPairMapped())) {
            return;
        }

        if (samRecord1.isInterChromosomal()) {
            return;
        }

        int endPosterior1 = 0;
        int endPosterior2 = 0;
        try {
            endPosterior1 = samRecord1.getEndPosterior();
            endPosterior2 = samRecord2.getEndPosterior();
        } catch (NumberFormatException e) {
            System.err.println("Problem with line: " + line);
            e.printStackTrace();
            throw new RuntimeException("Could not parse input record " + line);
        }

        SAMRecord leftRead = samRecord1.getPosition() < samRecord2.getPosition() ?
                samRecord1 : samRecord2;
        SAMRecord rightRead = samRecord1.getPosition() < samRecord2.getPosition() ?
                samRecord2 : samRecord1;


        //int insertSize = leftRead.getInsertSize();

        int insertSize = rightRead.getPosition() - leftRead.getPosition();

        double deletionScore = computeDeletionScore(endPosterior1, endPosterior2, insertSize, targetIsize, targetIsizeSD, samRecord1.isProperPair());

        int genomeOffset = leftRead.getPosition() - leftRead.getPosition() % 500;

        insertSize = insertSize + leftRead.getPosition() % 500 + 500 - rightRead.getPosition() % 500;

        if (insertSize > 1000000) {
            System.err.println("Pair " + fields1[0]  + ": Insert size would be greater than 100,000 - skipping");
            return;
        }

        for (int i = 0; i <= insertSize; i = i + 500) {
            Text outKey = new Text(samRecord1.getReferenceName() + "\t" + (genomeOffset + i));
            DoubleWritable outVal = new DoubleWritable(deletionScore);
            output.collect(outKey, outVal);
        }
    }

    public static double computeDeletionScore(double endPosterior1, double endPosterior2, int insertSize, Double targetIsize, Double targetIsizeSD, boolean properPair) {
        NormalDistribution insertSizeDist = new NormalDistribution(targetIsize, targetIsizeSD);
        // deletion score = endPosterior1 * endPosterior2 * P(X < insertSize - 2 * targetIsizeSD)
        double deletionProb;
        deletionProb = insertSizeDist.cumulativeProbability(insertSize - 2 * targetIsizeSD);

        return (1 - Math.pow(10, endPosterior1 / -10.0)) *
                (1 - Math.pow(10.0, endPosterior2 / -10.0)) *
                (properPair ? 1.0 : 1e-7) *
                deletionProb;

    }

}
