package edu.ohsu.sonmezsysbio.svpipeline.reducer;

import edu.ohsu.sonmezsysbio.svpipeline.io.GenomicLocation;
import edu.ohsu.sonmezsysbio.svpipeline.io.ReadPairInfo;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/6/12
 * Time: 1:35 PM
 */
public class IncrementalDelBeliefUpdateReadPairInfoReducer extends MapReduceBase implements Reducer<GenomicLocation, ReadPairInfo, Text, DoubleWritable> {

    private double targetIsize;
    private double targetIsizeSD;
    private boolean matePairs;
    private String faidxFileName;
    private Map<Short, String> chromosomesByKey;

    public void reduce(GenomicLocation key, Iterator<ReadPairInfo> values, OutputCollector<Text, DoubleWritable> output, Reporter reporter) throws IOException {

        LogNormalDistribution logNormalDistribution = new LogNormalDistribution(6, 0.6);
        NormalDistribution normalDistribution = new NormalDistribution(targetIsize, targetIsizeSD);

        double pDeletion = Math.log(2432.0 / 2700000000.0);
        double pNoDeletion = Math.log(1 - 2432.0 / 2700000000.0);


        while (values.hasNext()) {
            ReadPairInfo readPairInfo = values.next();
            int insertSize = readPairInfo.insertSize;
            double pMappingCorrect = readPairInfo.pMappingCorrect;

            double pISgivenDeletion = Math.log(logNormalDistribution.density(insertSize));         // todo add fragment size
            double pISgivenNoDeletion = Math.log(normalDistribution.density(insertSize));
            // todo
            // need to cap p(IS | NoDel) because it goes to infinity as the insert size gets large
            if (insertSize > targetIsize + 30 * targetIsizeSD) {
                pISgivenNoDeletion = Math.log(normalDistribution.density(targetIsize + 30 * targetIsizeSD));
            }

            double pMappingIncorrect = Math.log(1 - Math.exp(pMappingCorrect));

            double normalization = logAdd(pDeletion + pISgivenDeletion, pNoDeletion + pISgivenNoDeletion);
            double pDeletionGivenIS = pDeletion + pISgivenDeletion - normalization;
            double pNoDeletionGivenIS = pNoDeletion + pISgivenNoDeletion - normalization;

            pDeletion = logAdd(pDeletionGivenIS + pMappingCorrect, pDeletion + pMappingIncorrect);
            pNoDeletion = logAdd(pNoDeletionGivenIS + pMappingCorrect, pNoDeletion + pMappingIncorrect);

        }
        double lr = pDeletion / pNoDeletion;
        output.collect(new Text(chromosomesByKey.get(key.chromosome) + "\t" + key.pos), new DoubleWritable(lr));
    }

    // from https://facwiki.cs.byu.edu/nlp/index.php/Log_Domain_Computations
    public static double logAdd(double logX, double logY) {
        // 1. make X the max
        if (logY > logX) {
            double temp = logX;
            logX = logY;
            logY = temp;
        }
        // 2. now X is bigger
        if (logX == Double.NEGATIVE_INFINITY) {
            return logX;
        }
        // 3. how far "down" (think decibels) is logY from logX?
        //    if it's really small (20 orders of magnitude smaller), then ignore
        double negDiff = logY - logX;
        if (negDiff < -20) {
            return logX;
        }
        // 4. otherwise use some nice algebra to stay in the log domain
        //    (except for negDiff)
        return logX + java.lang.Math.log(1.0 + java.lang.Math.exp(negDiff));
    }

    @Override
    public void configure(JobConf job) {
        super.configure(job);
        targetIsize = Double.parseDouble(job.get("pileupDeletionScore.targetIsize"));
        targetIsizeSD = Double.parseDouble(job.get("pileupDeletionScore.targetIsizeSD"));

        matePairs = Boolean.parseBoolean(job.get("pileupDeletionScore.isMatePairs"));
        faidxFileName = job.get("alignment.faidx");
        try {
            chromosomesByKey = readFaidxFile(faidxFileName);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    private Map<Short, String> readFaidxFile(String faidxFileName) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(faidxFileName));
        return readFaidx(reader);
    }

    public Map<Short, String > readFaidx(BufferedReader bufferedReader) throws IOException {
        String line;
        short chrIdx = 0;
        Map<Short, String> chrTable = new HashMap<Short, String>();
        while ((line = bufferedReader.readLine()) != null) {
            String chromosomeName = line.split("\\s+")[0];
            chrTable.put(chrIdx, chromosomeName);
            chrIdx++;
        }
        return chrTable;
    }

}
