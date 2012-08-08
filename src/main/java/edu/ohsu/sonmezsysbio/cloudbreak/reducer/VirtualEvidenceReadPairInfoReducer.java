package edu.ohsu.sonmezsysbio.cloudbreak.reducer;

import edu.ohsu.sonmezsysbio.cloudbreak.ReadGroupInfo;
import edu.ohsu.sonmezsysbio.cloudbreak.file.ReadGroupInfoFileHelper;
import edu.ohsu.sonmezsysbio.cloudbreak.io.ReadPairInfo;
import edu.ohsu.sonmezsysbio.svpipeline.io.GenomicLocation;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.mapred.*;

import java.io.IOException;
import java.util.Iterator;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/6/12
 * Time: 1:35 PM
 */
public class VirtualEvidenceReadPairInfoReducer extends MapReduceBase implements Reducer<GenomicLocation, ReadPairInfo, GenomicLocation, DoubleWritable> {


    private Map<Short,ReadGroupInfo> readGroupInfos;

    public Map<Short, ReadGroupInfo> getReadGroupInfos() {
        return readGroupInfos;
    }

    public void setReadGroupInfos(Map<Short, ReadGroupInfo> readGroupInfos) {
        this.readGroupInfos = readGroupInfos;
    }

    public void reduce(GenomicLocation key, Iterator<ReadPairInfo> values, OutputCollector<GenomicLocation, DoubleWritable> output, Reporter reporter) throws IOException {
        System.err.println(key);
        LogNormalDistribution logNormalDistribution = new LogNormalDistribution(6, 0.6);

        double deletionPrior = Math.log(2432.0 / 2700000000.0);
        double noDeletionPrior = Math.log(1 - 2432.0 / 2700000000.0);

        System.err.println("del prior: " + deletionPrior);
        System.err.println("no del prior: " + noDeletionPrior);

        double pDeletionFactorProduct = deletionPrior;
        double pNoDeletionFactorProduct = deletionPrior;

        while (values.hasNext()) {
            ReadPairInfo readPairInfo = values.next();
            int insertSize = readPairInfo.insertSize;
            System.err.println("insert size: " + insertSize);
            double pMappingCorrect = readPairInfo.pMappingCorrect;
            System.err.println("p mapping correct: " + pMappingCorrect);
            short readGroupId = readPairInfo.readGroupId;

            ReadGroupInfo readGroupInfo = readGroupInfos.get(readGroupId);
            int targetIsize = readGroupInfo.isize;
            int targetIsizeSD = readGroupInfo.isizeSD;

            NormalDistribution fragmentSizeDistribution = new NormalDistribution(targetIsize, targetIsizeSD);
            double pISgivenDeletion = Math.log(logNormalDistribution.density(insertSize));         // todo add fragment size
            System.err.println("pISgivenDeletion: " + pISgivenDeletion);
            double pISgivenNoDeletion = Math.log(fragmentSizeDistribution.density(insertSize));
            // todo
            // need to cap p(IS | NoDel) because it goes to infinity as the insert size gets large
            if (insertSize > targetIsize + 30 * targetIsizeSD) {
                pISgivenNoDeletion = Math.log(fragmentSizeDistribution.density(targetIsize + 30 * targetIsizeSD));
            }
            System.err.println("pISgivenNoDeletion: " + pISgivenNoDeletion);

            double pMappingIncorrect = Math.log(1 - Math.exp(pMappingCorrect));

            double normalization = logAdd(deletionPrior + pISgivenDeletion, noDeletionPrior + pISgivenNoDeletion);
            double pDeletionGivenIS = deletionPrior + pISgivenDeletion - normalization;
            System.err.println("pDeletionGivenIS: " + pDeletionGivenIS);

            double pNoDeletionGivenIS = noDeletionPrior + pISgivenNoDeletion - normalization;
            System.err.println("pNoDeletionGivenIS: " + pNoDeletionGivenIS);

            double pDeletion = logAdd(pDeletionGivenIS + pMappingCorrect, deletionPrior + pMappingIncorrect);
            double pNoDeletion = logAdd(pNoDeletionGivenIS + pMappingCorrect, noDeletionPrior + pMappingIncorrect);

            pDeletionFactorProduct = pDeletionFactorProduct + pDeletion;
            pNoDeletionFactorProduct = pNoDeletionFactorProduct + pNoDeletion;
        }

        double lr = pDeletionFactorProduct - pNoDeletionFactorProduct;
        System.err.println("lr: " + lr);
        output.collect(key, new DoubleWritable(lr));
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
        return logX + Math.log(1.0 + Math.exp(negDiff));
    }

    @Override
    public void configure(JobConf job) {
        super.configure(job);

        String readGroupInfoFile = job.get("read.group.info.file");
        ReadGroupInfoFileHelper readGroupInfoFileHelper = new ReadGroupInfoFileHelper();
        try {
            readGroupInfos = readGroupInfoFileHelper.readReadGroupsById(readGroupInfoFile);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
