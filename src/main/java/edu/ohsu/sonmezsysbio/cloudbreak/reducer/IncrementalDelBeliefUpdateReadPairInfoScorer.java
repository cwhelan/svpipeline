package edu.ohsu.sonmezsysbio.cloudbreak.reducer;

import edu.ohsu.sonmezsysbio.cloudbreak.ReadGroupInfo;
import edu.ohsu.sonmezsysbio.cloudbreak.io.ReadPairInfo;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Iterator;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 8/22/12
 * Time: 11:14 AM
 */
public class IncrementalDelBeliefUpdateReadPairInfoScorer implements ReadPairInfoScorer {
    public double reduceReadPairInfos(Iterator<ReadPairInfo> values, Map<Short, ReadGroupInfo> readGroupInfos) {
        LogNormalDistribution logNormalDistribution = new LogNormalDistribution(6, 0.6);

        double pDeletion = Math.log(2432.0 / 2700000000.0);
        double pNoDeletion = Math.log(1 - 2432.0 / 2700000000.0);


        while (values.hasNext()) {
            ReadPairInfo readPairInfo = values.next();
            int insertSize = readPairInfo.insertSize;
            double pMappingCorrect = readPairInfo.pMappingCorrect;
            short readGroupId = readPairInfo.readGroupId;

            ReadGroupInfo readGroupInfo = readGroupInfos.get(readGroupId);
            int targetIsize = readGroupInfo.isize;
            int targetIsizeSD = readGroupInfo.isizeSD;
            boolean matePairs = readGroupInfo.matePair;

            NormalDistribution normalDistribution = new NormalDistribution(targetIsize, targetIsizeSD);
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
        return pDeletion - pNoDeletion;
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

}