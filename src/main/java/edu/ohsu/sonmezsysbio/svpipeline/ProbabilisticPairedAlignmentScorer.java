package edu.ohsu.sonmezsysbio.svpipeline;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 3/12/12
 * Time: 3:11 PM
 */
public class ProbabilisticPairedAlignmentScorer extends PairedAlignmentScorer {
    @Override
    public double computeDeletionScore(int insertSize, Double targetIsize, Double targetIsizeSD, Double pMappingCorrect) {
        // calculate:
        //
        // p(Del | IS) = p( IS | Del) p(Del) / p (IS)
        //
        // p(NoDel | IS) = p( IS | NoDel) p(NoDel) / p (IS)
        //
        // ((p(IS | Del) p (Del)) / (p(IS | NoDel) p(NoDel))) (((( * p (IS | AQ) ???? )))
        //
        // log ( ((p(IS | Del) p (Del)) / (p(IS | NoDel) p(NoDel))) * p (IS | AQ) )
        // log ((p(IS | Del) p (Del)) - log (p(IS | NoDel) p(NoDel)) + log ( p (IS | AQ))
        //
        // p (Del) = # of deletions per individual / genome size = 2432 / 3137161264
        // p (IS | Del) = density at IS in gamma(0.2, 40000)
        // p (noDel) =  1 - p(Del)
        // p (IS | NoDel) = density at IS in N(isize, isizeSD)

        //GammaDistribution gammaDistribution = new GammaDistribution(0.2, 40000);
        LogNormalDistribution logNormalDistribution = new LogNormalDistribution(6, 0.6);
        NormalDistribution normalDistribution = new NormalDistribution(targetIsize, targetIsizeSD);

        double pDeletion = Math.log(2432.0 / 2700000000.0);
        double pNoDeletion = Math.log(1 - 2432.0 / 2700000000.0);

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

        double pDeletionNew = logAdd(pDeletionGivenIS + pMappingCorrect, pDeletion + pMappingIncorrect);
        double pNoDeletionNew = logAdd(pNoDeletionGivenIS + pMappingCorrect, pNoDeletion + pMappingIncorrect);

        double likelihoodRatio = pDeletionNew
                                - pNoDeletionNew;

//        System.out.println("-----");
//        System.out.println("pMappingCorrect:\t" + pMappingCorrect);
//        System.out.println("pMappingIncorrect:\t" + pMappingIncorrect);
//        System.out.println("pISgivenDeletion:\t" + pISgivenDeletion);
//        System.out.println("pISgivenNoDeletion:\t" + pISgivenNoDeletion);
//        System.out.println("pDeletionGivenIS:\t" + pDeletionGivenIS);
//        System.out.println("pNoDeletionGivenIS:\t" + pNoDeletionGivenIS);
//        System.out.println("pDeletion:\t" + pDeletion);
//        System.out.println("pNoDeletion:\t" + pNoDeletion);
//
//        System.out.println("pDeletionGivenIS + pMappingCorrect:\t" + (pDeletionGivenIS + pMappingCorrect));
//        System.out.println("pNoDeletionGivenIS + pMappingCorrect:\t" + (pNoDeletionGivenIS + pMappingCorrect));
//
//        System.out.println("pDeletionNew:\t" + pDeletionNew);
//        System.out.println("pNoDeletionNew:\t" + pNoDeletionNew);
//        System.out.println("likelihoodRatio:\t" + likelihoodRatio);
//        System.out.println("-----");

        return likelihoodRatio;
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
