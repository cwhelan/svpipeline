package edu.ohsu.sonmezsysbio.svpipeline;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 3/12/12
 * Time: 3:11 PM
 */
public class ProbabilisticPairedAlignmentScorer extends PairedAlignmentScorer {
    @Override
    public double computeDeletionScore(int codedEndPosterior1, int codedEndPosterior2, int insertSize, Double targetIsize, Double targetIsizeSD) {
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

        GammaDistribution gammaDistribution = new GammaDistribution(0.2, 40000);
        NormalDistribution normalDistribution = new NormalDistribution(targetIsize, targetIsizeSD);

        double pDeletion = Math.log(2432.0 / 3137161264.0);
        double pISgivenDeletion = Math.log(gammaDistribution.density(insertSize));         // todo add fragment size

        double pNoDeletion = Math.log(1 - 2432.0 / 3137161264.0);
        double pISgivenNoDeletion = Math.log(normalDistribution.density(insertSize));

        // todo
        // need to cap p(IS | NoDel) because it goes to infinity as the insert size gets large
        if (insertSize > targetIsize + 10 * targetIsizeSD) {
            pISgivenNoDeletion = Math.log(normalDistribution.density(targetIsize + 10 * targetIsizeSD));
        }
        
        double endPosterior1 = Math.log(decodePosterior(codedEndPosterior1));
        double endPosterior2 = Math.log(decodePosterior(codedEndPosterior2));
        
        double likelihoodRatio = pDeletion + pISgivenDeletion + endPosterior1 + endPosterior2
                - pNoDeletion - pISgivenNoDeletion;

        return likelihoodRatio;
    }

    private double decodePosterior(double codedPosterior) {
        return codedPosterior == 0 ? 0.0001 : 1 - Math.pow(10.0, codedPosterior / -10.0);
    }
}
