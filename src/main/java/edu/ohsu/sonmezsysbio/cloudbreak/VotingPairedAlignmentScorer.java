package edu.ohsu.sonmezsysbio.cloudbreak;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.exception.MathIllegalArgumentException;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 3/12/12
 * Time: 1:27 PM
 */
public class VotingPairedAlignmentScorer extends PairedAlignmentScorer {
    public double computeDeletionScore(int insertSize, Double targetIsize, Double targetIsizeSD, Double pMappingCorrect) {
        //System.err.println("target isize: " + targetIsize + ", sd " + targetIsizeSD);
        NormalDistribution insertSizeDist = new NormalDistribution(targetIsize, targetIsizeSD);
        // deletion score = codedEndPosterior1 * codedEndPosterior2 * P(X < insertSize - 2 * targetIsizeSD)

        double deletionProb;
        try {
            deletionProb = insertSizeDist.cumulativeProbability(Math.max(0, insertSize - 1.5 * targetIsizeSD));
            //System.err.println("Deletion prob: " + deletionProb);
        } catch (MathIllegalArgumentException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

        double vote = deletionProb > 0.5 ? 1 : -1;


//        double endPosterior1 = decodePosterior(codedEndPosterior1);
//        double endPosterior2 = decodePosterior(codedEndPosterior2);
        //System.err.println("posteriors: " + endPosterior1 + "," + endPosterior2);

        //return deletionProb + endPosterior1 + endPosterior2;
        return vote * Math.exp(pMappingCorrect);
    }
}
