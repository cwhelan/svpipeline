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
        NormalDistribution insertSizeDist = new NormalDistribution(targetIsize, targetIsizeSD);

        double deletionProb;
        try {
            deletionProb = insertSizeDist.cumulativeProbability(Math.max(0, insertSize - 1.5 * targetIsizeSD));
            //System.err.println("Deletion prob: " + deletionProb);
        } catch (MathIllegalArgumentException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

        double vote = deletionProb > 0.5 ? 1 : -1;

        return vote * Math.exp(pMappingCorrect);
    }
}
