package edu.ohsu.sonmezsysbio.svpipeline;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 2/28/12
 * Time: 2:54 PM
 */
public class SingleEndAlignmentScorer {
    public static boolean validateInsertSize(int insertSize, String readPairId, Integer maxInsertSize1) {
        boolean validInsertSize = true;
        if (insertSize > maxInsertSize1) {
            System.err.println("Pair " + readPairId + ": Insert size would be greater than " + maxInsertSize1 + " - skipping");
            validInsertSize = false;
        }
        if (! validInsertSize) return true;
        return false;
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

    public static boolean validateMappingOrientations(NovoalignNativeRecord record1, NovoalignNativeRecord record2, boolean matePairs) {
        if (matePairs) {
            if (record1.isForward() && ! record2.isForward()) {
            } else if (!record1.isForward() && record2.isForward()) {
            } else {
                return false;
            }
        } else {
            if (record1.isForward() && ! record2.isForward()) {
                if (record1.getPosition() > record2.getPosition()) return false;
            } else if (!record1.isForward() && record2.isForward()) {
                if (record1.getPosition() < record2.getPosition()) return false;
            } else {
                return false;
            }
        }
        return true;
    }

    public static boolean isMatePairNotSmallFragment(NovoalignNativeRecord record1, NovoalignNativeRecord record2) {
        boolean matePair = false;
        if (record1.isForward() && ! record2.isForward()) {
            if (record1.getPosition() - record2.getPosition() > 0) matePair = true;
            if (record1.getPosition() - record2.getPosition() < 0 &&
                    record1.getPosition() - record2.getPosition() > -500) matePair = false;
        } else if (!record1.isForward() && record2.isForward()) {
            if (record1.getPosition() - record2.getPosition() < 0) matePair = true;
            if (record1.getPosition() - record2.getPosition() > 0 &&
                    record1.getPosition() - record2.getPosition() < 500) matePair = false;
        }
        return matePair;
    }
}
