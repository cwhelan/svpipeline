package edu.ohsu.sonmezsysbio.cloudbreak;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 2/28/12
 * Time: 2:54 PM
 */
public abstract class PairedAlignmentScorer {
    public boolean validateInsertSize(int insertSize, String readPairId, Integer maxInsertSize1) {
        if (insertSize == 0) return false;
        if (insertSize > maxInsertSize1) {
            // System.err.println("Pair " + readPairId + ": Insert size would be greater than " + maxInsertSize1 + " - skipping");
            return false;
        }
        return true;
    }

    public abstract double computeDeletionScore(int insertSize, Double targetIsize, Double targetIsizeSD, Double pMappingCorrect);

    public boolean validateMappingOrientations(AlignmentRecord record1, AlignmentRecord record2, boolean matePairs) {
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

    public boolean isMatePairNotSmallFragment(AlignmentRecord record1, AlignmentRecord record2) {
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
