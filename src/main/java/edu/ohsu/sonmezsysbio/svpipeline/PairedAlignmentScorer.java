package edu.ohsu.sonmezsysbio.svpipeline;

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
            System.err.println("Pair " + readPairId + ": Insert size would be greater than " + maxInsertSize1 + " - skipping");
            return false;
        }
        return true;
    }

    public abstract double computeDeletionScore(int insertSize, Double targetIsize, Double targetIsizeSD, Double pMappingCorrect);

    public boolean validateMappingOrientations(NovoalignNativeRecord record1, NovoalignNativeRecord record2, boolean matePairs) {
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

    public boolean isMatePairNotSmallFragment(NovoalignNativeRecord record1, NovoalignNativeRecord record2) {
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

    public double probabilityMappingIsCorrect(int codedEndPosterior1, int codedEndPosterior2) {
        double endPosterior1 = Math.log(decodePosterior(codedEndPosterior1));
        double endPosterior2 = Math.log(decodePosterior(codedEndPosterior2));

        return endPosterior1 + endPosterior2;
    }

    double decodePosterior(double codedPosterior) {
        return codedPosterior == 0 ? 0.0001 : 1 - Math.pow(10.0, codedPosterior / -10.0);
    }
}
