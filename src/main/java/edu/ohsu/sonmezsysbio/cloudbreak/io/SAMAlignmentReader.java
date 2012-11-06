package edu.ohsu.sonmezsysbio.cloudbreak.io;

import edu.ohsu.sonmezsysbio.cloudbreak.AlignmentRecord;
import edu.ohsu.sonmezsysbio.cloudbreak.ReadPairAlignments;
import edu.ohsu.sonmezsysbio.cloudbreak.SAMRecord;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 10/15/12
 * Time: 1:36 PM
 */
public class SAMAlignmentReader extends BaseAlignmentReader {

    public AlignmentRecord parseRecord(String alignmentRecord) {
        return parseRecord(alignmentRecord.split("\t"));
    }

    public AlignmentRecord parseRecord(String[] fields) {
        return SAMRecord.parseSamRecord(fields);
    }

    public double probabilityMappingIsCorrect(AlignmentRecord record1, AlignmentRecord record2, ReadPairAlignments readPairAlignments) {
        double r1normalization = sumAlignmentScores(readPairAlignments.getRead1Alignments());
        double r2normalization = sumAlignmentScores(readPairAlignments.getRead2Alignments());
        return Math.log(((double) record1.getAlignmentScore()) / r1normalization) + Math.log(((double) record2.getAlignmentScore()) / r2normalization);
    }

    private double sumAlignmentScores(List<AlignmentRecord> alignments) {
        double sumAlignmentScores = 0;
        for (AlignmentRecord record : alignments) {
            SAMRecord samRecord = (SAMRecord) record;
            sumAlignmentScores += samRecord.getAlignmentScore();
        }
        return sumAlignmentScores;
    }

}
