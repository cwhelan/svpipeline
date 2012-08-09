package edu.ohsu.sonmezsysbio.cloudbreak.io;

import edu.ohsu.sonmezsysbio.cloudbreak.AlignmentRecord;
import edu.ohsu.sonmezsysbio.cloudbreak.MrfastAlignmentRecord;
import edu.ohsu.sonmezsysbio.cloudbreak.ReadPairAlignments;

import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 8/8/12
 * Time: 2:43 PM
 */
public class MrfastAlignmentReader extends BaseAlignmentReader {
    public AlignmentRecord parseRecord(String[] fields) {
        MrfastAlignmentRecord record = new MrfastAlignmentRecord();

        // return readId + "\t" + orientation + "\t" + chrom + "\t" + position + "\t" + nm + "\t" + sequenceLength;
        record.setReadId(fields[0]);
        record.setForward("F".equals(fields[1]));
        record.setChromosomeName(fields[2]);
        record.setPosition(Integer.parseInt(fields[3]));
        record.setMismatches(Integer.parseInt(fields[4]));
        record.setSequenceLength(Integer.parseInt(fields[5]));
        return record;

    }

    public double probabilityMappingIsCorrect(AlignmentRecord record1, AlignmentRecord record2, ReadPairAlignments readPairAlignments) {
        Map<AlignmentRecord, Double> read1AlignmentScores = computeAlignmentScores(readPairAlignments.getRead1Alignments());
        return 0;
    }

    static Map<AlignmentRecord, Double> computeAlignmentScores(List<AlignmentRecord> read1Alignments) {
        return null;
    }

}
