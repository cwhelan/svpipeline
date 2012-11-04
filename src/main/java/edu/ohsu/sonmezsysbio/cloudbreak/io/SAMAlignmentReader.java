package edu.ohsu.sonmezsysbio.cloudbreak.io;

import edu.ohsu.sonmezsysbio.cloudbreak.AlignmentRecord;
import edu.ohsu.sonmezsysbio.cloudbreak.ReadPairAlignments;
import edu.ohsu.sonmezsysbio.cloudbreak.SAMRecord;

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

    public double probabilityMappingIsCorrect(AlignmentRecord record1, AlignmentRecord record2) {
        return 0;
    }

    public void resetForReadPairAlignemnts(ReadPairAlignments readPairAlignments) {
    }
}
