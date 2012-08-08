package edu.ohsu.sonmezsysbio.cloudbreak.io;

import edu.ohsu.sonmezsysbio.cloudbreak.AlignmentRecord;
import edu.ohsu.sonmezsysbio.cloudbreak.NovoalignNativeRecord;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 8/8/12
 * Time: 1:56 PM
 */
public class NovoalignAlignmentReader implements AlignmentReader {
    public AlignmentRecord parseRecord(String[] fields) {
        NovoalignNativeRecord record = new NovoalignNativeRecord();
        record.setReadId(fields[0]);
        record.setSequence(fields[2]);
        record.setMappingStatus(fields[4]);
        if (record.isMapped()) {
            String recordReferenceName = fields[7];
            // cut off the ">" that starts the chromosome name
            if (recordReferenceName.startsWith(">")) {
                recordReferenceName = recordReferenceName.substring(1);
            }
            record.setChromsomeName(recordReferenceName);

            record.setPosition(Integer.parseInt(fields[8]));
            record.setPosteriorProb(Double.parseDouble(fields[6]));
            record.setForward("F".equals(fields[9]));
        }

        return record;

    }

    public List<AlignmentRecord> parseAlignmentsIntoRecords(String[] alignments) {
        List<AlignmentRecord> read1AlignmentList = new ArrayList<AlignmentRecord>();
        for (String alignment : alignments) {
            String[] fields1 = alignment.split("\t");
            AlignmentRecord record1 = parseRecord(fields1);
            read1AlignmentList.add(record1);
        }
        return read1AlignmentList;
    }
}
