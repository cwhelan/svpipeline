package edu.ohsu.sonmezsysbio.cloudbreak.io;

import edu.ohsu.sonmezsysbio.cloudbreak.AlignmentRecord;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 8/8/12
 * Time: 2:42 PM
 */
public abstract class BaseAlignmentReader implements AlignmentReader {
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
