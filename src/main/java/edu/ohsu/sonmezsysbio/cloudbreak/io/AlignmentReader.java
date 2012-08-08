package edu.ohsu.sonmezsysbio.cloudbreak.io;

import edu.ohsu.sonmezsysbio.cloudbreak.AlignmentRecord;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 8/8/12
 * Time: 2:05 PM
 */
public interface AlignmentReader {
    AlignmentRecord parseRecord(String[] fields);

    List<AlignmentRecord> parseAlignmentsIntoRecords(String[] alignments);
}
