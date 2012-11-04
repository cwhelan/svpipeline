package edu.ohsu.sonmezsysbio.cloudbreak.io;

import edu.ohsu.sonmezsysbio.cloudbreak.AlignmentRecord;
import edu.ohsu.sonmezsysbio.cloudbreak.Cloudbreak;
import edu.ohsu.sonmezsysbio.cloudbreak.ReadPairAlignments;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 8/8/12
 * Time: 2:42 PM
 */
public abstract class BaseAlignmentReader implements AlignmentReader {

    public ReadPairAlignments parsePairAlignmentLine(String line) {
        String[] reads = line.split(Cloudbreak.READ_SEPARATOR);
        String read1AlignmentsString = reads[0];
        String[] read1Alignments = read1AlignmentsString.split(Cloudbreak.ALIGNMENT_SEPARATOR);

        List<AlignmentRecord> read1AlignmentRecords = parseAlignmentsIntoRecords(read1Alignments);

        String read2AlignmentsString = reads[1];
        String[] read2Alignments = read2AlignmentsString.split(Cloudbreak.ALIGNMENT_SEPARATOR);
        List<AlignmentRecord> read2AlignmentRecords = parseAlignmentsIntoRecords(read2Alignments);
        return new ReadPairAlignments(read1AlignmentRecords, read2AlignmentRecords);
    }

    public List<AlignmentRecord> parseAlignmentsIntoRecords(String[] alignments) {
        List<AlignmentRecord> read1AlignmentList = new ArrayList<AlignmentRecord>();
        for (String alignment : alignments) {
            AlignmentRecord record1 = parseRecord(alignment);
            read1AlignmentList.add(record1);
        }
        return read1AlignmentList;
    }
}
