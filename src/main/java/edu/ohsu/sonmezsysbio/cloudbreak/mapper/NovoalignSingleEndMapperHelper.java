package edu.ohsu.sonmezsysbio.cloudbreak.mapper;

import edu.ohsu.sonmezsysbio.cloudbreak.NovoalignNativeRecord;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 2/28/12
 * Time: 3:25 PM
 */
public class NovoalignSingleEndMapperHelper {
    public static List<NovoalignNativeRecord> parseAlignmentsIntoRecords(String[] alignments) {
        List<NovoalignNativeRecord> read1AlignmentList = new ArrayList<NovoalignNativeRecord>();
        for (String alignment : alignments) {
            String[] fields1 = alignment.split("\t");
            NovoalignNativeRecord record1 = NovoalignNativeRecord.parseRecord(fields1);
            read1AlignmentList.add(record1);
        }
        return read1AlignmentList;
    }
}
