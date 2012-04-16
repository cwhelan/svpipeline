package edu.ohsu.sonmezsysbio.svpipeline.mapper;

import edu.ohsu.sonmezsysbio.svpipeline.NovoalignNativeRecord;
import edu.ohsu.sonmezsysbio.svpipeline.SVPipeline;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;

import java.io.IOException;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/16/12
 * Time: 1:44 PM
 */
public class SingleEndAlignmentSummaryMapper extends MapReduceBase implements Mapper<Text, Text, Text, Text> {
    public void map(Text key, Text value, OutputCollector output, Reporter reporter) throws IOException {
        String line = value.toString();
        int firstTabIndex = line.indexOf('\t');
        String lineValues = line.substring(firstTabIndex + 1);

        String[] readAligments = lineValues.split(SVPipeline.READ_SEPARATOR);
        String read1AlignmentsString = readAligments[0];
        String[] read1Alignments = read1AlignmentsString.split(SVPipeline.ALIGNMENT_SEPARATOR);
        List<NovoalignNativeRecord> read1AlignmentRecords = NovoalignSingleEndMapperHelper.parseAlignmentsIntoRecords(read1Alignments);

        String read2AlignmentsString = readAligments[1];
        String[] read2Alignments = read2AlignmentsString.split(SVPipeline.ALIGNMENT_SEPARATOR);
        List<NovoalignNativeRecord> read2AlignmentRecords = NovoalignSingleEndMapperHelper.parseAlignmentsIntoRecords(read2Alignments);

        emitReadPairInfoForAllPairs(key, read1AlignmentRecords, read2AlignmentRecords, output);

    }
    private void emitReadPairInfoForAllPairs(Text key, List<NovoalignNativeRecord> read1AlignmentRecords, List<NovoalignNativeRecord> read2AlignmentRecords, OutputCollector<Text, Text> output) throws IOException {
        int pairs = 0;
        for (NovoalignNativeRecord record1 : read1AlignmentRecords) {
            for (NovoalignNativeRecord record2 : read2AlignmentRecords) {
                pairs += 1;
            }
        }
        output.collect(key, new Text("1\t" + pairs));
    }
}