package edu.ohsu.sonmezsysbio.cloudbreak.mapper;

import edu.ohsu.sonmezsysbio.cloudbreak.NovoalignNativeRecord;
import edu.ohsu.sonmezsysbio.cloudbreak.SVPipeline;
import org.apache.hadoop.io.LongWritable;
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
public class SingleEndAlignmentSummaryMapper extends MapReduceBase implements Mapper<LongWritable, Text, Text, Text> {

    Text outKey = new Text("k");

    public void map(LongWritable key, Text value, OutputCollector output, Reporter reporter) throws IOException {
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

        output.collect(outKey, new Text("1\t" + (read1AlignmentRecords.size() * read2AlignmentRecords.size())));

    }

}