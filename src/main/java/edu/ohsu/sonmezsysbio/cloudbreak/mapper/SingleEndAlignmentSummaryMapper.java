package edu.ohsu.sonmezsysbio.cloudbreak.mapper;

import edu.ohsu.sonmezsysbio.cloudbreak.CloudbreakMapReduceBase;
import edu.ohsu.sonmezsysbio.cloudbreak.ReadPairAlignments;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.Mapper;
import org.apache.hadoop.mapred.OutputCollector;
import org.apache.hadoop.mapred.Reporter;

import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/16/12
 * Time: 1:44 PM
 */
public class SingleEndAlignmentSummaryMapper extends CloudbreakMapReduceBase implements Mapper<Text, Text, Text, Text> {

    Text outKey = new Text("k");

    public void map(Text key, Text value, OutputCollector output, Reporter reporter) throws IOException {
        String line = value.toString();
        ReadPairAlignments readPairAlignments = parsePairAlignmentLine(line);

        output.collect(outKey, new Text("1\t" + (readPairAlignments.getRead1Alignments().size() * readPairAlignments.getRead2Alignments().size())));

    }

}