package edu.ohsu.sonmezsysbio.svpipeline.reducer;

import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;

import java.io.IOException;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/16/12
 * Time: 1:46 PM
 */
public class SingleEndAlignmentSummaryReducer extends MapReduceBase implements Reducer<Text, Text, Text, Text> {
    public void reduce(Text key, Iterator<Text> values, OutputCollector<Text, Text> output, Reporter reporter) throws IOException {
        Long[] totals = new Long[2];
        for (Text val = values.next(); values.hasNext();) {
            String[] fields = val.toString().split("\t");
            for (int i = 0; i < fields.length; i++) {
                totals[i] += Long.parseLong(fields[i]);
            }
        }
        output.collect(key, new Text(totals[0] + "\t" + totals[1]));
    }

}
