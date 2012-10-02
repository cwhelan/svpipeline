package edu.ohsu.sonmezsysbio.cloudbreak.mapper;

import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;

import java.io.*;
import java.util.zip.GZIPOutputStream;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 7/17/12
 * Time: 5:27 PM
 */
public abstract class SingleEndAlignmentMapper extends MapReduceBase implements Mapper<IntWritable, Text, Text, Text> {
    private String localDir;
    protected Writer s1FileWriter;
    protected File s1File;
    protected OutputCollector<Text, Text> output;
    protected Reporter reporter;

    @Override
    public void configure(JobConf job) {
        super.configure(job);
        this.localDir = job.get("mapred.child.tmp");
        try {
            s1File = new File(localDir + "/temp1_sequence.fastq.gz").getAbsoluteFile();
            s1File.createNewFile();
            s1FileWriter = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(s1File)));


        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void map(IntWritable key, Text value, OutputCollector<Text, Text> output, Reporter reporter) throws IOException {
        if (this.output == null) {
            this.output = output;
        }
        if (this.reporter == null) {
            this.reporter = reporter;
        }

        String line = value.toString();
        String[] fields = line.split("\t");

        if (fields[1].length() != fields[3].length()) {
            System.err.println("Warning; mismatching seq and qual lengths in record " + key.toString() + "!");
            System.err.println("Seq:");
            System.err.println(fields[1]);
            System.err.println("Qual:");
            System.err.println(fields[3]);
            System.err.println("DONE WARNING");
        }
        s1FileWriter.write(fields[0] + "\n");
        s1FileWriter.write(fields[1] + "\n");
        s1FileWriter.write(fields[2] + "\n");
        s1FileWriter.write(fields[3] + "\n");

        reporter.progress();
        //System.out.println("Done with map method, real work will happen in close");
    }

    public OutputCollector<Text, Text> getOutput() {
        return output;
    }

    public void setOutput(OutputCollector<Text, Text> output) {
        this.output = output;
    }

    public Reporter getReporter() {
        return reporter;
    }

    public void setReporter(Reporter reporter) {
        this.reporter = reporter;
    }
}
