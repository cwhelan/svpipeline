package edu.ohsu.sonmezsysbio.svpipeline.mapper;

import edu.ohsu.sonmezsysbio.svpipeline.NovoalignNativeRecord;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;

import javax.sound.sampled.Line;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/21/11
 * Time: 5:36 PM
 */
public class NovoalignSingleEndMapper extends MapReduceBase implements Mapper<LongWritable, Text, Text, Text> {

    private OutputCollector<Text, Text> output;
    private String localDir;
    private Writer s1FileWriter;
    private File s1File;
    private String reference;
    private Reporter reporter;
    private static boolean done = false;

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

    @Override
    public void configure(JobConf job) {
        super.configure(job);

        System.err.println("Current dir: " + new File(".").getAbsolutePath());

        this.localDir = job.get("mapred.child.tmp");
        try {
            s1File = new File(localDir + "/temp1_sequence.txt").getAbsoluteFile();
            s1File.createNewFile();
            s1FileWriter = new BufferedWriter(new FileWriter(s1File));

            reference = job.get("novoalign.reference");

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void map(LongWritable key, Text value, OutputCollector<Text, Text> output, Reporter reporter) throws IOException {
        if (this.output == null) {
            this.output = output;
        }
        if (this.reporter == null) {
            this.reporter = reporter;
        }

        String line = value.toString();
        String[] fields = line.split("\t");

        s1FileWriter.write(fields[0] + "\n");
        s1FileWriter.write(fields[1] + "\n");
        s1FileWriter.write(fields[2] + "\n");
        s1FileWriter.write(fields[3] + "\n");

        reporter.progress();
        //System.out.println("Done with map method, real work will happen in close");
    }

    class ProgressReporter implements Runnable {
        public void run() {
            while (! done) {
                try {
                    Thread.sleep(10000l);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
                reporter.progress();
            }

        }
    }

    @Override
    public void close() throws IOException {
        super.close();

        s1FileWriter.close();

        //Thread progressThread = new Thread(new ProgressReporter());
        //progressThread.start();
        String[] commandLine = buildCommandLine(reference, s1File.getPath());
        System.err.println("Executing command: " + Arrays.toString(commandLine));
        Process p = Runtime.getRuntime().exec(commandLine);


        BufferedReader stdInput = new BufferedReader(new
                         InputStreamReader(p.getInputStream()));

        readAlignments(stdInput, p.getErrorStream());
        done = true;

    }

    protected void readAlignments(BufferedReader stdInput, InputStream errorStream) throws IOException {
        String outLine;
        while ((outLine = stdInput.readLine()) != null) {
            //System.err.println("LINE: " + outLine);
            if (outLine.startsWith("#"))  continue;
            if (outLine.startsWith("Novoalign") || outLine.startsWith("Exception")) {
                String error = printErrorStream(errorStream);
                throw new RuntimeException(error);
            }

            String readPairId = outLine.substring(0,outLine.indexOf('\t')-2);
            NovoalignNativeRecord alignment = NovoalignNativeRecord.parseRecord(outLine.split("\t"));

            if (! alignment.isMapped()) {
                continue;
            }

            output.collect(new Text(readPairId), new Text(outLine));

        }
    }

    private String printErrorStream(InputStream errorStream) throws IOException {
        String outLine;BufferedReader stdErr = new BufferedReader(new
                InputStreamReader(errorStream));
        String firstErrorLine = null;
        while ((outLine = stdErr.readLine()) != null) {
            if (firstErrorLine == null) firstErrorLine = outLine;
            System.err.println(outLine);
        }
        return firstErrorLine;
    }

    protected static String[] buildCommandLine(String reference, String path1) {
        String[] commandArray = {
                "/g/whelanch/software/bin/" + "novoalign",
                "-d", reference,
                "-c", "1",
                "-f", path1,
                "-F", "ILMFQ",
                "-k", "-K", "calfile.txt",
                "-r", "Ex", "10", "-t", "250"
                //"-a", "GATCGGAAGAGCGGTTCAGCA", "GATCGGAAGAGCGTCGTGTAGGGA",
        };
//        String args = String.format("-d %s -c 1 -f %s %s -F ILMFQ -k -K calfile.txt -i MP %s,%s 150,50" +
//                " -a GATCGGAAGAGCGGTTCAGCA GATCGGAAGAGCGTCGTGTAGGGA " +
//                " -r %s -oSAM \"@RG\tID:RGID\tPU:ILLUMINA\tLB:%s\tSM:%s\" ",
//                reference, path1, path2, targetIsize, targetIsizeSD, repeatReport, libraryName, libraryName); // todo: read group ID?
//
//        return "/g/whelanch/software/bin/" + "novoalign " + args; //todo: unhardcode path
        return commandArray;
    }
}
