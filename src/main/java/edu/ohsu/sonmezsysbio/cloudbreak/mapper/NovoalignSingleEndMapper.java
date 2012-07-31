package edu.ohsu.sonmezsysbio.cloudbreak.mapper;

import edu.ohsu.sonmezsysbio.cloudbreak.NovoalignNativeRecord;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;

import java.io.*;
import java.util.Arrays;
import java.util.zip.GZIPOutputStream;

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
    private String threshold;
    private String baseQualityFormat;
    private String novoalignExecutable;

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
            s1File = new File(localDir + "/temp1_sequence.fastq.gz").getAbsoluteFile();
            s1File.createNewFile();
            s1FileWriter = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(s1File)));

            reference = job.get("novoalign.reference");
            threshold = job.get("novoalign.threshold");
            baseQualityFormat = job.get("novoalign.quality.format");
            novoalignExecutable = job.get("novoalign.executable");

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
        if (! s1File.exists()) {
            System.err.println("file does not exist: " + s1File.getPath());
        } else {
            System.err.println("read file length: " + s1File.length());
        }

        File indexFile = new File(reference);
        if (! indexFile.exists()) {
            System.err.println("index file does not exist: " + indexFile.getPath());
        } else {
            System.err.println("index file length: " + indexFile.length());
        }

        String[] commandLine = buildCommandLine(novoalignExecutable, reference, s1File.getPath(), threshold, baseQualityFormat);
        System.err.println("Executing command: " + Arrays.toString(commandLine));
        Process p = Runtime.getRuntime().exec(commandLine);
        System.err.println("Exec'd");

        BufferedReader stdInput = new BufferedReader(new
                         InputStreamReader(p.getInputStream()));

        readAlignments(stdInput, p.getErrorStream());
        done = true;

    }

    protected void readAlignments(BufferedReader stdInput, InputStream errorStream) throws IOException {
        String outLine;
        while ((outLine = stdInput.readLine()) != null) {
            // System.err.println("LINE: " + outLine);
            if (outLine.startsWith("#"))  {
                System.err.println("COMMENT LINE: " + outLine);
                continue;
            }
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

        String errLine;
        BufferedReader errorReader = new BufferedReader(new InputStreamReader(errorStream));
        while ((errLine = errorReader.readLine()) != null) {
            System.err.println("ERROR: " + errLine);
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

    protected static String[] buildCommandLine(String novoalignExecutable, String reference, String path1, String threshold, String baseQualityFormat) {
        String[] commandArray = {
                "./" + novoalignExecutable,
                "-d", reference,
                "-c", "1",
                "-f", path1,
                "-F", baseQualityFormat,
                "-k", "-K", "calfile.txt", "-q", "5",
                "-r", "Ex", "10", "-t", threshold, "-x", "10"
        };
        return commandArray;
    }
}
