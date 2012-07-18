package edu.ohsu.sonmezsysbio.cloudbreak.mapper;

import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.JobConf;

import java.io.*;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 7/17/12
 * Time: 5:46 PM
 */
public class MrFastSingleEndMapper extends SingleEndAlignmentMapper {

    private String reference;

    @Override
    public void configure(JobConf job) {
        super.configure(job);
        System.err.println("Current dir: " + new File(".").getAbsolutePath());
        reference = job.get("mrfast.reference");
    }

    @Override
    public void close() throws IOException {
        super.close();

        s1FileWriter.close();

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

        String[] commandLine = buildCommandLine(reference, s1File.getPath());
        System.err.println("Executing command: " + Arrays.toString(commandLine));
        Process p = Runtime.getRuntime().exec(commandLine);
        System.err.println("Exec'd");

        try {
            p.waitFor();
        } catch (InterruptedException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

        BufferedReader stdInput = new BufferedReader(new
                InputStreamReader(new GZIPInputStream(new FileInputStream(new File("output")))));

        readAlignments(stdInput, p.getErrorStream());
    }

    protected void readAlignments(BufferedReader stdInput, InputStream errorStream) throws IOException {
        String outLine;
        while ((outLine = stdInput.readLine()) != null) {
            String readPairId = outLine.substring(0,outLine.indexOf('\t')-2);
            getOutput().collect(new Text(readPairId), new Text(outLine));
        }
        String errLine;
        BufferedReader errorReader = new BufferedReader(new InputStreamReader(errorStream));
        while ((errLine = errorReader.readLine()) != null) {
            System.err.println("ERROR: " + errLine);
        }

    }

    protected static String[] buildCommandLine(String reference, String path1) {
        String[] commandArray = {
                "/g/whelanch/software/bin/" + "mrfast",
                "--search", reference,
                "--seq", path1,
                "--outcomp"
        };
        return commandArray;
    }

}
