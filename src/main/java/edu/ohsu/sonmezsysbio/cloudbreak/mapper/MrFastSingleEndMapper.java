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
    private String mrfastExecutable;

    @Override
    public void configure(JobConf job) {
        super.configure(job);
        System.err.println("Current dir: " + new File(".").getAbsolutePath());
        reference = job.get("mrfast.reference");
        mrfastExecutable = job.get("mrfast.executable");
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

        String[] commandLine = buildCommandLine(mrfastExecutable, reference, s1File.getPath());
        System.err.println("Executing command: " + Arrays.toString(commandLine));
        Process p = Runtime.getRuntime().exec(commandLine);
        System.err.println("Exec'd");

        try {
            p.waitFor();
        } catch (InterruptedException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
        System.err.println("process finished with exit code: " + p.exitValue());

        BufferedReader stdInput = new BufferedReader(new
                InputStreamReader(new GZIPInputStream(new FileInputStream(new File("output.gz")))));

        readAlignments(stdInput, p.getErrorStream());
    }

    protected void readAlignments(BufferedReader stdInput, InputStream errorStream) throws IOException {

        String errLine;
        BufferedReader errorReader = new BufferedReader(new InputStreamReader(errorStream));
        while ((errLine = errorReader.readLine()) != null) {
            System.err.println("ERROR: " + errLine);
        }

        String outLine;
        while ((outLine = stdInput.readLine()) != null) {
            String readPairId = outLine.substring(0,outLine.indexOf('\t')-2);
            String condensedAlignmentLine = condenseAlignmentLine(outLine);
            getOutput().collect(new Text(readPairId), new Text(condensedAlignmentLine));
        }

    }

    /**
     * MRfast SAM output is very long so only pick out the important fields: alignment location and number
     * of mismatches
     * @param outLine
     * @return
     */
    private String condenseAlignmentLine(String outLine) {
        String[] fields = outLine.split("\t");
        String readId = fields[0];
        String orientation = "0".equals(fields[1]) ? "F" : "R";
        String chrom = fields[2];
        String position = fields[3];
        String sequenceLength = String.valueOf(fields[9].length());
        String nm = "NA";
        for (int i = 4; i < fields.length; i++) {
            if (fields[i].startsWith("NM:i:")) {
                nm = fields[i].substring(5);
                break;
            }
        }
        return readId + "\t" + orientation + "\t" + chrom + "\t" + position + "\t" + nm + "\t" + sequenceLength;
    }

    protected static String[] buildCommandLine(String mrfastExecutable, String reference, String path1) {
        String[] commandArray = {
                "./" + mrfastExecutable,
                "--search", reference,
                "--seq", path1,
                "--outcomp",
                "--seqcomp"
        };
        return commandArray;
    }

}
