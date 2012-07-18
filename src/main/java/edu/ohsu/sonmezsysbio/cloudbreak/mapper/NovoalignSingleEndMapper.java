package edu.ohsu.sonmezsysbio.cloudbreak.mapper;

import edu.ohsu.sonmezsysbio.cloudbreak.NovoalignNativeRecord;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.JobConf;

import java.io.*;
import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/21/11
 * Time: 5:36 PM
 */
public class NovoalignSingleEndMapper extends SingleEndAlignmentMapper {

    private String reference;
    private String threshold;
    private String baseQualityFormat;

    @Override
    public void configure(JobConf job) {
        super.configure(job);

        System.err.println("Current dir: " + new File(".").getAbsolutePath());

        reference = job.get("novoalign.reference");
        threshold = job.get("novoalign.threshold");
        baseQualityFormat = job.get("novoalign.quality.format");

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

        String[] commandLine = buildCommandLine(reference, s1File.getPath(), threshold, baseQualityFormat);
        System.err.println("Executing command: " + Arrays.toString(commandLine));
        Process p = Runtime.getRuntime().exec(commandLine);
        System.err.println("Exec'd");

        BufferedReader stdInput = new BufferedReader(new
                         InputStreamReader(p.getInputStream()));

        readAlignments(stdInput, p.getErrorStream());
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

            getOutput().collect(new Text(readPairId), new Text(outLine));

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

    protected static String[] buildCommandLine(String reference, String path1, String threshold, String baseQualityFormat) {
        String[] commandArray = {
                "/g/whelanch/software/bin/" + "novoalign",
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
