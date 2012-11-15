package edu.ohsu.sonmezsysbio.cloudbreak.mapper;

import edu.ohsu.sonmezsysbio.cloudbreak.AlignmentRecord;
import edu.ohsu.sonmezsysbio.cloudbreak.io.SAMAlignmentReader;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.JobConf;
import org.apache.hadoop.mapred.OutputCollector;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/21/11
 * Time: 5:36 PM
 */
public class RazerS3SingleEndMapper extends SingleEndAlignmentMapper {

    private static Logger logger = Logger.getLogger(RazerS3SingleEndMapper.class);

    //{ logger.setLevel(Level.DEBUG); }

    private OutputCollector<Text, Text> output;
    private String localDir;
    private String reference;
    private String numReports;
    private String razerS3Executable;

    @Override
    protected boolean getCompressTempReadFile() {
        return false;
    }

    @Override
    public void configure(JobConf job) {
        super.configure(job);

        logger.debug("Current dir: " + new File(".").getAbsolutePath());

        this.localDir = job.get("mapred.child.tmp");
        reference = job.get("razers3.reference");
        numReports = job.get("razers3.num.reports");
        razerS3Executable = job.get("razers3.executable");

    }

    @Override
    public void close() throws IOException {
        super.close();

        s1FileWriter.close();

        if (! s1File.exists()) {
            logger.error("file does not exist: " + s1File.getPath());
        } else {
            logger.info("read file length: " + s1File.length());
        }

        String referenceBaseName = new File(reference).getName();
        String[] commandLine = buildCommandLine(razerS3Executable, referenceBaseName, s1File.getPath(), numReports);
        logger.debug("Executing command: " + Arrays.toString(commandLine));
        Process p = Runtime.getRuntime().exec(commandLine);
        logger.debug("Exec'd");

        try {
            p.waitFor();
        } catch (InterruptedException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
        logger.debug("done");

        BufferedReader stdInput = new BufferedReader(new FileReader("map.result"));
        readAlignments(stdInput, p.getErrorStream());
    }

    protected void readAlignments(BufferedReader stdInput, InputStream errorStream) throws IOException {
        String outLine;
        SAMAlignmentReader alignmentReader = new SAMAlignmentReader();
        while ((outLine = stdInput.readLine()) != null) {
            if (logger.isDebugEnabled()) {
                logger.debug("LINE: " + outLine);
            }
            if (outLine.startsWith("@"))  {
                logger.debug("SAM HEADER LINE: " + outLine);
                continue;
            }

            String readPairId = outLine.substring(0,outLine.indexOf('\t')-2);
            AlignmentRecord alignment = alignmentReader.parseRecord(outLine);

            if (! alignment.isMapped()) {
                continue;
            }

            getOutput().collect(new Text(readPairId), new Text(outLine));

        }

        String errLine;
        BufferedReader errorReader = new BufferedReader(new InputStreamReader(errorStream));
        while ((errLine = errorReader.readLine()) != null) {
            logger.error("ERROR: " + errLine);
        }
    }

    private String printErrorStream(InputStream errorStream) throws IOException {
        String outLine;BufferedReader stdErr = new BufferedReader(new
                InputStreamReader(errorStream));
        String firstErrorLine = null;
        while ((outLine = stdErr.readLine()) != null) {
            if (firstErrorLine == null) firstErrorLine = outLine;
            logger.error(outLine);
        }
        return firstErrorLine;
    }

    protected static String[] buildCommandLine(String razers3Executable, String referenceBaseName, String path1, String numReports) {
        String[] commandArray = {
                "./" + razers3Executable,
                "-o", "map.result", "-of", "sam", "-rr", "100", "-i", "96", "-m", numReports,
                referenceBaseName,
                path1
        };
        return commandArray;
    }
}
