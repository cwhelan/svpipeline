package edu.ohsu.sonmezsysbio.svpipeline.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;

import java.io.*;
import java.util.zip.GZIPInputStream;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 2/4/12
 * Time: 3:38 PM
 */
@Parameters(separators = "=", commandDescription = "Load paired fastq files into HDFS")
public class CommandReadPairedEndFilesIntoHDFS implements SVPipelineCommand {

    @Parameter(names = {"--HDFSDataDir"}, required = true)
    String hdfsDataDir;

    @Parameter(names = {"--fastqFile1"}, required = true)
    String readFile1;

    @Parameter(names = {"--fastqFile2"}, required = true)
    String readFile2;

    private int numRecords;

    public void copyReadFilesToHdfs() throws IOException {
        Configuration config = new Configuration();

        FileSystem hdfs = FileSystem.get(config);

        FSDataOutputStream outputStream = hdfs.create(new Path(hdfsDataDir + "/" + "novoIn.txt"));
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(outputStream));

        try {
            readFile(writer, readFile1);
            readFile(writer, readFile2);
        } finally {
            writer.close();
        }

    }

    private void readFile(BufferedWriter writer, String pathname) throws IOException {
        BufferedReader inputReader1 = null;

        if (pathname.endsWith("gz")) {
            inputReader1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(pathname))));
        } else {
            inputReader1 = new BufferedReader(new FileReader(new File(pathname)));
        }

        numRecords = 0;
        try {
            String convertedFastqLine = readFastqEntry(inputReader1);
            while (convertedFastqLine != null) {
                writer.write(convertedFastqLine);
                convertedFastqLine = readFastqEntry(inputReader1);
                numRecords++;
            }
        } finally {
            inputReader1.close();
        }
    }

    private String readFastqEntry(BufferedReader inputReader1) throws IOException {
        String read1 = inputReader1.readLine();
        if (read1 == null) {
            return null;
        }

        String seq1 = inputReader1.readLine();
        String sep1 = inputReader1.readLine();
        String qual1 = inputReader1.readLine();

        StringBuffer lineBuffer = new StringBuffer();
        lineBuffer.append(read1).append("\t").append(seq1).append("\t").append(sep1).append("\t").append(qual1);
        lineBuffer.append("\n");
        return lineBuffer.toString();
    }

    public void run(Configuration conf) throws Exception {
        copyReadFilesToHdfs();
        System.out.println("Loaded " + numRecords + " records.");
    }
}
