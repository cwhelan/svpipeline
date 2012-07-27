package edu.ohsu.sonmezsysbio.cloudbreak.command;

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
public class CommandReadPairedEndFilesIntoHDFS implements CloudbreakCommand {

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
            readFile(writer, readFile1, readFile2);
        } finally {
            writer.close();
        }

    }

    private void readFile(BufferedWriter writer, String pathname1, String pathname2) throws IOException {
        BufferedReader inputReader1;
        BufferedReader inputReader2 = null;

        inputReader1 = openFile(pathname1);
        inputReader2 = openFile(pathname2);

        numRecords = 0;
        try {
            String convertedFastqLine = readFastqEntries(inputReader1, inputReader2);
            while (convertedFastqLine != null) {
                writer.write(convertedFastqLine);
                convertedFastqLine = readFastqEntries(inputReader1, inputReader2);
                numRecords++;
            }
        } finally {
            inputReader1.close();
        }
    }

    private BufferedReader openFile(String pathname) throws IOException {
        BufferedReader inputReader1;
        if (pathname.endsWith("gz")) {
            inputReader1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(pathname))));
        } else {
            inputReader1 = new BufferedReader(new FileReader(new File(pathname)));
        }
        return inputReader1;
    }

    private String readFastqEntries(BufferedReader inputReader1, BufferedReader inputReader2) throws IOException {
        String read1 = inputReader1.readLine();
        if (read1 == null) {
            return null;
        }

        String seq1 = inputReader1.readLine();
        String sep1 = inputReader1.readLine();
        String qual1 = inputReader1.readLine();

        String read2 = inputReader2.readLine();
        if (read2 == null) {
            return null;
        }

        String seq2 = inputReader2.readLine();
        String sep2 = inputReader2.readLine();
        String qual2 = inputReader2.readLine();

        String readPrefix = greatestCommonPrefix(read1, read2);

        StringBuffer lineBuffer = new StringBuffer();

        lineBuffer.append(readPrefix);
        lineBuffer.append("/1");
        lineBuffer.append("\t").append(seq1).append("\t").append(sep1).append("\t").append(qual1);
        lineBuffer.append("\n");

        lineBuffer.append(readPrefix);
        lineBuffer.append("/2");
        lineBuffer.append("\t").append(seq2).append("\t").append(sep2).append("\t").append(qual2);
        lineBuffer.append("\n");

        return lineBuffer.toString();
    }

    protected static String greatestCommonPrefix(String read1, String read2) {
        int i = 0;
        while (i < Math.max(read1.length(), read2.length()) && read1.charAt(i) == read2.charAt(i)) {
            i = i + 1;
        }
        return read1.substring(0,i);
    }

    public void run(Configuration conf) throws Exception {
        copyReadFilesToHdfs();
        System.out.println("Loaded " + numRecords + " records.");
    }
}
