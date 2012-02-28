package edu.ohsu.sonmezsysbio.svpipeline.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.ohsu.sonmezsysbio.svpipeline.mapper.NovoalignMatePairMapper;
import edu.ohsu.sonmezsysbio.svpipeline.SVPipeline;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;

import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/18/11
 * Time: 2:01 PM
 */
@Parameters(separators = "=", commandDescription = "Run a novoalign mate pair alignment")
public class CommandNovoalignMatePair implements SVPipelineCommand {

    @Parameter(names = {"--HDFSDir"}, required = true)
    String hdfsDir;

    @Parameter(names = {"--fastqFile1"}, required = true)
    String readFile1;

    @Parameter(names = {"--fastqFile2"}, required = true)
    String readFile2;

    @Parameter(names = {"--reference"}, required = true)
    String reference;

    @Parameter(names = {"--repeatReport"}, required = false)
    String repeatReport = "Random";

    @Parameter(names = {"--targetIsize"}, required = true)
    int targetIsize;

    @Parameter(names = {"--targetIsizeSD"}, required = true)
    int targetIsizeSD;

    @Parameter(names = {"--libraryName"}, required = true)
    String libraryName;

    private int numPairs;

    public void copyReadFilesToHdfs() throws IOException {
        Configuration config = new Configuration();

        FileSystem hdfs = FileSystem.get(config);

        FSDataOutputStream outputStream = hdfs.create(new Path(hdfsDir + "/" + "novoMatePairIn.txt"));
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(outputStream));

        BufferedReader inputReader1 = new BufferedReader(new FileReader(new File(readFile1)));
        BufferedReader inputReader2 = new BufferedReader(new FileReader(new File(readFile2)));

        numPairs = 0;
        try {
            String convertedFastqLine = readFastqPair(inputReader1, inputReader2);
            while (convertedFastqLine != null) {
                writer.write(convertedFastqLine);
                convertedFastqLine = readFastqPair(inputReader1, inputReader2);
                numPairs++;
            }
        } finally {
            writer.close();
            inputReader1.close();
            inputReader2.close();
        }

    }

    private String readFastqPair(BufferedReader inputReader1, BufferedReader inputReader2) throws IOException {
        String read1 = inputReader1.readLine();
        if (read1 == null) {
            if (inputReader2.readLine() != null) {
                throw new IOException("Read file 2 was longer than read file 1");
            }
            return null;
        }

        String seq1 = inputReader1.readLine();
        String sep1 = inputReader1.readLine();
        String qual1 = inputReader1.readLine();

        String read2 = inputReader2.readLine();
        String seq2 = inputReader2.readLine();
        String sep2 = inputReader2.readLine();
        String qual2 = inputReader2.readLine();

        StringBuffer lineBuffer = new StringBuffer();
        lineBuffer.append(read1).append("\t").append(seq1).append("\t").append(sep1).append("\t").append(qual1);
        lineBuffer.append("\t");
        lineBuffer.append(read2).append("\t").append(seq2).append("\t").append(sep2).append("\t").append(qual2);
        lineBuffer.append("\n");
        return lineBuffer.toString();
    }

    public void runHadoopJob(Configuration configuration) throws IOException {
        JobConf conf = new JobConf(configuration);

        conf.setJobName("Novoalign Mate Pair Alignment");
        conf.setJarByClass(SVPipeline.class);
        FileInputFormat.addInputPath(conf, new Path(hdfsDir));
        Path outputDir = new Path("/user/whelanch/tmp/svout/");
        FileSystem.get(conf).delete(outputDir);

        FileOutputFormat.setOutputPath(conf, outputDir);

        conf.setInputFormat(TextInputFormat.class);

        conf.set("mapred.task.timeout", "3600000");
        conf.set("novoalign.reference", reference);
        conf.set("novoalign.repeatReport", repeatReport);
        conf.set("novoalign.targetIsize", String.valueOf(targetIsize));
        conf.set("novoalign.targetIsizeSD", String.valueOf(targetIsizeSD));
        conf.set("novoalign.libraryName", libraryName);

        conf.setMapperClass(NovoalignMatePairMapper.class);
        conf.setMapOutputKeyClass(Text.class);
        conf.setMapOutputValueClass(Text.class);
        conf.setNumMapTasks(numPairs / 60000 + 1);

        conf.setOutputKeyClass(Text.class);
        conf.setOutputValueClass(DoubleWritable.class);
        conf.setCompressMapOutput(true);

        conf.setNumReduceTasks(0);

        JobClient.runJob(conf);

    }

    public void run(Configuration conf) throws IOException {
        copyReadFilesToHdfs();
        runHadoopJob(conf);

    }
}
