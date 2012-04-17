package edu.ohsu.sonmezsysbio.svpipeline.command;

import com.beust.jcommander.Parameter;
import edu.ohsu.sonmezsysbio.svpipeline.SVPipeline;
import edu.ohsu.sonmezsysbio.svpipeline.mapper.SingleEndAlignmentSummaryMapper;
import edu.ohsu.sonmezsysbio.svpipeline.reducer.SingleEndAlignmentSummaryReducer;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;

import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/16/12
 * Time: 1:39 PM
 */
public class CommandSummarizeAlignments implements SVPipelineCommand {
    @Parameter(names = {"--inputHDFSDir"}, required = true)
    String inputHDFSDir;

    @Parameter(names = {"--outputHDFSDir"}, required = true)
    String outputHDFSDir;

    public void run(Configuration conf) throws Exception {
        runHadoopJob(conf);
    }

    private void runHadoopJob(Configuration configuration) throws IOException {
        JobConf conf = new JobConf(configuration);

        conf.setJobName("Summarize Single End Alignments");
        conf.setJarByClass(SVPipeline.class);
        FileInputFormat.addInputPath(conf, new Path(inputHDFSDir));
        Path outputDir = new Path(outputHDFSDir);
        FileSystem.get(conf).delete(outputDir);

        FileOutputFormat.setOutputPath(conf, outputDir);

        conf.setInputFormat(TextInputFormat.class);

        conf.setMapperClass(SingleEndAlignmentSummaryMapper.class);
        conf.setMapOutputKeyClass(Text.class);
        conf.setMapOutputValueClass(Text.class);

        conf.setCombinerClass(SingleEndAlignmentSummaryReducer.class);

        conf.setReducerClass(SingleEndAlignmentSummaryReducer.class);

        conf.setOutputKeyClass(Text.class);
        conf.setOutputValueClass(Text.class);

        conf.setCompressMapOutput(true);

        JobClient.runJob(conf);

    }
}