package edu.ohsu.sonmezsysbio.svpipeline.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.ohsu.sonmezsysbio.svpipeline.mapper.MappedPairToDeletionScoreMapper;
import edu.ohsu.sonmezsysbio.svpipeline.SVPipeline;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;
import org.apache.hadoop.mapred.lib.IdentityReducer;
import org.apache.hadoop.mapred.lib.KeyFieldBasedComparator;

import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/23/11
 * Time: 10:02 AM
 */
@Parameters(separators = "=", commandDescription = "Calculate Deletion Scores Across the Genome")
public class CommandPileupDeletionScores implements SVPipelineCommand {

    @Parameter(names = {"--inputHDFSDir"}, required = true)
    String inputHDFSDir;

    @Parameter(names = {"--outputHDFSDir"}, required = true)
    String ouptutHDFSDir;

    @Parameter(names = {"--targetIsize"}, required = true)
    int targetIsize;

    @Parameter(names = {"--targetIsizeSD"}, required = true)
    int targetIsizeSD;

    public void run(Configuration conf) throws IOException {
        runHadoopJob(conf);
    }

    private void runHadoopJob(Configuration configuration) throws IOException {
        JobConf conf = new JobConf(configuration);

        conf.setJobName("Pileup Deletion Score");
        conf.setJarByClass(SVPipeline.class);
        FileInputFormat.addInputPath(conf, new Path(inputHDFSDir));
        Path outputDir = new Path(ouptutHDFSDir);
        FileSystem.get(conf).delete(outputDir);

        FileOutputFormat.setOutputPath(conf, outputDir);

        conf.setInputFormat(TextInputFormat.class);

        conf.set("pileupDeletionScore.targetIsize", String.valueOf(targetIsize));
        conf.set("pileupDeletionScore.targetIsizeSD", String.valueOf(targetIsizeSD));

        conf.setMapperClass(MappedPairToDeletionScoreMapper.class);
        conf.setMapOutputKeyClass(Text.class);
        conf.setMapOutputValueClass(DoubleWritable.class);

        //conf.setReducerClass(DeletionScorePileupReducer.class);
        conf.setReducerClass(IdentityReducer.class);

        conf.setOutputKeyClass(Text.class);
        conf.setOutputValueClass(DoubleWritable.class);
        conf.setKeyFieldComparatorOptions("-k 1,1 -k 2,2n");
        conf.setOutputKeyComparatorClass(KeyFieldBasedComparator.class);

        conf.setCompressMapOutput(true);

        JobClient.runJob(conf);

    }
}
