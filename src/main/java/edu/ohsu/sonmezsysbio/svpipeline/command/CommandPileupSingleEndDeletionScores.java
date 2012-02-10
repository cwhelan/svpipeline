package edu.ohsu.sonmezsysbio.svpipeline.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.ohsu.sonmezsysbio.svpipeline.SVPipeline;
import edu.ohsu.sonmezsysbio.svpipeline.mapper.SingleEndAlignmentsToDeletionScoreMapper;
import edu.ohsu.sonmezsysbio.svpipeline.reducer.SingleEndDeletionScorePileupReducer;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;

import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/23/11
 * Time: 10:02 AM
 */
@Parameters(separators = "=", commandDescription = "Calculate Deletion Scores Across the Genome")
public class CommandPileupSingleEndDeletionScores implements SVPipelineCommand {

    // todo: rename
    @Parameter(names = {"--inputHDFSDir"}, required = true)
    String inputHDFSDir;

    @Parameter(names = {"--outputHDFSDir"}, required = true)
    String ouptutHDFSDir;

    @Parameter(names = {"--targetIsize"}, required = true)
    int targetIsize;

    @Parameter(names = {"--targetIsizeSD"}, required = true)
    int targetIsizeSD;

    @Parameter(names = {"--isMatePairs"})
    boolean matePairs = false;

    public void run() throws IOException {
        runHadoopJob();
    }

    private void runHadoopJob() throws IOException {
        JobConf conf = new JobConf();

        conf.setJobName("Pileup Single End Deletion Score");
        conf.setJarByClass(SVPipeline.class);
        FileInputFormat.addInputPath(conf, new Path(inputHDFSDir));
        Path outputDir = new Path(ouptutHDFSDir);
        FileSystem.get(conf).delete(outputDir);

        FileOutputFormat.setOutputPath(conf, outputDir);

        conf.setInputFormat(TextInputFormat.class);

        conf.set("pileupDeletionScore.targetIsize", String.valueOf(targetIsize));
        conf.set("pileupDeletionScore.targetIsizeSD", String.valueOf(targetIsizeSD));
        conf.set("pileupDeletionScore.isMatePairs", String.valueOf(matePairs));

        conf.setMapperClass(SingleEndAlignmentsToDeletionScoreMapper.class);
        conf.setMapOutputKeyClass(Text.class);
        conf.setMapOutputValueClass(DoubleWritable.class);

        conf.setReducerClass(SingleEndDeletionScorePileupReducer.class);
        //conf.setReducerClass(IdentityReducer.class);

        conf.setOutputKeyClass(Text.class);
        conf.setOutputValueClass(DoubleWritable.class);
        conf.setCompressMapOutput(true);

        JobClient.runJob(conf);

    }
}
