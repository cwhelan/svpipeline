package edu.ohsu.sonmezsysbio.svpipeline.command;

import com.beust.jcommander.Parameter;
import edu.ohsu.sonmezsysbio.svpipeline.SVPipeline;
import edu.ohsu.sonmezsysbio.svpipeline.mapper.SingleEndAlignmentsToDeletionScoreMapper;
import edu.ohsu.sonmezsysbio.svpipeline.reducer.SingleEndDeletionScorePileupReducer;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;
import org.apache.hadoop.mapred.lib.KeyFieldBasedComparator;

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

    private void runHadoopJob(Configuration configuration) {
        JobConf conf = new JobConf(configuration);

        conf.setJobName("Pileup Single End Deletion Score");
        conf.setJarByClass(SVPipeline.class);
        FileInputFormat.addInputPath(conf, new Path(inputHDFSDir));
        Path outputDir = new Path(outputHDFSDir);
        FileSystem.get(conf).delete(outputDir);

        FileOutputFormat.setOutputPath(conf, outputDir);

        conf.setInputFormat(TextInputFormat.class);

        conf.set("pileupDeletionScore.targetIsize", String.valueOf(targetIsize));
        conf.set("pileupDeletionScore.targetIsizeSD", String.valueOf(targetIsizeSD));
        conf.set("pileupDeletionScore.isMatePairs", String.valueOf(matePairs));
        conf.set("pileupDeletionScore.maxInsertSize", String.valueOf(maxInsertSize));

        conf.setMapperClass(SingleEndAlignmentsToDeletionScoreMapper.class);
        conf.setMapOutputKeyClass(Text.class);
        conf.setMapOutputValueClass(DoubleWritable.class);

        conf.setCombinerClass(SingleEndDeletionScorePileupReducer.class);

        conf.setReducerClass(SingleEndDeletionScorePileupReducer.class);
        //conf.setReducerClass(IdentityReducer.class);

        conf.setOutputKeyClass(Text.class);
        conf.setOutputValueClass(DoubleWritable.class);
        conf.setKeyFieldComparatorOptions("-k 1,1 -k 2,2n");
        conf.setOutputKeyComparatorClass(KeyFieldBasedComparator.class);

        conf.setCompressMapOutput(true);

        JobClient.runJob(conf);

    }
}
