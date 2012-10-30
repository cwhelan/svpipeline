package edu.ohsu.sonmezsysbio.cloudbreak.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.ohsu.sonmezsysbio.cloudbreak.Cloudbreak;
import edu.ohsu.sonmezsysbio.cloudbreak.mapper.SingleEndAlignmentsToBedSpansMapper;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;
import org.apache.hadoop.mapred.lib.IdentityReducer;

import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/23/11
 * Time: 10:02 AM
 */
@Parameters(separators = "=", commandDescription = "Dump all spanning read pairs with their deletion scores to BED format (debugging)")
public class CommandDumpReadsWithScores implements CloudbreakCommand {

    @Parameter(names = {"--inputHDFSDir"}, required = true)
    String inputHDFSDir;

    @Parameter(names = {"--outputHDFSDir"}, required = true)
    String ouptutHDFSDir;

    @Parameter(names = {"--region"}, required = true)
    String region;

    @Parameter(names = {"--isMatePairs"})
    boolean matePairs = false;

    @Parameter(names = {"--maxInsertSize"})
    int maxInsertSize = 500000;

    @Parameter(names = {"--aligner"})
    String aligner = Cloudbreak.ALIGNER_NOVOALIGN;

    public void run(Configuration configuration) throws IOException {
        runHadoopJob(configuration);
    }

    private void runHadoopJob(Configuration configuration) throws IOException {
        JobConf conf = new JobConf(configuration);

        conf.setJobName("Pileup Deletion Score");
        conf.setJarByClass(Cloudbreak.class);
        FileInputFormat.addInputPath(conf, new Path(inputHDFSDir));
        Path outputDir = new Path(ouptutHDFSDir);
        FileSystem.get(conf).delete(outputDir);

        FileOutputFormat.setOutputPath(conf, outputDir);

        conf.setInputFormat(SequenceFileInputFormat.class);

        conf.set("cloudbreak.aligner", aligner);
        conf.set("pileupDeletionScore.isMatePairs", String.valueOf(matePairs));
        conf.set("pileupDeletionScore.maxInsertSize", String.valueOf(maxInsertSize));
        conf.set("pileupDeletionScore.region", region);

        conf.setMapperClass(SingleEndAlignmentsToBedSpansMapper.class);
        conf.setMapOutputKeyClass(Text.class);
        conf.setMapOutputValueClass(Text.class);

        conf.setReducerClass(IdentityReducer.class);

        conf.setOutputKeyClass(Text.class);
        conf.setOutputValueClass(Text.class);

        JobClient.runJob(conf);

    }
}
