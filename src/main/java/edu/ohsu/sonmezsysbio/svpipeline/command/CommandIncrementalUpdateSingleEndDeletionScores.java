package edu.ohsu.sonmezsysbio.svpipeline.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.ohsu.sonmezsysbio.svpipeline.SVPipeline;
import edu.ohsu.sonmezsysbio.svpipeline.io.GenomicLocation;
import edu.ohsu.sonmezsysbio.svpipeline.io.ReadPairInfo;
import edu.ohsu.sonmezsysbio.svpipeline.mapper.SingleEndAlignmentsToReadPairInfoMapper;
import edu.ohsu.sonmezsysbio.svpipeline.reducer.IncrementalDelBeliefUpdateReadPairInfoReducer;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.filecache.DistributedCache;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.mapred.*;

import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/06/12
 * Time: 2:53 PM
 */
@Parameters(separators = "=", commandDescription = "Calculate Deletion Scores Across the Genome via Incremental Belief Update")
public class CommandIncrementalUpdateSingleEndDeletionScores implements SVPipelineCommand {

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

    @Parameter(names = {"--maxInsertSize"})
    int maxInsertSize = 500000;

    @Parameter(names = {"--faidx"}, required=true)
    String faidxFileName;

    @Parameter(names = {"--chrFilter"})
    String chrFilter;

    @Parameter(names = {"--startFilter"})
    Long startFilter;

    @Parameter(names = {"--endFilter"})
    Long endFilter;

    @Parameter(names = {"--excludePairsMappingIn"})
    String exclusionRegionsFileName;

    public void run(Configuration conf) throws IOException, URISyntaxException {
        runHadoopJob(conf);
    }

    private void runHadoopJob(Configuration configuration) throws IOException, URISyntaxException {
        JobConf conf = new JobConf(configuration);

        conf.setJobName("Incremental Update Single End Deletion Score");
        conf.setJarByClass(SVPipeline.class);
        FileInputFormat.addInputPath(conf, new Path(inputHDFSDir));
        Path outputDir = new Path(ouptutHDFSDir);
        FileSystem.get(conf).delete(outputDir);

        FileOutputFormat.setOutputPath(conf, outputDir);


        File faidxFile = new File(faidxFileName);
        String faidxFileBasename = faidxFile.getName();
        String faidxFileDir = faidxFile.getParent();

        DistributedCache.addCacheFile(new URI(faidxFileDir + "/" + faidxFileBasename + "#" + faidxFileBasename),
                conf);

        conf.set("alignment.faidx", faidxFileBasename);

        if (exclusionRegionsFileName != null) {
            File exclusionRegionsFile = new File(exclusionRegionsFileName);
            String exclusionRegionsFileBasename = exclusionRegionsFile.getName();
            String exclusionRegionsFileDir = exclusionRegionsFile.getParent();

            DistributedCache.addCacheFile(new URI(exclusionRegionsFileDir + "/" + exclusionRegionsFileBasename + "#" + exclusionRegionsFileBasename),
                    conf);

            conf.set("alignment.exclusionRegions", exclusionRegionsFileBasename);
        }

        DistributedCache.createSymlink(conf);

        conf.set("pileupDeletionScore.targetIsize", String.valueOf(targetIsize));
        conf.set("pileupDeletionScore.targetIsizeSD", String.valueOf(targetIsizeSD));
        conf.set("pileupDeletionScore.isMatePairs", String.valueOf(matePairs));
        conf.set("pileupDeletionScore.maxInsertSize", String.valueOf(maxInsertSize));

        if (chrFilter != null) {
            conf.set("alignments.filterchr", chrFilter);
            conf.set("alignments.filterstart", startFilter.toString());
            conf.set("alignments.filterend", endFilter.toString());
        }

        conf.setInputFormat(TextInputFormat.class);

        conf.setMapperClass(SingleEndAlignmentsToReadPairInfoMapper.class);
        conf.setMapOutputKeyClass(GenomicLocation.class);
        conf.setMapOutputValueClass(ReadPairInfo.class);

        conf.setReducerClass(IncrementalDelBeliefUpdateReadPairInfoReducer.class);
        //conf.setReducerClass(IdentityReducer.class);

        conf.setOutputKeyClass(GenomicLocation.class);
        conf.setOutputValueClass(DoubleWritable.class);
        conf.setOutputFormat(SequenceFileOutputFormat.class);

        conf.setCompressMapOutput(true);

        JobClient.runJob(conf);

    }
}
