package edu.ohsu.sonmezsysbio.svpipeline.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.ohsu.sonmezsysbio.svpipeline.ReadGroupInfo;
import edu.ohsu.sonmezsysbio.svpipeline.SVPipeline;
import edu.ohsu.sonmezsysbio.svpipeline.file.DFSFacade;
import edu.ohsu.sonmezsysbio.svpipeline.file.ReadGroupInfoFileHelper;
import edu.ohsu.sonmezsysbio.svpipeline.io.GenomicLocation;
import edu.ohsu.sonmezsysbio.svpipeline.io.ReadPairInfo;
import edu.ohsu.sonmezsysbio.svpipeline.mapper.SingleEndAlignmentsToReadPairInfoMapper;
import edu.ohsu.sonmezsysbio.svpipeline.partitioner.GenomicLocationPartitioner;
import edu.ohsu.sonmezsysbio.svpipeline.reducer.IncrementalDelBeliefUpdateReadPairInfoReducer;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.filecache.DistributedCache;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hdfs.DistributedFileSystem;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.mapred.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/06/12
 * Time: 2:53 PM
 */
@Parameters(separators = "=", commandDescription = "Calculate Deletion Scores Across the Genome via Incremental Belief Update")
public class CommandIncrementalUpdateSingleEndDeletionScores implements SVPipelineCommand {

    @Parameter(names = {"--inputFileDescriptor"}, required = true)
    String inputFileDescriptor;

    @Parameter(names = {"--outputHDFSDir"}, required = true)
    String outputHDFSDir;

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

    @Parameter(names = {"--resolution"})
    int resolution = SVPipeline.DEFAULT_RESOLUTION;

    @Parameter(names = {"--mapabilityWeighting"})
    String mapabilityWeightingFileName;

    public void run(Configuration conf) throws IOException, URISyntaxException {
        runHadoopJob(conf);
    }

    private void runHadoopJob(Configuration configuration) throws IOException, URISyntaxException {
        JobConf conf = new JobConf(configuration);

        conf.setJobName("Incremental Update Single End Deletion Score");
        conf.setJarByClass(SVPipeline.class);

        ReadGroupInfoFileHelper readGroupInfoFileHelper = new ReadGroupInfoFileHelper();
        FileSystem dfs = DistributedFileSystem.get(conf);

        Map<Short, ReadGroupInfo> readGroupInfoMap = readGroupInfoFileHelper.readReadGroupsById(
                new BufferedReader(
                        new InputStreamReader(
                                new DFSFacade(dfs, conf).openPath(new Path(inputFileDescriptor)))));

        for (ReadGroupInfo readGroupInfo : readGroupInfoMap.values()) {
            FileInputFormat.addInputPath(conf, new Path(readGroupInfo.hdfsPath));
        }

        Path outputDir = new Path(outputHDFSDir);

        FileSystem.get(conf).delete(outputDir);

        FileOutputFormat.setOutputPath(conf, outputDir);

        setupDistributedCacheFile(conf, inputFileDescriptor, "read.group.info.file");
        setupDistributedCacheFile(conf, faidxFileName, "alignment.faidx");

        if (exclusionRegionsFileName != null) {
            setupDistributedCacheFile(conf, exclusionRegionsFileName, "alignment.exclusionRegions");
        }

        if (mapabilityWeightingFileName != null) {
            setupDistributedCacheFile(conf, mapabilityWeightingFileName, "alignment.mapabilityWeighting");
        }

        DistributedCache.createSymlink(conf);

        conf.set("svpipeline.resolution", String.valueOf(resolution));

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
        conf.setPartitionerClass(GenomicLocationPartitioner.class);

        conf.setReducerClass(IncrementalDelBeliefUpdateReadPairInfoReducer.class);
        //conf.setReducerClass(IdentityReducer.class);

        conf.setOutputKeyClass(GenomicLocation.class);
        conf.setOutputValueClass(DoubleWritable.class);
        conf.setOutputFormat(SequenceFileOutputFormat.class);

        conf.setCompressMapOutput(true);

        JobClient.runJob(conf);

    }

    private void setupDistributedCacheFile(JobConf conf, String fileName, String confPropertyName) throws URISyntaxException {
        File file = new File(fileName);
        String fileBasename = file.getName();
        String fileDir = file.getParent();

        DistributedCache.addCacheFile(new URI(fileDir + "/" + fileBasename + "#" + fileBasename),
                conf);

        conf.set(confPropertyName, fileBasename);
    }
}
