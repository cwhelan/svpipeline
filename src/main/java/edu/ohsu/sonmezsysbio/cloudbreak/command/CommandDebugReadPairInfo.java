package edu.ohsu.sonmezsysbio.cloudbreak.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.ohsu.sonmezsysbio.cloudbreak.SVPipeline;
import edu.ohsu.sonmezsysbio.cloudbreak.io.GenomicLocation;
import edu.ohsu.sonmezsysbio.cloudbreak.io.ReadPairInfo;
import edu.ohsu.sonmezsysbio.cloudbreak.mapper.SingleEndAlignmentsToReadPairInfoMapper;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.filecache.DistributedCache;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapred.*;
import org.apache.hadoop.mapred.lib.IdentityReducer;

import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/9/12
 * Time: 11:12 AM
 */
@Parameters(separators = "=", commandDescription = "View read pair infos")
public class CommandDebugReadPairInfo implements SVPipelineCommand {

    @Parameter(names = {"--inputHDFSDir"}, required = true)
    String inputHDFSDir;

    @Parameter(names = {"--outputHDFSDir"}, required = true)
    String ouptutHDFSDir;

    @Parameter(names = {"--maxInsertSize"})
    int maxInsertSize = 500000;

    @Parameter(names = {"--faidx"}, required=true)
    String faidxFileName;

    @Parameter(names = {"--chrFilter"}, required = true)
    String chrFilter;

    @Parameter(names = {"--startFilter"}, required = true)
    Long startFilter;

    @Parameter(names = {"--endFilter"}, required = true)
    Long endFilter;

    @Parameter(names = {"--resolution"})
    final int resolution = SVPipeline.DEFAULT_RESOLUTION;

    @Parameter(names = {"--excludePairsMappingIn"})
    String exclusionRegionsFileName;

    @Parameter(names = {"--mapabilityWeighting"})
    String mapabilityWeightingFileName;


    public void run(Configuration conf) throws Exception {
        runHadoopJob(conf);
    }

    private void runHadoopJob(Configuration configuration) throws IOException, URISyntaxException {
        JobConf conf = new JobConf(configuration);

        conf.setJobName("Debug Read Pair Info");
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

        if (exclusionRegionsFileName != null) {
            File exclusionRegionsFile = new File(exclusionRegionsFileName);
            String exclusionRegionsFileBasename = exclusionRegionsFile.getName();
            String exclusionRegionsFileDir = exclusionRegionsFile.getParent();

            URI uri = new URI(exclusionRegionsFileDir + "/" + exclusionRegionsFileBasename + "#" + exclusionRegionsFileBasename);
            System.err.println("URI: " + uri);
            DistributedCache.addCacheFile(uri,
                    conf);

            conf.set("alignment.exclusionRegions", exclusionRegionsFileBasename);
        }

        if (mapabilityWeightingFileName != null) {
            File mapabilityWeightingFile = new File(mapabilityWeightingFileName);
            String mapabilityWeightingFileBasename = mapabilityWeightingFile.getName();
            String mapabilityWeightingFileDir = mapabilityWeightingFile.getParent();

            URI uri = new URI(mapabilityWeightingFileDir + "/" + mapabilityWeightingFileBasename + "#" + mapabilityWeightingFileBasename);
            System.err.println("URI: " + uri);
            DistributedCache.addCacheFile(uri,
                    conf);

            conf.set("alignment.mapabilityWeighting", mapabilityWeightingFileBasename);
        }

        DistributedCache.createSymlink(conf);

        conf.setInputFormat(TextInputFormat.class);

        conf.set("cloudbreak.resolution", String.valueOf(resolution));

        conf.set("alignment.faidx", faidxFileBasename);
        conf.set("pileupDeletionScore.maxInsertSize", String.valueOf(maxInsertSize));
        conf.set("alignments.filterchr", chrFilter);
        conf.set("alignments.filterstart", startFilter.toString());
        conf.set("alignments.filterend", endFilter.toString());

        conf.setMapperClass(SingleEndAlignmentsToReadPairInfoMapper.class);
        conf.setMapOutputKeyClass(GenomicLocation.class);
        conf.setMapOutputValueClass(ReadPairInfo.class);

        conf.setReducerClass(IdentityReducer.class);

        conf.setOutputKeyClass(GenomicLocation.class);
        conf.setOutputValueClass(ReadPairInfo.class);
        conf.setOutputFormat(SequenceFileOutputFormat.class);

        JobClient.runJob(conf);

    }
}
