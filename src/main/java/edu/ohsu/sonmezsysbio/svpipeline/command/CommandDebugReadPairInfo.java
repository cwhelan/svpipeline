package edu.ohsu.sonmezsysbio.svpipeline.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.ohsu.sonmezsysbio.svpipeline.SVPipeline;
import edu.ohsu.sonmezsysbio.svpipeline.io.GenomicLocation;
import edu.ohsu.sonmezsysbio.svpipeline.io.ReadPairInfo;
import edu.ohsu.sonmezsysbio.svpipeline.mapper.SingleEndAlignmentsToReadPairInfoMapper;
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
        DistributedCache.createSymlink(conf);

        conf.setInputFormat(TextInputFormat.class);
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

        JobClient.runJob(conf);

    }
}
