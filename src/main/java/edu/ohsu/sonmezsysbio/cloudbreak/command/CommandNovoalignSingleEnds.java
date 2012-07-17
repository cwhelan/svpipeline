package edu.ohsu.sonmezsysbio.cloudbreak.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.ohsu.sonmezsysbio.cloudbreak.Cloudbreak;
import edu.ohsu.sonmezsysbio.cloudbreak.reducer.NovoalignSingleEndAlignmentsToPairsReducer;
import edu.ohsu.sonmezsysbio.cloudbreak.mapper.NovoalignSingleEndMapper;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.filecache.DistributedCache;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;

import java.io.*;
import java.net.URI;
import java.net.URISyntaxException;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/18/11
 * Time: 2:01 PM
 */
@Parameters(separators = "=", commandDescription = "Run a novoalign mate pair alignment")
public class CommandNovoalignSingleEnds implements CloudbreakCommand {

    @Parameter(names = {"--HDFSDataDir"}, required = true)
    String hdfsDataDir;

    @Parameter(names = {"--HDFSAlignmentsDir"}, required = true)
    String hdfsAlignmentsDir;

    @Parameter(names = {"--reference"}, required = true)
    String reference;
    
    @Parameter(names = {"--threshold"}, required = true)
    String threshold;

    @Parameter(names = {"--qualityFormat"})
    String qualityFormat = "ILMFQ";

    public void runHadoopJob(Configuration configuration) throws IOException, URISyntaxException {
        JobConf conf = new JobConf(configuration);

        conf.setJobName("Single End Alignment");
        conf.setJarByClass(Cloudbreak.class);
        FileInputFormat.addInputPath(conf, new Path(hdfsDataDir));
        Path outputDir = new Path(hdfsAlignmentsDir);
        FileSystem.get(conf).delete(outputDir);

        FileOutputFormat.setOutputPath(conf, outputDir);

        conf.setInputFormat(TextInputFormat.class);

        File referenceFile = new File(reference);
        String referenceBasename = referenceFile.getName();
        String referenceDir = referenceFile.getParent();

        DistributedCache.addCacheFile(new URI(referenceDir + "/" + referenceBasename + "#" + referenceBasename),
                conf);
        DistributedCache.createSymlink(conf);
        conf.set("mapred.task.timeout", "3600000");
        conf.set("novoalign.reference", referenceBasename);
        conf.set("novoalign.threshold", threshold);
        conf.set("novoalign.quality.format", qualityFormat);

        conf.setMapperClass(NovoalignSingleEndMapper.class);
        conf.setMapOutputKeyClass(Text.class);
        conf.setMapOutputValueClass(Text.class);

        //todo: do I need this? Is there a way to configure hadoop to do it
        //conf.setNumMapTasks(numRecords / 60000 + 1);

        conf.setOutputKeyClass(Text.class);
        conf.setOutputValueClass(DoubleWritable.class);
        conf.setCompressMapOutput(true);

        conf.setReducerClass(NovoalignSingleEndAlignmentsToPairsReducer.class);

        JobClient.runJob(conf);

    }

    public void run(Configuration conf) throws Exception {
        runHadoopJob(conf);
    }
}
