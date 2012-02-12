package edu.ohsu.sonmezsysbio.svpipeline.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.ohsu.sonmezsysbio.svpipeline.reducer.NovoalignSingleEndAlignmentsToPairsReducer;
import edu.ohsu.sonmezsysbio.svpipeline.mapper.NovoalignSingleEndMapper;
import edu.ohsu.sonmezsysbio.svpipeline.SVPipeline;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.filecache.DistributedCache;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;

import java.io.*;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/18/11
 * Time: 2:01 PM
 */
@Parameters(separators = "=", commandDescription = "Run a novoalign mate pair alignment")
public class CommandNovoalignSingleEnds implements SVPipelineCommand {

    @Parameter(names = {"--HDFSDataDir"}, required = true)
    String hdfsDataDir;

    @Parameter(names = {"--HDFSAlignmentsDir"}, required = true)
    String hdfsAlignmentsDir;

    @Parameter(names = {"--reference"}, required = true)
    String reference;

    public void runHadoopJob() throws IOException, URISyntaxException {
        JobConf conf = new JobConf();

        conf.setJobName("Single End Alignment");
        conf.setJarByClass(SVPipeline.class);
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
        conf.set("mapred.output.compress", "true");

        conf.setMapperClass(NovoalignSingleEndMapper.class);
        conf.setMapOutputKeyClass(Text.class);
        conf.setMapOutputValueClass(Text.class);

        //todo: do I need this? Is there a way to configure hadoop to do it
        //conf.setNumMapTasks(numRecords / 60000 + 1);

        conf.setOutputKeyClass(Text.class);
        conf.setOutputValueClass(DoubleWritable.class);
        //conf.setCompressMapOutput(true);

        conf.setReducerClass(NovoalignSingleEndAlignmentsToPairsReducer.class);

        JobClient.runJob(conf);

    }

    public void run() throws Exception {
        runHadoopJob();
    }
}
