package edu.ohsu.sonmezsysbio.cloudbreak.command;

import com.beust.jcommander.Parameter;
import edu.ohsu.sonmezsysbio.cloudbreak.Cloudbreak;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapred.*;
import org.apache.hadoop.mapred.lib.IdentityMapper;
import org.apache.hadoop.mapred.lib.IdentityReducer;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 2/6/12
 * Time: 10:29 PM
 */
public class CommandSplitDeflatedOutput implements CloudbreakCommand {

    @Parameter(names = {"--inputHDFSDir"}, required = true)
    String inputHDFSDir;

    @Parameter(names = {"--outputHDFSDir"}, required = true)
    String ouptutHDFSDir;
        
    @Parameter(names = {"--numSplits"}, required = true)
    int numSplits;

    public void run(Configuration configuration) throws Exception {
        JobConf conf = new JobConf(configuration);

        conf.setJobName("Split deflated file");
        conf.setJarByClass(Cloudbreak.class);
        FileInputFormat.addInputPath(conf, new Path(inputHDFSDir));
        Path outputDir = new Path(ouptutHDFSDir);
        FileSystem.get(conf).delete(outputDir);

        FileOutputFormat.setOutputPath(conf, outputDir);

        conf.setInputFormat(TextInputFormat.class);

        conf.setMapperClass(IdentityMapper.class);

        //conf.setReducerClass(DeletionScorePileupReducer.class);
        conf.setReducerClass(IdentityReducer.class);
        conf.setNumReduceTasks(numSplits);
        conf.setCompressMapOutput(true);

        JobClient.runJob(conf);

    }
}
