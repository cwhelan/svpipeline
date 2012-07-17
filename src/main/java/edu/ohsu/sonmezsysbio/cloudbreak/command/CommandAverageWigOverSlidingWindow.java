package edu.ohsu.sonmezsysbio.cloudbreak.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.ohsu.sonmezsysbio.cloudbreak.SVPipeline;
import edu.ohsu.sonmezsysbio.cloudbreak.file.WigFileHelper;
import org.apache.hadoop.conf.Configuration;

import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 6/6/11
 * Time: 3:42 PM
 */
@Parameters(separators = "=", commandDescription = "compute a transformed wig file averaged over a sliding window")
public class CommandAverageWigOverSlidingWindow implements CloudbreakCommand {

    @Parameter(names = {"--InFile"}, required = true)
    String inFile;

    @Parameter(names = {"--OutFile"}, required = true)
    String outFile;

    @Parameter(names = {"--resolution"})
    int resolution = SVPipeline.DEFAULT_RESOLUTION;

    public void run(Configuration conf) throws IOException {
        BufferedReader inFileReader = new BufferedReader(new FileReader(new File(inFile)));
        BufferedWriter outFileWriter = new BufferedWriter(new FileWriter(new File(outFile)));
        try {
            WigFileHelper.averageWigOverSlidingWindow(resolution, SVPipeline.WINDOW_SIZE_IN_LINES, inFileReader, outFileWriter);
        } finally {
            inFileReader.close();
            outFileWriter.close();
        }

    }


}
