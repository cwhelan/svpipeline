package edu.ohsu.sonmezsysbio.svpipeline.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.ohsu.sonmezsysbio.svpipeline.WigFileHelper;
import org.apache.hadoop.conf.Configuration;

import java.io.*;
import java.util.zip.GZIPInputStream;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 3/22/12
 * Time: 10:33 PM
 */
@Parameters(separators = "=", commandDescription = "Extract positive regions from a WIG file into a BED file")
public class CommandExtractPositiveRegionsFromWig implements SVPipelineCommand {

    @Parameter(names = {"--name"}, required = true)
    String name;
    
    @Parameter(names = {"--inputWigFile"}, required = true)
    String inputWigFile;

    @Parameter(names = {"--outputBedFile"}, required = true)
    String outputBedFile;

    public void run(Configuration conf) throws Exception {
        BufferedReader wigFileReader;
        if (inputWigFile.endsWith(".gz")) {
            wigFileReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(new File(inputWigFile)))));
        } else {
            wigFileReader = new BufferedReader(new FileReader(new File(inputWigFile)));
        }
        BufferedWriter bedFileWriter = new BufferedWriter(new FileWriter(new File(outputBedFile)));
        try {
            WigFileHelper.exportPositiveRegionsFromWig(name, wigFileReader, bedFileWriter);
        } finally {
            wigFileReader.close();
            bedFileWriter.close();
        }

    }
}
