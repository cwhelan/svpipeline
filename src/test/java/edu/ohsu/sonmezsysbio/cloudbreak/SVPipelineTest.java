package edu.ohsu.sonmezsysbio.cloudbreak;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import org.junit.Test;

import static org.junit.Assert.*;


/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/18/11
 * Time: 2:12 PM
 */
public class SVPipelineTest {
    @Test
    public void testBuildJCommander() throws Exception {
        String[] argv = { "novoalignMatePair" };
        JCommander jc = SVPipeline.buildJCommander();
        try {
            jc.parse(argv);
            fail("validation should have failed because no files provided");
        } catch (ParameterException e) {
        }

        String[] argv2 = { "novoalignMatePair",
                "--HDFSDir", "/user/whelanch/tmp/svpipelinetest",
                "--fastqFile1", "s_1_1_sequence.txt",
                "--fastqFile2", "s_1_2_sequence.txt",
                "--reference", "ref",
                "--repeatReport", "Random",
                "--targetIsize", "3000",
                "--targetIsizeSD", "300",
                "--libraryName", "TC1"};
        jc.parse(argv2);

    }
}
