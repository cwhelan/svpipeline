package edu.ohsu.sonmezsysbio.svpipeline.file;

import edu.ohsu.sonmezsysbio.svpipeline.file.FaidxFileHelper;
import edu.ohsu.sonmezsysbio.svpipeline.file.WigFileHelper;
import org.junit.Test;

import java.io.*;

import static junit.framework.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/16/12
 * Time: 10:06 AM
 */
public class WigFileHelperTest {
    @Test
    public void testExportPositiveRegionsFromWig() throws Exception {

        FaidxFileHelper faidx = new FaidxFileHelper("foo") {
            @Override
            public Long getLengthForChromName(String name) throws IOException {
                if ("chr1".equals(name)) return 4500l;
                if ("chr2".equals(name)) return 5000l;
                return null;
            }
        };

        // test chromosome transition
        String outputPrefix = "test";
        String wigFile = "variableStep chrom=chr1 span=1000\n" +
                "1000\t2.0\n" +
                "2000\t4.0\n" +
                "3000\t5.0\n" +
                "variableStep chrom=chr2 span=1000\n" +
                "1000\t7.0\n" +
                "2000\t6.0\n" +
                "3000\t3.0\n";
        BufferedReader wigFileReader = new BufferedReader(new StringReader(wigFile));

        StringWriter stringWriter = new StringWriter();
        BufferedWriter bedFileWriter = new BufferedWriter(stringWriter);
        double threshold = 4.5;

        WigFileHelper.exportPositiveRegionsFromWig(outputPrefix, wigFileReader, bedFileWriter, threshold, faidx, 1);
        bedFileWriter.close();

        String expectedOutput =
                "track name = \"test peaks over " + threshold + "\"\n" +
                "chr1\t3000\t3999\t1\t0\n" +
                "chr2\t1000\t2999\t2\t0\n";
        assertEquals(expectedOutput, stringWriter.getBuffer().toString());
    }

    @Test
    public void testExportPositiveRegionsFromWigWithMedianFilter() throws Exception {
        FaidxFileHelper faidx = new FaidxFileHelper("foo") {
            @Override
            public Long getLengthForChromName(String name) throws IOException {
                if ("chr1".equals(name)) return 6500l;
                return null;
            }
        };
        String outputPrefix = "test";
        String wigFile = "variableStep chrom=chr1 span=1000\n" +
                "1000\t2.0\n" +
                "2000\t5.0\n" +
                "3000\t5.0\n" +
                "4000\t3.0\n" +
                "5000\t6.0\n" +
                "6000\t3.0\n";
        BufferedReader wigFileReader = new BufferedReader(new StringReader(wigFile));

        StringWriter stringWriter = new StringWriter();
        BufferedWriter bedFileWriter = new BufferedWriter(stringWriter);
        double threshold = 4.5;

        WigFileHelper.exportPositiveRegionsFromWig(outputPrefix, wigFileReader, bedFileWriter, threshold, faidx, 3);
        bedFileWriter.close();

        String expectedOutput =
                "track name = \"test peaks over " + threshold + "\"\n" +
                        "chr1\t2000\t4999\t1\t0\n";
        assertEquals(expectedOutput, stringWriter.getBuffer().toString());

    }
}