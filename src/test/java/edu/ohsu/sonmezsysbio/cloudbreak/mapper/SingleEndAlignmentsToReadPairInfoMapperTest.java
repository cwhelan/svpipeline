package edu.ohsu.sonmezsysbio.cloudbreak.mapper;

import edu.ohsu.sonmezsysbio.cloudbreak.Cloudbreak;
import edu.ohsu.sonmezsysbio.cloudbreak.file.FaidxFileHelper;
import edu.ohsu.sonmezsysbio.cloudbreak.file.GFFFileHelper;
import edu.ohsu.sonmezsysbio.cloudbreak.ProbabilisticPairedAlignmentScorer;
import edu.ohsu.sonmezsysbio.cloudbreak.io.GenomicLocationWithQuality;
import edu.ohsu.sonmezsysbio.cloudbreak.io.NovoalignAlignmentReader;
import edu.ohsu.sonmezsysbio.svpipeline.io.GenomicLocation;
import edu.ohsu.sonmezsysbio.cloudbreak.io.ReadPairInfo;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.OutputCollector;
import org.apache.hadoop.mapred.Reporter;
import org.junit.Before;
import org.junit.Test;
import org.mockito.Mockito;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/6/12
 * Time: 4:58 PM
 */
public class SingleEndAlignmentsToReadPairInfoMapperTest {

    private SingleEndAlignmentsToReadPairInfoMapper mapper;

    @Before
    public void setup() {
        mapper = new SingleEndAlignmentsToReadPairInfoMapper();
        mapper.setScorer(new ProbabilisticPairedAlignmentScorer());
        mapper.setAlignmentReader(new NovoalignAlignmentReader());
        mapper.setFaix(new FaidxFileHelper("foo") {
            @Override
            public Short getKeyForChromName(String name) throws IOException {
                assertEquals("10", name);
                return (short) 9;
            }
        });

    }

    @Test
    public void testMapPairedEnd() throws Exception {

        String key = "@ERR000545.10000001 EAS139_44:1:93:532:453";
        String value = "@ERR000545.10000001 EAS139_44:1:93:532:453/1\tS\tCAAAAACCACTTGTACTCCAAAAGCTATTGAAGTTTAAGTTAAAATAAAAA\t<??>>?;<>=@?=?>8@<<9=98=:@>>>=:>>:6?7>9:?<46:;9;.:9\tR\t114\t0.00000\t>10\t43049466\tR\t.\t.\t.\t10A>C 22C>A 46-A\tSVP_READ\t@ERR000545.10000001 EAS139_44:1:93:532:453/2\tS\tTTATTGCACTTACCATGACTGTCTTCTGAAATGCATCTCAACCCTTGAATA\t;<8>>:$>=?@?>>:>=>:9=8<>1;8:<=>>9=:=7><>=;=;=>=72:>\tU\t34\t145.2313\t>10\t43039500\tF\t.\t.\t.\t3G>A";

        mapper.setReadGroupId((short) 3);

        MockOutputCollector collector = new MockOutputCollector();
        Reporter reporter = Mockito.mock(Reporter.class);

        mapper.map(new Text(key), new Text(value), collector, reporter);

        int idx = 0;
        for (int i = 43039500; i <= 43049500; i = i + Cloudbreak.DEFAULT_RESOLUTION) {
            assertEquals((short) 9, collector.keys.get(idx).chromosome);
            assertEquals(i, collector.keys.get(idx).pos);

            assertEquals(10017, collector.values.get(idx).insertSize);
            assertEquals(-9.2103, collector.values.get(idx).pMappingCorrect, .0001);
            assertEquals((short) 3, (short) collector.values.get(idx).readGroupId);
            idx++;
        }

        assertEquals(idx, collector.keys.size());
//
//        String complexInputLine ="@ERR000545.10000241 EAS139_44:1:93:535:874\t@ERR000545.10000241 EAS139_44:1:93:535:874/1\tS\tTAGGTAATGTTTGGGAGGGAGTTGTTGGTTTTGTTGATTTATTATATCTTG\t??==>79>=?@@==>9>>=5>:?<>?<<<>??>;?=<=>?7?>;?8=4>>8\tR\t0\t60\t>10\t81550461\tR\t.\t.\t.\t" +
//                "SVP_ALIGNMENT\t@ERR000545.10000241 EAS139_44:1:93:535:874/1\tS\tTAGGTAATGTTTGGGAGGGAGTTGTTGGTTTTGTTGATTTATTATATCTTG\t??==>79>=?@@==>9>>=5>:?<>?<<<>??>;?=<=>?7?>;?8=4>>8\tR\t60\t0\t>10\t89072620\tR\t.\t.\t.\t48A>C 49T>C\t" +
//                "SVP_READ\t@ERR000545.10000241 EAS139_44:1:93:535:874/2\tS\tTATTATATTTTTTAATTTGACAGAGTAGTGCAGGCAATAATGAAATGGTAT\t<>==>A?@=?@A@>=??>>;>>;-<==;:<<<:>;==>>8?<::<>;80;>\tR\t0\t3\t>10\t89072431\tF\t.\t.\t.\t" +
//                "SVP_ALIGNMENT\t@ERR000545.10000241 EAS139_44:1:93:535:874/2\tS\tTATTATATTTTTTAATTTGACAGAGTAGTGCAGGCAATAATGAAATGGTAT\t<>==>A?@=?@A@>=??>>;>>;-<==;:<<<:>;==>>8?<::<>;80;>\tR\t0\t3\t>10\t81550273\tF\t.\t.\t.";
//
//        collector = Mockito.mock(OutputCollector.class);
//        reporter = Mockito.mock(Reporter.class);
//
//        mapper.map(new LongWritable(1), new Text(complexInputLine), collector, reporter);
//
//        for (int i = 81550200; i <= 81550500; i = i + Cloudbreak.DEFAULT_RESOLUTION) {
//            verify(collector).collect(new Text("10\t" + i),
//                    new DoubleWritable(1));
//        }
//
//        for (int i = 89072400; i <= 89072700; i = i + Cloudbreak.DEFAULT_RESOLUTION) {
//            verify(collector).collect(new Text("10\t" + i),
//                    new DoubleWritable(1));
//        }
//
//        verifyNoMoreInteractions(collector);
//

    }

    @Test
    public void testMapPairedEndWithSegmentalDuplications() throws Exception {

        String key = "@ERR000545.10000001 EAS139_44:1:93:532:453";
        String value = "@ERR000545.10000001 EAS139_44:1:93:532:453/1\tS\tCAAAAACCACTTGTACTCCAAAAGCTATTGAAGTTTAAGTTAAAATAAAAA\t<??>>?;<>=@?=?>8@<<9=98=:@>>>=:>>:6?7>9:?<46:;9;.:9\tR\t114\t0.00000\t>10\t43049466\tR\t.\t.\t.\t10A>C 22C>A 46-A\tSVP_READ\t@ERR000545.10000001 EAS139_44:1:93:532:453/2\tS\tTTATTGCACTTACCATGACTGTCTTCTGAAATGCATCTCAACCCTTGAATA\t;<8>>:$>=?@?>>:>=>:9=8<>1;8:<=>>9=:=7><>=;=;=>=72:>\tU\t34\t145.2313\t>10\t43039500\tF\t.\t.\t.\t3G>A";

        mapper.setExclusionRegions(new GFFFileHelper() {
            @Override
            public boolean doesLocationOverlap(String chrom, int start, int end) throws Exception {
                if ("10".equals(chrom) &&
                        ((start > 43039000 && start < 43040000) && (end > 43039000 && end < 43040000)) ||
                        ((start > 43049000 && start < 43050000) && (end > 43049000 && end < 43050000))) {
                    return true;
                }
                return false;
            }
        });

        MockOutputCollector collector = new MockOutputCollector();
        Reporter reporter = Mockito.mock(Reporter.class);

        mapper.map(new Text(key), new Text(value), collector, reporter);

        assertEquals(0, collector.keys.size());

        key = "@ERR000545.10000001 EAS139_44:1:93:532:453";
        value = "@ERR000545.10000001 EAS139_44:1:93:532:453/1\tS\tCAAAAACCACTTGTACTCCAAAAGCTATTGAAGTTTAAGTTAAAATAAAAA\t<??>>?;<>=@?=?>8@<<9=98=:@>>>=:>>:6?7>9:?<46:;9;.:9\tR\t114\t0.00000\t>10\t49466\tR\t.\t.\t.\t10A>C 22C>A 46-A\tSVP_READ\t@ERR000545.10000001 EAS139_44:1:93:532:453/2\tS\tTTATTGCACTTACCATGACTGTCTTCTGAAATGCATCTCAACCCTTGAATA\t;<8>>:$>=?@?>>:>=>:9=8<>1;8:<=>>9=:=7><>=;=;=>=72:>\tU\t34\t145.2313\t>10\t39500\tF\t.\t.\t.\t3G>A";
        collector = new MockOutputCollector();
        mapper.map(new Text(key), new Text(value), collector, reporter);
        assertEquals(101, collector.keys.size());
    }

    private static class MockOutputCollector implements OutputCollector<GenomicLocationWithQuality, ReadPairInfo> {

        List<GenomicLocationWithQuality> keys = new ArrayList<GenomicLocationWithQuality>();
        List<ReadPairInfo> values = new ArrayList<ReadPairInfo>();

        public void collect(GenomicLocationWithQuality key, ReadPairInfo value) throws IOException {
            keys.add(key);
            values.add(value);
        }

        public void reset() {
            keys.clear();
            values.clear();
        }
    }

    @Test
    public void testMap_real1() throws Exception {
        String key = "@2_132385096_132385222_0_1_0_0_1:0:0_1:0:0_1bfnjd/";
        String val = "@2_132385096_132385222_0_1_0_0_1:0:0_1:0:0_1bfnjd//1\tS\tTAAAAAGCCGCGGCGACTAAAAGCCGCTGAGAGGGGGCAAAAAGCAGCGG\t66554410000////1.0000/----,/,.,.,,++++-----+*-****\tU\t25\t139.09267\t>2\t132512583\tF\t.\t.\t.\t18A>T\tSVP_READ\t@2_132385096_132385222_0_1_0_0_1:0:0_1:0:0_1bfnjd//2\tS\tCCCCTGCCCCGCCGCGGCTTTTTGCGGCTTTCCGCCCCGGCCGCCGCGGA\t33324110000////...000//---,,...,,,++++++++++*****,\tU\t23\t135.68983\t>2\t132512814\tR\t.\t.\t.\t1G>T";
        String va2 = "@2_132385096_132385222_0_1_0_0_1:0:0_1:0:0_1bf17d//1\tS\tTAAAAAGCCGCGGCGACTAAAAGCCGCTGAGAGGGGGCAAAAAGCAGCGG\t66554410000////1.0000/----,/,.,.,,++++-----+*-****\tU\t25\t139.09267\t>2\t132512583\tF\t.\t.\t.\t18A>T\tSVP_READ\t@2_132385096_132385222_0_1_0_0_1:0:0_1:0:0_1bf17d//2\tS\tCCCCTGCCCCGCCGCGGCTTTTTGCGGCTTTCCGCCCCGGCCGCCGCGGA\t33324110000////...000//---,,...,,,++++++++++*****,\t23\t135.68983\t>2\t132512814\tR\t.\t.\t.\t1G>T";

        mapper.setFaix(new FaidxFileHelper("foo") {
            @Override
            public Short getKeyForChromName(String name) throws IOException {
                assertEquals("2", name);
                return (short) 0;
            }
        });

        MockOutputCollector mockOutputCollector = new MockOutputCollector();
        mapper.map(new Text(key), new Text(val), mockOutputCollector, null);
        mapper.setChromosomeFilter("2");
        mapper.setStartFilter(132512600l);
        mapper.setEndFilter(132512800l);

        GenomicLocationWithQuality gl132512700 = new GenomicLocationWithQuality();
        gl132512700.chromosome = 0;
        gl132512700.pos = 132512700;
        gl132512700.pMappingCorrect = -3.9301895071730983E-14;

        assertTrue(mockOutputCollector.keys.contains(gl132512700));
        assertEquals(281, mockOutputCollector.values.get(0).insertSize);
    }

    @Test
    public void testGetInputPath() throws Exception {
        assertEquals("/user/whelanch/cloudbreak/jcvi_chr2_lc/se_alignments_t180/part-00000",
                SingleEndAlignmentsToReadPairInfoMapper.getInputPath("hdfs://bigbird51.csee.ogi.edu:50030/user/whelanch/cloudbreak/jcvi_chr2_lc/se_alignments_t180/part-00000"));
    }
}
