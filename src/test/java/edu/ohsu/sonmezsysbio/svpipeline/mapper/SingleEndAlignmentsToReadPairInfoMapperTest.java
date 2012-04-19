package edu.ohsu.sonmezsysbio.svpipeline.mapper;

import edu.ohsu.sonmezsysbio.svpipeline.FaidxFileHelper;
import edu.ohsu.sonmezsysbio.svpipeline.GFFFileHelper;
import edu.ohsu.sonmezsysbio.svpipeline.ProbabilisticPairedAlignmentScorer;
import edu.ohsu.sonmezsysbio.svpipeline.SVPipeline;
import edu.ohsu.sonmezsysbio.svpipeline.io.GenomicLocation;
import edu.ohsu.sonmezsysbio.svpipeline.io.ReadPairInfo;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.OutputCollector;
import org.apache.hadoop.mapred.Reporter;
import org.junit.Test;
import org.mockito.Mockito;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import static org.junit.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/6/12
 * Time: 4:58 PM
 */
public class SingleEndAlignmentsToReadPairInfoMapperTest {

    @Test
    public void testMapPairedEnd() throws Exception {

        String inputLine = "@ERR000545.10000001 EAS139_44:1:93:532:453\t@ERR000545.10000001 EAS139_44:1:93:532:453/1\tS\tCAAAAACCACTTGTACTCCAAAAGCTATTGAAGTTTAAGTTAAAATAAAAA\t<??>>?;<>=@?=?>8@<<9=98=:@>>>=:>>:6?7>9:?<46:;9;.:9\tR\t114\t0.00000\t>10\t43049466\tR\t.\t.\t.\t10A>C 22C>A 46-A\tSVP_READ\t@ERR000545.10000001 EAS139_44:1:93:532:453/2\tS\tTTATTGCACTTACCATGACTGTCTTCTGAAATGCATCTCAACCCTTGAATA\t;<8>>:$>=?@?>>:>=>:9=8<>1;8:<=>>9=:=7><>=;=;=>=72:>\tU\t34\t145.2313\t>10\t43039500\tF\t.\t.\t.\t3G>A";

        SingleEndAlignmentsToReadPairInfoMapper mapper = new SingleEndAlignmentsToReadPairInfoMapper();

        mapper.setFaix(new FaidxFileHelper("foo") {
            @Override
            public Short getKeyForChromName(String name) throws IOException {
                assertEquals("10", name);
                return (short) 9;
            }
        });
        mapper.setScorer(new ProbabilisticPairedAlignmentScorer());

        MockOutputCollector collector = new MockOutputCollector();
        Reporter reporter = Mockito.mock(Reporter.class);

        mapper.map(new LongWritable(1), new Text(inputLine), collector, reporter);

        int idx = 0;
        for (int i = 43039500; i <= 43049500; i = i + SVPipeline.RESOLUTION) {
            assertEquals((short) 9, collector.keys.get(idx).chromosome);
            assertEquals(i, collector.keys.get(idx).pos);

            assertEquals(10017, collector.values.get(idx).insertSize);
            assertEquals(-9.2103, collector.values.get(idx).pMappingCorrect, .0001);

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
//        for (int i = 81550200; i <= 81550500; i = i + SVPipeline.RESOLUTION) {
//            verify(collector).collect(new Text("10\t" + i),
//                    new DoubleWritable(1));
//        }
//
//        for (int i = 89072400; i <= 89072700; i = i + SVPipeline.RESOLUTION) {
//            verify(collector).collect(new Text("10\t" + i),
//                    new DoubleWritable(1));
//        }
//
//        verifyNoMoreInteractions(collector);
//

    }

    @Test
    public void testMapPairedEndWithSegmentalDuplications() throws Exception {

        String inputLine = "@ERR000545.10000001 EAS139_44:1:93:532:453\t@ERR000545.10000001 EAS139_44:1:93:532:453/1\tS\tCAAAAACCACTTGTACTCCAAAAGCTATTGAAGTTTAAGTTAAAATAAAAA\t<??>>?;<>=@?=?>8@<<9=98=:@>>>=:>>:6?7>9:?<46:;9;.:9\tR\t114\t0.00000\t>10\t43049466\tR\t.\t.\t.\t10A>C 22C>A 46-A\tSVP_READ\t@ERR000545.10000001 EAS139_44:1:93:532:453/2\tS\tTTATTGCACTTACCATGACTGTCTTCTGAAATGCATCTCAACCCTTGAATA\t;<8>>:$>=?@?>>:>=>:9=8<>1;8:<=>>9=:=7><>=;=;=>=72:>\tU\t34\t145.2313\t>10\t43039500\tF\t.\t.\t.\t3G>A";

        SingleEndAlignmentsToReadPairInfoMapper mapper = new SingleEndAlignmentsToReadPairInfoMapper();

        mapper.setFaix(new FaidxFileHelper("foo") {
            @Override
            public Short getKeyForChromName(String name) throws IOException {
                assertEquals("10", name);
                return (short) 9;
            }
        });
        mapper.setScorer(new ProbabilisticPairedAlignmentScorer());
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

        mapper.map(new LongWritable(1), new Text(inputLine), collector, reporter);

        assertEquals(0, collector.keys.size());

        inputLine = "@ERR000545.10000001 EAS139_44:1:93:532:453\t@ERR000545.10000001 EAS139_44:1:93:532:453/1\tS\tCAAAAACCACTTGTACTCCAAAAGCTATTGAAGTTTAAGTTAAAATAAAAA\t<??>>?;<>=@?=?>8@<<9=98=:@>>>=:>>:6?7>9:?<46:;9;.:9\tR\t114\t0.00000\t>10\t49466\tR\t.\t.\t.\t10A>C 22C>A 46-A\tSVP_READ\t@ERR000545.10000001 EAS139_44:1:93:532:453/2\tS\tTTATTGCACTTACCATGACTGTCTTCTGAAATGCATCTCAACCCTTGAATA\t;<8>>:$>=?@?>>:>=>:9=8<>1;8:<=>>9=:=7><>=;=;=>=72:>\tU\t34\t145.2313\t>10\t39500\tF\t.\t.\t.\t3G>A";
        collector = new MockOutputCollector();
        mapper.map(new LongWritable(1), new Text(inputLine), collector, reporter);
        assertEquals(101, collector.keys.size());
    }

        private static class MockOutputCollector implements OutputCollector<GenomicLocation, ReadPairInfo> {

        List<GenomicLocation> keys = new ArrayList<GenomicLocation>();
        List<ReadPairInfo> values = new ArrayList<ReadPairInfo>();

        public void collect(GenomicLocation key, ReadPairInfo value) throws IOException {
            keys.add(key);
            values.add(value);
        }

        public void reset() {
            keys.clear();
            values.clear();
        }
    }

}
