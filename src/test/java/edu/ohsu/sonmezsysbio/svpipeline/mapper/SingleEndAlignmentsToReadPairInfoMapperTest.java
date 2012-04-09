package edu.ohsu.sonmezsysbio.svpipeline.mapper;

import edu.ohsu.sonmezsysbio.svpipeline.ProbabilisticPairedAlignmentScorer;
import edu.ohsu.sonmezsysbio.svpipeline.SVPipeline;
import edu.ohsu.sonmezsysbio.svpipeline.io.GenomicLocation;
import edu.ohsu.sonmezsysbio.svpipeline.io.ReadPairInfo;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.OutputCollector;
import org.apache.hadoop.mapred.Reporter;
import org.hamcrest.BaseMatcher;
import org.hamcrest.Description;
import org.hamcrest.Matcher;
import org.junit.Test;
import org.mockito.Mockito;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertEquals;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.verifyNoMoreInteractions;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/6/12
 * Time: 4:58 PM
 */
public class SingleEndAlignmentsToReadPairInfoMapperTest {
    @Test
    public void testReadFaidx() throws IOException {
        String faidx = "1       249250621       52      60      61\n" +
                "2       243199373       253404903       60      61\n" +
                "3       198022430       500657651       60      61\n" +
                "4       191154276       701980507       60      61\n" +
                "5       180915260       896320740       60      61\n" +
                "6       171115067       1080251307      60      61\n" +
                "7       159138663       1254218344      60      61\n" +
                "8       146364022       1416009371      60      61\n" +
                "9       141213431       1564812846      60      61\n" +
                "10      135534747       1708379889      60      61\n" +
                "11      135006516       1846173603      60      61\n" +
                "12      133851895       1983430282      60      61\n" +
                "13      115169878       2119513096      60      61\n" +
                "14      107349540       2236602526      60      61\n" +
                "15      102531392       2345741279      60      61\n" +
                "16      90354753        2449981581      60      61\n" +
                "17      81195210        2541842300      60      61\n" +
                "18      78077248        2624390817      60      61\n" +
                "19      59128983        2703769406      60      61\n" +
                "20      63025520        2763883926      60      61\n" +
                "21      48129895        2827959925      60      61\n" +
                "22      51304566        2876892038      60      61\n" +
                "X       155270560       2929051733      60      61\n" +
                "Y       59373566        3086910193      60      61\n" +
                "MT      16569   3147273397      70      71\n" +
                "GL000234.1      40531   3147923585      60      61\n";
        SingleEndAlignmentsToReadPairInfoMapper mapper = new SingleEndAlignmentsToReadPairInfoMapper();
        Map chromTable = mapper.readFaidx(new BufferedReader(new StringReader(faidx)));
        assertEquals((short) 1, chromTable.get("2"));
        assertEquals((short) 23, chromTable.get("Y"));
        assertEquals((short) 25, chromTable.get("GL000234.1"));
    }

    @Test
    public void testMapPairedEnd() throws Exception {

        String inputLine = "@ERR000545.10000001 EAS139_44:1:93:532:453\t@ERR000545.10000001 EAS139_44:1:93:532:453/1\tS\tCAAAAACCACTTGTACTCCAAAAGCTATTGAAGTTTAAGTTAAAATAAAAA\t<??>>?;<>=@?=?>8@<<9=98=:@>>>=:>>:6?7>9:?<46:;9;.:9\tR\t114\t0\t>10\t43049466\tR\t.\t.\t.\t10A>C 22C>A 46-A\tSVP_READ\t@ERR000545.10000001 EAS139_44:1:93:532:453/2\tS\tTTATTGCACTTACCATGACTGTCTTCTGAAATGCATCTCAACCCTTGAATA\t;<8>>:$>=?@?>>:>=>:9=8<>1;8:<=>>9=:=7><>=;=;=>=72:>\tU\t34\t145\t>10\t43039500\tF\t.\t.\t.\t3G>A";

        SingleEndAlignmentsToReadPairInfoMapper mapper = new SingleEndAlignmentsToReadPairInfoMapper();
        HashMap<String, Short> chromosomeKeys = new HashMap<String, Short>();
        chromosomeKeys.put("10", (short) 9);
        mapper.setChromosomeKeys(chromosomeKeys);
        mapper.setScorer(new ProbabilisticPairedAlignmentScorer());

        MockOutputCollector collector = new MockOutputCollector();
        Reporter reporter = Mockito.mock(Reporter.class);

        mapper.map(new LongWritable(1), new Text(inputLine), collector, reporter);

        int idx = 0;
        for (int i = 43039500; i <= 43049500; i = i + SVPipeline.RESOLUTION) {
            assertEquals((short) 9, collector.keys.get(idx).chromosome);
            assertEquals(i, collector.keys.get(idx).pos);

            assertEquals(10051, collector.values.get(idx).insertSize);
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
