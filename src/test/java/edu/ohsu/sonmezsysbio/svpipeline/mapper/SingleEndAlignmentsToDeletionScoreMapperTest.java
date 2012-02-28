package edu.ohsu.sonmezsysbio.svpipeline.mapper;

import edu.ohsu.sonmezsysbio.svpipeline.mapper.SingleEndAlignmentsToDeletionScoreMapper;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.OutputCollector;
import org.apache.hadoop.mapred.Reporter;
import org.junit.Test;
import org.mockito.Mockito;

import static org.junit.Assert.assertTrue;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.verifyNoMoreInteractions;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 6/5/11
 * Time: 9:10 PM
 */
public class SingleEndAlignmentsToDeletionScoreMapperTest {

    @Test
    public void testMapMatePairs() throws Exception {
        //String inputLine = "@SimSeq_100064\t@SimSeq_100064/1\tS\tGGAATGCTGAACCGAAGCCACACATATATTACATTGCTGATATTTCCATG\t22335345345444554445454566655554555335344444432221\tU\t159\t40\t>2\t41823857\tF\t.\t.\t.\t6A>G 13A>C 15T>A 19T>C 36A>G 42T>A 45A>T\tSEP\t@SimSeq_100064/2\tS\tATATCCGTATGTCAGATTTTATAATCTAGATATTTGATGGCTTTTTTTTA\t43555556665767677777777876D76777676566554665655543\tU\t79\t126\t>2\t41821123\tR\t.\t.\t.\t20T>A 22A>C 27C>T";
        String inputLine = "@SimSeq_100064\t@SimSeq_100064/1\tS\tGGAATGCTGAACCGAAGCCACACATATATTACATTGCTGATATTTCCATG\t22335345345444554445454566655554555335344444432221\tU\t159\t40\t>2\t41823857\tF\t.\t.\t.\t6A>G 13A>C 15T>A 19T>C 36A>G 42T>A 45A>T\tSVP_READ\t@SimSeq_100064/2\tS\tATATCCGTATGTCAGATTTTATAATCTAGATATTTGATGGCTTTTTTTTA\t43555556665767677777777876D76777676566554665655543\tU\t79\t126\t>2\t41821123\tR\t.\t.\t.\t20T>A 22A>C 27C>T";

        SingleEndAlignmentsToDeletionScoreMapper mapper = new SingleEndAlignmentsToDeletionScoreMapper();
        mapper.setTargetIsize(3000.0);
        mapper.setTargetIsizeSD(300.0);
        mapper.setMatePairs(true);

        OutputCollector<Text, DoubleWritable> collector = Mockito.mock(OutputCollector.class);
        Reporter reporter = Mockito.mock(Reporter.class);

        mapper.map(new LongWritable(1), new Text(inputLine), collector, reporter);

        for (int i = 41821100; i <= 41823900; i = i + 100) {
            verify(collector).collect(new Text("2\t" + i),
                    new DoubleWritable(-0.9998999999997488));
        }

        verifyNoMoreInteractions(collector);


    }

    @Test
    public void testMapPairedEnd() throws Exception {

        String inputLine = "@ERR000545.10000001 EAS139_44:1:93:532:453\t@ERR000545.10000001 EAS139_44:1:93:532:453/1\tS\tCAAAAACCACTTGTACTCCAAAAGCTATTGAAGTTTAAGTTAAAATAAAAA\t<??>>?;<>=@?=?>8@<<9=98=:@>>>=:>>:6?7>9:?<46:;9;.:9\tR\t114\t0\t>10\t43049466\tR\t.\t.\t.\t10A>C 22C>A 46-A\tSVP_READ\t@ERR000545.10000001 EAS139_44:1:93:532:453/2\tS\tTTATTGCACTTACCATGACTGTCTTCTGAAATGCATCTCAACCCTTGAATA\t;<8>>:$>=?@?>>:>=>:9=8<>1;8:<=>>9=:=7><>=;=;=>=72:>\tU\t34\t145\t>10\t43039500\tF\t.\t.\t.\t3G>A";

        SingleEndAlignmentsToDeletionScoreMapper mapper = new SingleEndAlignmentsToDeletionScoreMapper();
        mapper.setTargetIsize(200.0);
        mapper.setTargetIsizeSD(35.0);
        mapper.setMatePairs(false);

        OutputCollector<Text, DoubleWritable> collector = Mockito.mock(OutputCollector.class);
        Reporter reporter = Mockito.mock(Reporter.class);

        mapper.map(new LongWritable(1), new Text(inputLine), collector, reporter);

        for (int i = 43039500; i <= 43049500; i = i + 100) {
            verify(collector).collect(new Text("10\t" + i),
                    new DoubleWritable(9.999999999999969E-5));
        }

        verifyNoMoreInteractions(collector);

        String complexInputLine ="@ERR000545.10000241 EAS139_44:1:93:535:874\t@ERR000545.10000241 EAS139_44:1:93:535:874/1\tS\tTAGGTAATGTTTGGGAGGGAGTTGTTGGTTTTGTTGATTTATTATATCTTG\t??==>79>=?@@==>9>>=5>:?<>?<<<>??>;?=<=>?7?>;?8=4>>8\tR\t0\t60\t>10\t81550461\tR\t.\t.\t.\tSVP_ALIGNMENT\t@ERR000545.10000241 EAS139_44:1:93:535:874/1\tS\tTAGGTAATGTTTGGGAGGGAGTTGTTGGTTTTGTTGATTTATTATATCTTG\t??==>79>=?@@==>9>>=5>:?<>?<<<>??>;?=<=>?7?>;?8=4>>8\tR\t60\t0\t>10\t89072620\tR\t.\t.\t.\t48A>C 49T>C\tSVP_READ\t@ERR000545.10000241 EAS139_44:1:93:535:874/2\tS\tTATTATATTTTTTAATTTGACAGAGTAGTGCAGGCAATAATGAAATGGTAT\t<>==>A?@=?@A@>=??>>;>>;-<==;:<<<:>;==>>8?<::<>;80;>\tR\t0\t3\t>10\t89072431\tF\t.\t.\t.\tSVP_ALIGNMENT\t@ERR000545.10000241 EAS139_44:1:93:535:874/2\tS\tTATTATATTTTTTAATTTGACAGAGTAGTGCAGGCAATAATGAAATGGTAT\t<>==>A?@=?@A@>=??>>;>>;-<==;:<<<:>;==>>8?<::<>;80;>\tR\t0\t3\t>10\t81550273\tF\t.\t.\t.";

        collector = Mockito.mock(OutputCollector.class);
        reporter = Mockito.mock(Reporter.class);

        mapper.map(new LongWritable(1), new Text(complexInputLine), collector, reporter);

        for (int i = 81550200; i <= 81550500; i = i + 100) {
            verify(collector).collect(new Text("10\t" + i),
                    new DoubleWritable(-0.4988122675599614));
        }

        for (int i = 89072400; i <= 89072700; i = i + 100) {
            verify(collector).collect(new Text("10\t" + i),
                    new DoubleWritable(-4.988127663727278E-5));
        }

        verifyNoMoreInteractions(collector);


    }


    @Test
    public void testComputeDeletionScore() throws Exception {
        // bigger insert size: more likely deletion
        assertTrue(SingleEndAlignmentsToDeletionScoreMapper.computeDeletionScore(185, 117, 5114, 3000.0, 300.0) >
                   SingleEndAlignmentsToDeletionScoreMapper.computeDeletionScore(185, 117, 3114, 3000.0, 300.0));
        // higher quality: more likely deletion
        assertTrue(SingleEndAlignmentsToDeletionScoreMapper.computeDeletionScore(185, 117, 3114, 3000.0, 300.0) <
                   SingleEndAlignmentsToDeletionScoreMapper.computeDeletionScore(185, 20, 3114, 3000.0, 300.0));
    }

}
