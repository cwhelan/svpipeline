package edu.ohsu.sonmezsysbio.svpipeline;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 3/12/12
 * Time: 1:56 PM
 */
public class NovoalignNativeRecordTest {
    @Test
    public void testParseRecord() {
        String record = "@ERR000545.10000001 EAS139_44:1:93:532:453/1\tS\tCAAAAACCACTTGTACTCCAAAAGCTATTGAAGTTTAAGTTAAAATAAAAA\t<??>>?;<>=@?=?>8@<<9=98=:@>>>=:>>:6?7>9:?<46:;9;.:9\tR\t114\t0\t>10\t43049466\tR\t.\t.\t.\t10A>C 22C>A 46-A";
        String[] fields = record.split("\t");
        NovoalignNativeRecord novoalignNativeRecord = NovoalignNativeRecord.parseRecord(fields);
        assertEquals(0, novoalignNativeRecord.getPosteriorProb());
        
        record = "@ERR000545.10000001 EAS139_44:1:93:532:453/2\tS\tTTATTGCACTTACCATGACTGTCTTCTGAAATGCATCTCAACCCTTGAATA\t;<8>>:$>=?@?>>:>=>:9=8<>1;8:<=>>9=:=7><>=;=;=>=72:>\tU\t34\t145\t>10\t43039500\tF\t.\t.\t.\t3G>A";
        fields = record.split("\t");
        novoalignNativeRecord = NovoalignNativeRecord.parseRecord(fields);
        assertEquals(145, novoalignNativeRecord.getPosteriorProb());
    }
}
