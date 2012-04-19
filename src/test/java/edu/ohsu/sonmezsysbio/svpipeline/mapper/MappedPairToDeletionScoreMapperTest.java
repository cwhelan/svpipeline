package edu.ohsu.sonmezsysbio.svpipeline.mapper;

import edu.ohsu.sonmezsysbio.svpipeline.SAMRecord;
import edu.ohsu.sonmezsysbio.svpipeline.mapper.MappedPairToDeletionScoreMapper;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.OutputCollector;
import org.apache.hadoop.mapred.Reporter;
import org.junit.Test;
import org.mockito.Mockito;

import static org.junit.Assert.*;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.verifyNoMoreInteractions;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/23/11
 * Time: 11:28 AM
 */
public class MappedPairToDeletionScoreMapperTest {
    @Test
    public void testParseSamRecord() throws Exception {
        String samLine = "SimSeq_1a\t403\t2\t107742996\t0\t1S49M\t=\t107745559\t2515\t"+
                "GAGAATCCCACGTAGCTAAGATAGTAGGCGTGTACCCCCATACTCAGCTA\t+,+#--,,,-,,.-,,...,..--,.---,.,..,,,,,---,-,.-,--\t" +
                "PG:Z:novoalign\tRG:Z:RG1\tPU:Z:ILLUMINA\tLB:Z:LC_1\tZO:Z:-+ AS:i:170\tUQ:i:170\tNM:i:10\tMD:Z:2C0T5A6G0G3T0A0C8G9C6\tPQ:i:340\tSM:i:0\tAM:i:0\tZS:Z:R\tZN:i:2";
        SAMRecord samRecord = SAMRecord.parseSamRecord(samLine.split("\t"));
        assertTrue(samRecord.isPairMapped());
        assertFalse(samRecord.isReverseComplemented());
        assertEquals("2", samRecord.getReferenceName());
        assertEquals("=", samRecord.getPairReferenceName());
        assertFalse(samRecord.isInterChromosomal());
        assertEquals(2515, samRecord.getInsertSize());
        assertEquals(107742996, samRecord.getPosition());
        assertEquals(340, samRecord.getPairPosterior());
    }

    @Test
    public void testMap() throws Exception {
        String inputLine = "SimSeq_253e3e-0\tSimSeq_253e3e\t99\t2\t160043044\t150\t48M2S\t=\t160039861\t-3134\tCTTCACAGGCTGAGGGATGAGAAACACTTGGGCCTGGGAGGTCGAGGTCG\t78:9;;<:<<;==<=<=:<?5==.:<<;<;<;;<;;<;:;;;<;;;7666\tPG:Z:novoalign\tRG:Z:RGID\tPU:Z:ILLUMINA\tLB:Z:SIMTEST\tAS:i:163\tUQ:i:163\tNM:i:5\tMD:Z:3G1G17T6A11G5\tPQ:i:342\tSM:i:150\tAM:i:3\tSEP\tSimSeq_253e3e\t147\t2\t160039861\t150\t50M\t=\t160043044\t3134\tGGTTCGTCAGGTTGGGGGAGAAGCTGACAGTTGAAGATATTGCATCTCTG\t66899::;;;:6H;;<;8=:<=;<<=<<>;<<;==;E=<=<;;:;:9885\tPG:Z:novoalign\tRG:Z:RGID\tPU:Z:ILLUMINA\tLB:Z:SIMTEST\tAS:i:179\tUQ:i:179\tNM:i:6\tMD:Z:2G3G8A22C1A5A3\tPQ:i:342\tSM:i:3\tAM:i:3";

        MappedPairToDeletionScoreMapper mapper = new MappedPairToDeletionScoreMapper();
        mapper.setTargetIsize(3000.0);
        mapper.setTargetIsizeSD(300.0);

        OutputCollector<Text, DoubleWritable> collector = Mockito.mock(OutputCollector.class);
        Reporter reporter = Mockito.mock(Reporter.class);

        mapper.map(new LongWritable(1), new Text(inputLine), collector, reporter);

        verify(collector).collect(new Text("2\t160039500"),
                new DoubleWritable(0.08226443867766886));
        verify(collector).collect(new Text("2\t160040000"),
                new DoubleWritable(0.08226443867766886));
        verify(collector).collect(new Text("2\t160040500"),
                new DoubleWritable(0.08226443867766886));
        verify(collector).collect(new Text("2\t160041000"),
                new DoubleWritable(0.08226443867766886));
        verify(collector).collect(new Text("2\t160041500"),
                new DoubleWritable(0.08226443867766886));
        verify(collector).collect(new Text("2\t160042000"),
                new DoubleWritable(0.08226443867766886));
        verify(collector).collect(new Text("2\t160042500"),
                new DoubleWritable(0.08226443867766886));
        verify(collector).collect(new Text("2\t160043000"),
                new DoubleWritable(0.08226443867766886));
        verify(collector).collect(new Text("2\t160043500"),
                new DoubleWritable(0.08226443867766886));
        verifyNoMoreInteractions(collector);


    }

    @Test
    public void testComputeDeletionScore() throws Exception {
        // bigger insert size: more likely deletion
        assertTrue(MappedPairToDeletionScoreMapper.computeDeletionScore(185.2, 117.1, 5114, 3000.0, 300.0, true) >
                   MappedPairToDeletionScoreMapper.computeDeletionScore(185.2, 117.1, 3114, 3000.0, 300.0, true));
        // higher quality: more likely deletion
        assertTrue(MappedPairToDeletionScoreMapper.computeDeletionScore(185.2, 117.1, 3114, 3000.0, 300.0, true) >
                   MappedPairToDeletionScoreMapper.computeDeletionScore(185.2, 20.1231, 3114, 3000.0, 300.0, true));
    }
}
