package edu.ohsu.sonmezsysbio.cloudbreak;

import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/23/11
 * Time: 10:50 AM
 */
public class SAMRecordTest {
    @Test
    public void testFlags() throws Exception {
        SAMRecord samRecord = new SAMRecord();
        samRecord.setFlag(403);
        assertTrue(samRecord.isMapped());
        assertTrue(samRecord.isMateMapped());
        assertTrue(samRecord.isPairMapped());

        assertTrue(! samRecord.isReverseComplemented());

        samRecord.setFlag(355);
        assertTrue(samRecord.isReverseComplemented());
    }

    @Test
    public void testPairPosteriorTag() throws Exception {
        SAMRecord samRecord = new SAMRecord();
        samRecord.addTag("PQ:i:340");
        assertEquals(340, samRecord.getPairPosterior());
    }

    @Test
    public void testOrientationFlags() throws Exception {
        SAMRecord r1 = new SAMRecord();
        r1.setFlag(0x99);

        SAMRecord r2 = new SAMRecord();
        r2.setFlag(0x147);

        assertTrue(r1.isForward());
        assertTrue(! r2.isForward());
    }

}
