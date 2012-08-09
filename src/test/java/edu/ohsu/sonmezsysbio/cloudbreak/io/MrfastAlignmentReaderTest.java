package edu.ohsu.sonmezsysbio.cloudbreak.io;

import edu.ohsu.sonmezsysbio.cloudbreak.AlignmentRecord;
import edu.ohsu.sonmezsysbio.cloudbreak.MrfastAlignmentRecord;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 8/8/12
 * Time: 4:39 PM
 */
public class MrfastAlignmentReaderTest {

    @Test
    public void testComputeAlignmentScores() throws Exception {
        MrfastAlignmentRecord record1 = new MrfastAlignmentRecord();
        record1.setMismatches(0);

        List<AlignmentRecord> records = new ArrayList<AlignmentRecord>();
        records.add(record1);

        Map<AlignmentRecord, Double> alignmentScores = MrfastAlignmentReader.computeAlignmentScores(records);
        assertEquals(0, alignmentScores.get(record1), 0.000001);
    }
}
