package edu.ohsu.sonmezsysbio.svpipeline.reducer;

import edu.ohsu.sonmezsysbio.svpipeline.io.GenomicLocation;
import edu.ohsu.sonmezsysbio.svpipeline.io.ReadPairInfo;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.mapred.OutputCollector;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/8/12
 * Time: 9:53 AM
 */
public class IncrementalDelBeliefUpdateReadPairInfoReducerTest {
    @Test
    public void testReduce() throws Exception {
        GenomicLocation genomicLocation = new GenomicLocation((short) 1,10000);

        ReadPairInfo readPairInfo1 = new ReadPairInfo(3000, -0.69);
        ReadPairInfo readPairInfo2 = new ReadPairInfo(3000, -9.2103);

        MockOutputCollector outputCollector = new MockOutputCollector();

        IncrementalDelBeliefUpdateReadPairInfoReducer reducer = new IncrementalDelBeliefUpdateReadPairInfoReducer();
        reducer.setTargetIsize(200);
        reducer.setTargetIsizeSD(30);

        List<ReadPairInfo> readPairInfos = new ArrayList<ReadPairInfo>();
        readPairInfos.add(readPairInfo1);

        reducer.reduce(genomicLocation, readPairInfos.iterator(), outputCollector, null);

        assertEquals(1, outputCollector.keys.size());
        assertEquals(1, outputCollector.values.size());

        assertEquals((short) 1, outputCollector.keys.get(0).chromosome);
        assertEquals(10000, outputCollector.keys.get(0).pos);
        assertEquals(1.550446e-06, outputCollector.values.get(0).get(), 0.01);

        readPairInfos.add(readPairInfo2);
        outputCollector.reset();

        reducer.reduce(genomicLocation, readPairInfos.iterator(), outputCollector, null);

        assertEquals(1, outputCollector.keys.size());
        assertEquals(1, outputCollector.values.size());

        assertEquals((short) 1, outputCollector.keys.get(0).chromosome);
        assertEquals(10000, outputCollector.keys.get(0).pos);
        assertEquals(0.0065, outputCollector.values.get(0).get(), 0.1);
    }

    private static class MockOutputCollector implements OutputCollector<GenomicLocation, DoubleWritable> {

        List<GenomicLocation> keys = new ArrayList<GenomicLocation>();
        List<DoubleWritable> values = new ArrayList<DoubleWritable>();

        public void collect(GenomicLocation key, DoubleWritable value) throws IOException {
            keys.add(key);
            values.add(value);
        }

        public void reset() {
            keys.clear();
            values.clear();
        }
    }
}
