package edu.ohsu.sonmezsysbio.cloudbreak.reducer;

import edu.ohsu.sonmezsysbio.cloudbreak.ReadGroupInfo;
import edu.ohsu.sonmezsysbio.svpipeline.io.GenomicLocation;
import edu.ohsu.sonmezsysbio.cloudbreak.io.ReadPairInfo;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.mapred.OutputCollector;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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

        ReadPairInfo readPairInfo1 = new ReadPairInfo(3000, -0.69, (short) 0);
        ReadPairInfo readPairInfo2 = new ReadPairInfo(3000, -9.2103, (short) 0);

        MockOutputCollector outputCollector = new MockOutputCollector();

        Map<Short, ReadGroupInfo> readGroupInfos = new HashMap<Short, ReadGroupInfo>();
        ReadGroupInfo rg1 = new ReadGroupInfo();
        rg1.isize = 200;
        rg1.isizeSD = 30;
        rg1.matePair = false;
        readGroupInfos.put((short) 0, rg1);

        IncrementalDelBeliefUpdateReadPairInfoReducer reducer = new IncrementalDelBeliefUpdateReadPairInfoReducer();
        reducer.setReadGroupInfos(readGroupInfos);

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