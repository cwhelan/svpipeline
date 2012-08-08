package edu.ohsu.sonmezsysbio.cloudbreak.mapper;

import edu.ohsu.sonmezsysbio.cloudbreak.io.NovoalignAlignmentReader;
import org.junit.Before;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 2/28/12
 * Time: 4:18 PM
 */
public class SingleEndAlignmentsToBedSpansMapperTest {

    SingleEndAlignmentsToBedSpansMapper mapper;

    @Before
    public void setup() {
        mapper = new SingleEndAlignmentsToBedSpansMapper();
        mapper.setAlignmentReader(new NovoalignAlignmentReader());
    }

    @Test
    public void testParseRegion() {
        mapper.parseRegion("1:54817520-54915143");
        assertEquals("1", mapper.getChromosome());
        assertEquals(54817520, mapper.getRegionStart());
        assertEquals(54915143, mapper.getRegionEnd());
    }
}
