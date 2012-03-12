package edu.ohsu.sonmezsysbio.svpipeline.mapper;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 2/28/12
 * Time: 4:18 PM
 */
public class SingleEndAlignmentsToBedSpansMapperTest {
    @Test
    public void testParseRegion() {
        SingleEndAlignmentsToBedSpansMapper mapper = new SingleEndAlignmentsToBedSpansMapper();
        mapper.parseRegion("1:54817520-54915143");
        assertEquals("1", mapper.getChromosome());
        assertEquals(54817520, mapper.getRegionStart());
        assertEquals(54915143, mapper.getRegionEnd());
    }
}
