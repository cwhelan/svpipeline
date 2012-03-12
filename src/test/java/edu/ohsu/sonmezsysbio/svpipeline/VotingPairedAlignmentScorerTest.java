package edu.ohsu.sonmezsysbio.svpipeline;

import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 3/12/12
 * Time: 1:32 PM
 */
public class VotingPairedAlignmentScorerTest {

    @Test
    public void testComputeDeletionScore() throws Exception {
        VotingPairedAlignmentScorer scorer = new VotingPairedAlignmentScorer();
        // bigger insert size: more likely deletion
        assertTrue(scorer.computeDeletionScore(185, 117, 5114, 3000.0, 300.0) >
                scorer.computeDeletionScore(185, 117, 3114, 3000.0, 300.0));
        // higher quality: more likely deletion
        assertTrue(scorer.computeDeletionScore(185, 117, 3114, 3000.0, 300.0) <
                scorer.computeDeletionScore(185, 20, 3114, 3000.0, 300.0));

        assertEquals(9.999999999999969E-5, scorer.computeDeletionScore(0,145,10000,200.0,35.0), 0.000001);

        assertEquals(-0.4988122675599614, scorer.computeDeletionScore(60,3,189,200.0,35.0), 0.000001);

        assertEquals(-4.988127663727278E-5, scorer.computeDeletionScore(0,3,188,200.0,35.0), 0.000001);

        assertEquals(-0.9998999999997488, scorer.computeDeletionScore(40,126,2734,3000.0,300.0), 0.000001);

    }

}
