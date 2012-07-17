package edu.ohsu.sonmezsysbio.cloudbreak;

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
        assertTrue(scorer.computeDeletionScore(5114, 3000.0, 300.0, scorer.probabilityMappingIsCorrect(185, 117)) >
                scorer.computeDeletionScore(3114, 3000.0, 300.0, scorer.probabilityMappingIsCorrect(185, 117)));
        // higher quality: more likely deletion
        assertTrue(scorer.computeDeletionScore(3114, 3000.0, 300.0, scorer.probabilityMappingIsCorrect(185, 117)) <
                scorer.computeDeletionScore(3114, 3000.0, 300.0, scorer.probabilityMappingIsCorrect(185, 20)));

        assertEquals(9.999999999999969E-5, scorer.computeDeletionScore(10000,200.0,35.0, scorer.probabilityMappingIsCorrect(0, 145)), 0.000001);

        assertEquals(-0.4988122675599614, scorer.computeDeletionScore(189,200.0,35.0, scorer.probabilityMappingIsCorrect(60,3)), 0.000001);

        assertEquals(-4.988127663727278E-5, scorer.computeDeletionScore(188,200.0,35.0, scorer.probabilityMappingIsCorrect(0,3)),0.000001);

        assertEquals(-0.9998999999997488, scorer.computeDeletionScore(2734,3000.0,300.0, scorer.probabilityMappingIsCorrect(40, 126)), 0.000001);

    }

}
