package edu.ohsu.sonmezsysbio.svpipeline;

import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 3/15/12
 * Time: 12:54 PM
 */
public class ProbabilisticPairedAlignmentScorerTest {
    @Test
    public void testComputeDeletionScore() throws Exception {
        ProbabilisticPairedAlignmentScorer scorer = new ProbabilisticPairedAlignmentScorer();
        // bigger insert size: more likely deletion
        assertTrue(scorer.computeDeletionScore(185, 117, 5114, 3000.0, 300.0) >
                scorer.computeDeletionScore(185, 117, 3114, 3000.0, 300.0));
        // higher quality: more likely deletion
        assertTrue(scorer.computeDeletionScore(185, 117, 3114, 3000.0, 300.0) >
                scorer.computeDeletionScore(185, 0, 3114, 3000.0, 300.0));

        assertTrue(scorer.computeDeletionScore(185, 117, 5114, 200.0, 30.0) > 0);
    }

}
