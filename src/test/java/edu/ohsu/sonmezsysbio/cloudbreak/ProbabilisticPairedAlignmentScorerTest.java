package edu.ohsu.sonmezsysbio.cloudbreak;

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
        assertTrue(scorer.computeDeletionScore(5114, 3000.0, 300.0, scorer.probabilityMappingIsCorrect(185, 117)) >
                scorer.computeDeletionScore(3114, 3000.0, 300.0, scorer.probabilityMappingIsCorrect(185, 117)));
        // higher quality: less likely deletion since insert size close to mean
        assertTrue(scorer.computeDeletionScore(3114, 3000.0, 300.0, scorer.probabilityMappingIsCorrect(185, 117)) <
                scorer.computeDeletionScore(3114, 3000.0, 300.0, scorer.probabilityMappingIsCorrect(185, 0)));

        // higher quality: more likely deletion since insert size far from mean
        assertTrue(scorer.computeDeletionScore(7000, 3000.0, 300.0, scorer.probabilityMappingIsCorrect(185, 117)) >
                scorer.computeDeletionScore(7000, 3000.0, 300.0, scorer.probabilityMappingIsCorrect(185, 0)));

        // bigger insert size: more likely deletion
        assertTrue(scorer.computeDeletionScore(5114, 3000.0, 300.0, scorer.probabilityMappingIsCorrect(185, 117)) >
                scorer.computeDeletionScore(2000, 3000.0, 300.0, scorer.probabilityMappingIsCorrect(185, 117)));

        assertTrue(scorer.computeDeletionScore(5114, 200.0, 30.0, scorer.probabilityMappingIsCorrect(185, 117)) > 0);

    }

    @Test
    public void testProbabilityMappingIsCorrect() {
        ProbabilisticPairedAlignmentScorer scorer = new ProbabilisticPairedAlignmentScorer();
        assertTrue(scorer.probabilityMappingIsCorrect(185,117) > scorer.probabilityMappingIsCorrect(185,0));
    }
}
