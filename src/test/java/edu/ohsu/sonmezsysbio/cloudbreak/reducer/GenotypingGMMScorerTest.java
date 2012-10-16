package edu.ohsu.sonmezsysbio.cloudbreak.reducer;

import edu.ohsu.sonmezsysbio.cloudbreak.ReadGroupInfo;
import edu.ohsu.sonmezsysbio.cloudbreak.io.ReadPairInfo;
import org.junit.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 10/12/12
 * Time: 10:11 AM
 */

/**
 * These tests are mostly direct ports from my R prototype
 */
public class GenotypingGMMScorerTest {

    @Test
    public void testLogsumexp() throws Exception {
        double[] x = new double[] {-7.018195, -309.17497};
        GenotypingGMMScorer scorer = new GenotypingGMMScorer();
        assertEquals(-7.018195, scorer.logsumexp(x), 0.00001);
    }

    @Test
    public void testNnclean() throws Exception {
        double[] y = new double[] {260.0736, 197.4272,   194.8618,  1217.8588,
                1228.2190,  1151.7017,  4511.5326, 19719.9700, 19707.2091, 16788.1891,
                22556.1203,  9909.8687, 13709.7259,  9219.4829, 13076.1122,  8140.4713};
        GenotypingGMMScorer scorer = new GenotypingGMMScorer();
        double[] nonnoise = new double[] {260.0736, 197.4272,   194.8618,  1217.8588,
                1228.2190,  1151.7017};
        assertArrayEquals(nonnoise, scorer.nnclean(y, 30, 2), 0.000001);
    }

    @Test
    public void testNncleanNotEnoughValues() throws Exception {
        double[] y = new double[] {260.0736, 1217.8588};
        GenotypingGMMScorer scorer = new GenotypingGMMScorer();
        double[] nonnoise = new double[] {};
        assertArrayEquals(nonnoise, scorer.nnclean(y, 30, 2), 0.000001);
    }

    @Test
    public void testEstimateW() throws Exception {
        double[] y = new double[] {260.0736, 197.4272,   194.8618,  1217.8588,
                1228.2190,  1151.7017,  4511.5326, 19719.9700, 19707.2091, 16788.1891,
                22556.1203,  9909.8687, 13709.7259,  9219.4829, 13076.1122,  8140.4713};

        double sigma = 30;
        double[] initialW = new double[] {Math.log(.5),Math.log(.5)};
        double[] initialMu = new double[] {200, 1000};

        GenotypingGMMScorer scorer = new GenotypingGMMScorer();
        assertEquals(.5, Math.exp(scorer.estimateW(y, initialW, initialMu, sigma)), 0.0001);
    }

    public void testEstimateWNoValues() throws Exception {
        double[] y = new double[] {};

        double sigma = 30;
        double[] initialW = new double[] {Math.log(.5),Math.log(.5)};
        double[] initialMu = new double[] {200, 1000};

        GenotypingGMMScorer scorer = new GenotypingGMMScorer();
        assertEquals(1, Math.exp(scorer.estimateW(y, initialW, initialMu, sigma)), 0.0001);
    }

    public void testReduce_2_91700() throws Exception {
        double[] y = new double[] {152,
                216,
                194,
                169,
                202,
                237,
                178,
                202,
                208,
                210,
                247,
                227,
                182,
                207,
                191,
                147,
                251,
                236,
                209,
                200,
                204,
                184,
                158,
                152,
                255,
                227,
                185,
                183,
                204,
                226,
                276,
                207,
                223
        };

        ReadGroupInfo readGroupInfo = new ReadGroupInfo();
        readGroupInfo.isize = 200;
        readGroupInfo.isizeSD = 30;

        List<ReadPairInfo> rpis = new ArrayList<ReadPairInfo>();
        for (int i = 0; i < y.length; i++) {
            double isize = y[i];
            ReadPairInfo rpi = new ReadPairInfo();
            rpi.insertSize = (int) isize;
            rpi.readGroupId = 0;
            rpis.add(rpi);
        }

        GenotypingGMMScorer scorer = new GenotypingGMMScorer();

        Map<Short, ReadGroupInfo>  rgis = new HashMap<Short, ReadGroupInfo>();
        rgis.put((short) 0, readGroupInfo);
        double score = scorer.reduceReadPairInfos(rpis.iterator(), rgis);
        assertEquals(1, score, 0.00001);
    }
}
