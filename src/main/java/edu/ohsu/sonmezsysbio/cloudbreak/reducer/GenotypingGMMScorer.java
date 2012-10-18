package edu.ohsu.sonmezsysbio.cloudbreak.reducer;

import com.google.common.primitives.Doubles;
import edu.ohsu.sonmezsysbio.cloudbreak.ReadGroupInfo;
import edu.ohsu.sonmezsysbio.cloudbreak.io.ReadPairInfo;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import java.util.*;

import static org.apache.commons.math3.stat.StatUtils.mean;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 10/10/12
 * Time: 9:17 PM
 */
public class GenotypingGMMScorer implements ReadPairInfoScorer {

    private static org.apache.log4j.Logger log = Logger.getLogger(GenotypingGMMScorer.class);
    public static final int MAX_COVERAGE = 200;

    { log.setLevel(Level.INFO); }

    private double[] pointLikelihoods(double[] y, double mu, double sigma) {
        double[] pointLikelihoods = new double[y.length];
        for (int i = 0; i < y.length; i++) {
            pointLikelihoods[i] = lognormal(y[i], mu, sigma);
        }
        return pointLikelihoods;
    }

    private double lognormal(double y, double mu, double sigma) {
        return -1.0 * Math.log(sigma + Math.sqrt(2.0 * Math.PI)) + -1.0 / 2.0 * Math.pow(((y - mu) / sigma), 2);
    }

    double likelihood(double[] y, double[] w, double[] mu, double sigma) {
        double[][] pointLikelihoods = new double[w.length][y.length];
        for (int i = 0; i < w.length; i++) {
            pointLikelihoods[i] = pointLikelihoods(y, mu[i], sigma);
        }

        double sumWeightedLikelihoods = 0;
        int i = 0;
        for (i = 0; i < y.length; i++) {
            for (int j = 0; j < w.length; j++) {
                sumWeightedLikelihoods += Math.exp(w[j]) * pointLikelihoods[j][i];
            }
        }
        return sumWeightedLikelihoods / y.length;
    }

    double logsumexp(double[] x) {
        double m = Doubles.max(x);
        double sumexp = 0;
        for (int i = 0; i < x.length; i++) {
            sumexp += Math.exp(x[i] - m);
        }
        double s = m + Math.log(sumexp);
        return s;
    }

    private double[][] gamma(double[] y, double[] w, double[] mu, double sigma) {
        double[] m1likelihoods = new double[y.length];
        double[] m2likelihoods = new double[y.length];
        double[][] gamma = new double[y.length][2];
        for (int i = 0; i < y.length; i++) {
            m1likelihoods[i] = w[0] + lognormal(y[i], mu[0], sigma);
            m2likelihoods[i] = w[1] + lognormal(y[i], mu[1], sigma);
            log.debug("m1likelihoods[" + i +"] " + m1likelihoods[i]);
            log.debug("m2likelihoods[" + i +"] " + m2likelihoods[i]);
            double total = logsumexp(new double[] {m1likelihoods[i], m2likelihoods[i]});
            log.debug("total likelihood[" + i + "]: "+ total);
            gamma[i][0] = m1likelihoods[i] - total;
            gamma[i][1] = m2likelihoods[i] - total;
        }
        return gamma;
    }

    private double[] cacluateN(double[][] gamma) {
        double[] ns = new double[2];
        double[] temp = new double[gamma.length];
        for (int i = 0; i < gamma.length; i++) {
            temp[i] = gamma[i][0];
        }
        ns[0] = logsumexp(temp);
        for (int i = 0; i < gamma.length; i++) {
            temp[i] = gamma[i][1];
        }
        ns[1] = logsumexp(temp);
        return ns;
    }

    private double[] updateW(double[] ns, double[] y) {
        double[] w = new double[2];
        w[0] = ns[0] - Math.log(y.length);
        w[1] = ns[1] - Math.log(y.length);
        return w;
    }

    double updateMu2(double[][] gamma, double[] y, double[] n) {
        double numerator = 0;
        double[] lse = new double[y.length];
        for (int i = 0; i < y.length; i++) {
            lse[i] = gamma[i][1] + Math.log(y[i]);
        }
        numerator = logsumexp(lse);
        log.debug("update mu2 n " +  numerator);
        return Math.exp(numerator - n[1]);
    }

    private static class EMUpdates {
        double[] w;
        double mu2;

        @Override
        public String toString() {
            return "EMUpdates{" +
                    "w=" + w[0] + ", " + w[1] +
                    ", mu2=" + mu2 +
                    '}';
        }
    }

    private EMUpdates emStep(double[] y, double[] w, double[] mu, double sigma) {
        double[][] gamma = gamma(y, w, mu, sigma);
        if (log.isDebugEnabled()) {
            for (int i = 0; i < gamma.length; i++) {
                for (int j = 0; j < gamma[i].length; j++) {
                    log.debug("gamma[" + i + "][" + j + "] = "  + gamma[i][j]);
                }
            }
        }
        double[] n = cacluateN(gamma);
        log.debug("n[0] " + n[0]);
        log.debug("n[1] " + n[1]);
        double[] wprime = updateW(n, y);
        log.debug("wprime[0] " + wprime[0]);
        log.debug("wprime[1] " + wprime[1]);
        double mu2prime = updateMu2(gamma, y, n);
        log.debug(mu2prime);
        EMUpdates updates = new EMUpdates();
        updates.w = wprime;
        updates.mu2 = mu2prime;
        return updates;
    }

    public double[] nnclean(double[] y, double sigma, int m) {
        if (m >= y.length) return new double[]{};
        List<Double> ysWithCloseNeighbors = new ArrayList<Double>();

        // todo: there's a much better way to do this: sort ys and only calculate distances of neigbors
        double[][] dist = new double[y.length][y.length];
        for (int i = 0; i < y.length; i++) {
            for (int j = 0; j < y.length; j++) {
                dist[i][j] = Math.abs(y[i] - y[j]);
            }
        }
        for (int i = 0; i < y.length; i++) {
            Arrays.sort(dist[i]);
            double distanceToMthNeighbor = dist[i][m];
            if (distanceToMthNeighbor < 5 * sigma) {
                ysWithCloseNeighbors.add(y[i]);
            }
        }

        double[] result = new double[ysWithCloseNeighbors.size()];
        for (int i = 0; i < ysWithCloseNeighbors.size(); i++) {
            result[i] = ysWithCloseNeighbors.get(i);
        }
        return result;
    }

    public double estimateW(double[] y, double[] initialW, double initialMu1, double sigma) {
        int maxIterations = 10;
        if (log.isDebugEnabled()) {
            log.debug("ys:");
            for (int i = 0; i < y.length; i++) {
                log.debug(y[i]);
            }
        }
        double[] yclean = nnclean(y, sigma, 2);
        if (log.isDebugEnabled()) {
            log.debug("ycleans:");
            for (int i = 0; i < yclean.length; i++) {
                log.debug(yclean[i]);
            }
        }

        if (yclean.length == 0) {
            log.debug("not enough ycleans, returning 1");
            return -1;
        }
        double[] initialMu = new double[]{initialMu1,mean(yclean)};

        int i = 1;
        double[] w = initialW;
        double[] mu = initialMu;
        double l = likelihood(yclean, w, mu, sigma);
        log.debug("initial likelihood: " + l);
        while(true) {
            EMUpdates updates = emStep(yclean, w, mu, sigma);
            if (log.isDebugEnabled()) {
                log.debug("updates: " + updates.toString());
            }
            w = updates.w;
            mu[1] = updates.mu2;
            double lprime = likelihood(yclean, w, mu, sigma);
            log.debug("new likelihood: " + l);
            i += 1;
            if (Math.abs(l - lprime) < 0.0001 || i > maxIterations) {
                break;
            }
            l = lprime;
        }
        if (Math.abs(mu[1] - mu[0]) < 2 * sigma) {
            log.debug("means too close, returning 1");
            return 1;
        }
        log.debug("returning " + Math.exp(w[0]));
        return Math.exp(w[0]);
    }

    public double reduceReadPairInfos(Iterator<ReadPairInfo> values, Map<Short, ReadGroupInfo> readGroupInfos) {
        List<Double> insertSizes = new ArrayList<Double>();
        double maxSD = 0;
        if (readGroupInfos.values().size() > 1) {
            throw new UnsupportedOperationException("GMM Reducer can't work with more than one read group right now");
        }
        boolean first = true;
        int targetIsize = 0;
        double bestMappingQuality = 0;
        while (values.hasNext()) {
            ReadPairInfo rpi = values.next();
            if (! first & (bestMappingQuality - rpi.pMappingCorrect > 6)) {
                break;
            }
            int insertSize = rpi.insertSize;
            short readGroupId = rpi.readGroupId;
            ReadGroupInfo readGroupInfo = readGroupInfos.get(readGroupId);
            insertSizes.add((double) (insertSize));
            if (readGroupInfo.isizeSD > maxSD) {
                maxSD = readGroupInfo.isizeSD;
            }
            if (first) {
                targetIsize = readGroupInfo.isize;
                bestMappingQuality = rpi.pMappingCorrect;
                first = false;
            }
        }
        if (insertSizes.size() >= MAX_COVERAGE) {
            insertSizes = insertSizes.subList(0, MAX_COVERAGE);
        }
        double[] initialW = new double[]{Math.log(.5),Math.log(.5)};
        double[] insertSizeArray = new double[insertSizes.size()];
        for (int i = 0; i < insertSizes.size(); i++) {
            insertSizeArray[i] = insertSizes.get(i);
        }

        return estimateW(insertSizeArray, initialW, targetIsize, maxSD);
    }
}
