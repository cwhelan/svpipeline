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
public class GenotypingGMMScorer {

    private static org.apache.log4j.Logger log = Logger.getLogger(GenotypingGMMScorer.class);
    public static final int MAX_COVERAGE = 200;
    public static final int MAX_LOG_MAPQ_DIFF = 5;

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
        double[][] gamma = new double[y.length][w.length];
        double[][] likelihoods = new double[y.length][w.length];
        for (int i = 0; i < y.length; i++) {
            for (int j = 0; j < w.length; j++) {
                likelihoods[i][j] = w[j] + lognormal(y[i], mu[j], sigma);
                log.debug("likelihoods[" + i +"][" + j + "] " + likelihoods[i][j]);
            }
            double total = logsumexp(likelihoods[i]);
            log.debug("total likelihood[" + i + "]: "+ total);
            for (int j = 0; j < w.length; j++) {
                gamma[i][j] = likelihoods[i][j] - total;
            }
        }
        return gamma;
    }

    private double[] cacluateN(double[][] gamma) {
        double[] ns = new double[gamma[0].length];
        for (int j = 0; j < gamma[0].length; j++) {
            double[] temp = new double[gamma.length];
            for (int i = 0; i < gamma.length; i++) {
                temp[i] = gamma[i][j];
            }
            ns[j] = logsumexp(temp);
        }
        return ns;
    }

    private double[] updateW(double[] ns, double[] y) {
        double[] w = new double[ns.length];
        for (int j = 0; j < ns.length; j++) {
            w[j] = ns[j] - Math.log(y.length);
        }
        return w;
    }

    double updateMuForComponent(double[][] gamma, double[] y, double[] n, int component) {
        double numerator = 0;
        double[] lse = new double[y.length];
        for (int i = 0; i < y.length; i++) {
            lse[i] = gamma[i][component] + Math.log(y[i]);
        }
        numerator = logsumexp(lse);
        log.debug("update mu[" + component + "] n " +  numerator);
        return Math.exp(numerator - n[component]);
    }

    private static class EMUpdates {
        double[] w;
        double[] mu;

        @Override
        public String toString() {
            return "EMUpdates{" +
                    "w=" + w +
                    ", mu=" + mu +
                    '}';
        }
    }

    private EMUpdates emStep(double[] y, double[] w, double[] mu, double sigma, int[] freeMus) {
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
        EMUpdates updates = new EMUpdates();
        updates.mu = Arrays.copyOf(mu, mu.length);

        double[] wprime = updateW(n, y);
        updates.w = wprime;
        for (int j = 0; j < wprime.length; j++) {
            log.debug("wprime[" + j + "] " + wprime[j]);
        }

        for (int k = 0; k < freeMus.length; k++) {
            int j = freeMus[k];
            double mujprime = updateMuForComponent(gamma, y, n, j);
            log.debug("mu[" + j + "] prime: " + mujprime);
            updates.mu[j] = mujprime;
        }

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

    public GMMScorerResults estimateW(double[] y, double[] initialW, double initialMu1, double sigma) {
        GMMScorerResults results = new GMMScorerResults();
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
            results.w0 = -1;
            return results;
        }
        results.nodelOneComponentLikelihood = likelihood(yclean, new double[]{Math.log(1)}, new double[]{initialMu1}, sigma);
        double[] initialMu = new double[]{initialMu1,mean(yclean)};

        int i = 1;
        double[] w = initialW;
        double[] mu = initialMu;
        double l = likelihood(yclean, w, mu, sigma);
        log.debug("initial likelihood: " + l);
        while(true) {
            EMUpdates updates = emStep(yclean, w, mu, sigma, new int[] {1});
            if (log.isDebugEnabled()) {
                log.debug("updates: " + updates.toString());
            }
            w = updates.w;
            mu[1] = updates.mu[1];
            double lprime = likelihood(yclean, w, mu, sigma);
            log.debug("new likelihood: " + l);
            i += 1;
            if (Math.abs(l - lprime) < 0.0001 || i > maxIterations) {
                break;
            }
            l = lprime;
        }
        results.twoComponentLikelihood = l;
        results.likelihoodRatio = l - results.nodelOneComponentLikelihood;
        results.mu2 = mu[1];
        if (Math.abs(mu[1] - mu[0]) < 2 * sigma) {
            log.debug("means too close, returning 1");
            results.w0 = 1;
            return results;
        }
        log.debug("returning " + Math.exp(w[0]));
        results.w0 = Math.exp(w[0]);
        return results;
    }

    public GMMScorerResults reduceReadPairInfos(Iterator<ReadPairInfo> values, Map<Short, ReadGroupInfo> readGroupInfos) {
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
            if (! first & (bestMappingQuality - rpi.pMappingCorrect > MAX_LOG_MAPQ_DIFF)) {
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
