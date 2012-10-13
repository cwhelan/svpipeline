package edu.ohsu.sonmezsysbio.cloudbreak.reducer;

import com.google.common.primitives.Doubles;
import edu.ohsu.sonmezsysbio.cloudbreak.ReadGroupInfo;
import edu.ohsu.sonmezsysbio.cloudbreak.io.ReadPairInfo;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 10/10/12
 * Time: 9:17 PM
 */
public class GenotypingGMMScorer implements ReadPairInfoScorer {

    private double[] pointLikelihoods(double[] y, NormalDistribution d) {
        double[] pointLikelihoods = new double[y.length];
        for (int i = 0; i < y.length; i++) {
            pointLikelihoods[i] = Math.log(d.density(y[i]));
        }
        return pointLikelihoods;
    }

    private double likelihood(double[] y, double[] w, double[] mu, double sigma) {
        NormalDistribution d1 = new NormalDistribution(mu[0], sigma);
        NormalDistribution d2 = new NormalDistribution(mu[1], sigma);
        double[] d1PointLikelihoods = pointLikelihoods(y, d1);
        double[] d2PointLikelihoods = pointLikelihoods(y, d2);
        double sumWeightedLikelihoods = 0;
        int i = 0;
        for (i = 0; i < y.length; i++) {
            sumWeightedLikelihoods += Math.exp(w[0]) * d1PointLikelihoods[i] + Math.exp(w[1]) * d2PointLikelihoods[i];
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
        NormalDistribution d1 = new NormalDistribution(mu[0], sigma);
        NormalDistribution d2 = new NormalDistribution(mu[1], sigma);
        double[] m1likelihoods = new double[y.length];
        double[] m2likelihoods = new double[y.length];
        double[][] gamma = new double[y.length][2];
        for (int i = 0; i < y.length; i++) {
            m1likelihoods[i] = w[0] + Math.log(d1.density(y[i]));
            m2likelihoods[i] = w[1] + Math.log(d2.density(y[i]));
            double total = logsumexp(new double[] {m1likelihoods[i], m2likelihoods[i]});
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

    double updateMu2(double[][] gamma, double[] y) {
        double n = 0;
        double d = 0;
        for (int i = 0; i < y.length; i++) {
            n += Math.exp(gamma[i][1] + Math.log(y[i]));
            d += Math.exp(gamma[i][1]);
        }
        return n / d;
    }

    private static class EMUpdates {
        double[] w;
        double mu2;
    }

    private EMUpdates emStep(double[] y, double[] w, double[] mu, double sigma) {
        double[][] gamma = gamma(y, w, mu, sigma);
        double[] n = cacluateN(gamma);
        double[] wprime = updateW(n, y);
        double mu2prime = updateMu2(gamma, y);
        EMUpdates updates = new EMUpdates();
        updates.w = wprime;
        updates.mu2 = mu2prime;
        return updates;
    }

    public double[] nnclean(double[] y, double sigma, int m) {
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

    public double estimateW(double[] y, double[] initialW, double[] initialMu, double sigma) {
        int maxIterations = 10;
        double[] yclean = nnclean(y, sigma, 2);
        if (yclean.length == 0) {
            return 1;
        }
        int i = 1;
        double[] w = initialW;
        double[] mu = initialMu;
        double l = likelihood(yclean, w, mu, sigma);
        while(true) {
            EMUpdates updates = emStep(yclean, w, mu, sigma);
            w = updates.w;
            mu[1] = updates.mu2;
            double lprime = likelihood(yclean, w, mu, sigma);
            i += 1;
            if (Math.abs(l - lprime) < 0.0001 || i > maxIterations) {
                break;
            }
            l = lprime;
        }
        return w[0];
    }

    public double reduceReadPairInfos(Iterator<ReadPairInfo> values, Map<Short, ReadGroupInfo> readGroupInfos) {
        List<Double> insertSizes = new ArrayList<Double>();
        double maxSD = 0;
        while (values.hasNext()) {
            ReadPairInfo rpi = values.next();
            int insertSize = rpi.insertSize;
            short readGroupId = rpi.readGroupId;
            ReadGroupInfo readGroupInfo = readGroupInfos.get(readGroupId);
            insertSizes.add((double) (insertSize - readGroupInfo.isize));
            if (readGroupInfo.isizeSD > maxSD) {
                maxSD = readGroupInfo.isizeSD;
            }
        }
        double[] initialW = new double[]{.5,.5};
        double[] initialMu = new double[]{0,1000};
        double[] insertSizeArray = new double[insertSizes.size()];
        for (int i = 0; i < insertSizes.size(); i++) {
            insertSizeArray[i] = insertSizes.get(i);
        }
        return estimateW(insertSizeArray, initialW, initialMu, maxSD);
    }
}
