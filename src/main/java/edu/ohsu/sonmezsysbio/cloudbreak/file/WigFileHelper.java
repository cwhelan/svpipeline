package edu.ohsu.sonmezsysbio.cloudbreak.file;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 3/19/12
 * Time: 3:28 PM
 */
public class WigFileHelper {

    public static void averageWigOverSlidingWindow(int resolution, int windowSizeToAverageOver, BufferedReader inFileReader, BufferedWriter outFileWriter) throws IOException {

            String line;

            HashMap<Integer,Double> window;
            double windowTotal;
            int lastPos;

            window = new HashMap<Integer, Double>();
            windowTotal = 0;
            lastPos = 0;

            while ((line = inFileReader.readLine()) != null) {
                if (line.startsWith("track")) {
                    outFileWriter.write(line + "\n");
                    continue;
                }

                if (line.startsWith("variableStep")) {
                    outFileWriter.write(line + "\n");

                    window = new HashMap<Integer, Double>();
                    windowTotal = 0;
                    lastPos = 0;

                    continue;
                }

                String[] fields = line.split("\t");
                Integer pos = Integer.parseInt(fields[0]);
                if (pos - lastPos > resolution) {
                    // System.err.println("hit a gap between " + lastPos + " and " + pos);
                }
                while (pos - lastPos > resolution) {
                    lastPos = lastPos + resolution;
                    //System.err.println("adding "+ lastPos + ", " + 0);
                    window.put(lastPos,0d);
                    if (window.keySet().size() > windowSizeToAverageOver) {
                        windowTotal = writeVal(outFileWriter, window, windowTotal, lastPos, windowSizeToAverageOver, resolution);
                    }
                }
                lastPos = pos;
                Double val = Double.parseDouble(fields[1]);
                //System.err.println("adding "+ pos + ", " + val);
                window.put(pos,val);
                windowTotal += val;


                if (window.keySet().size() > windowSizeToAverageOver) {
                    windowTotal = writeVal(outFileWriter, window, windowTotal, pos, windowSizeToAverageOver, resolution);
                }
            }
    }

    private static double writeVal(BufferedWriter writer, HashMap<Integer, Double> window, double windowTotal, Integer pos, int windowToAverageOver, int resolution) throws IOException {
        Double avg = windowTotal / windowToAverageOver;

        if (! window.containsKey(pos -  (windowToAverageOver / 2) * resolution)) {
            System.err.println("Current position = " + pos + ", but did not have mid position " + (pos - (windowToAverageOver / 2) * resolution));
        }
        Double midVal = window.get(pos - (windowToAverageOver / 2) * resolution);

        Double modVal = midVal - avg;

        int positionToLeave = pos - windowToAverageOver * resolution;
        if (! window.containsKey(pos -  windowToAverageOver * resolution)) {
            System.err.println("Current position = " + pos + ", but did not have begin position " + (pos - windowToAverageOver * resolution));
            List<Integer> positions = new ArrayList<Integer>(window.keySet());
            Collections.sort(positions);
            System.err.println("lowest positions = " + positions.get(0) + ", " + positions.get(1));
        }

        writer.write(Integer.toString(pos - (windowToAverageOver / 2) * resolution)  + "\t" + modVal + "\n");


        Double leaving = window.get(positionToLeave);

        //Double leaving = windowValues.removeLast();
        windowTotal -= leaving;
        //System.err.println("removing " + positionToLeave + ", new total = " + windowTotal);
        window.remove(positionToLeave);
        return windowTotal;
    }

    public static void exportRegionsOverThresholdFromWig(String outputPrefix, BufferedReader averagedWigFileReader, BufferedWriter bedFileWriter, double threshold, FaidxFileHelper faidx, int medianFilterWindow) throws IOException {
        String trackName = outputPrefix + " peaks over " + threshold;
        bedFileWriter.write("track name = \"" + trackName + "\"\n");

        String line;
        String currentChromosome = "";

        double[] values = null;
        int resolution = 0;
        int peakNum = 1;

        while ((line = averagedWigFileReader.readLine()) != null) {
            if (line.startsWith("track")) {
                continue;
            }
            if (line.startsWith("variableStep")) {

                if (values != null) {
                    double[] filteredVals = medianFilterValues(values, medianFilterWindow, threshold);
                    peakNum = writePositiveRegions(filteredVals, bedFileWriter, currentChromosome, faidx, resolution, peakNum);
                }
                currentChromosome = line.split(" ")[1].split("=")[1];
                resolution = Integer.valueOf(line.split(" ")[2].split("=")[1]);
                int numTiles = (int) Math.ceil(((double) faidx.getLengthForChromName(currentChromosome)) / resolution);
                values = new double[numTiles];

            } else {
                String[] fields = line.split("\t");
                if (fields.length < 2) {
                    throw new RuntimeException("Failed to parse line: " + line);
                }
                long pos = Long.valueOf(fields[0]);
                if (pos > faidx.getLengthForChromName(currentChromosome)) continue;
                double val = Double.valueOf(fields[1]);
                int tileNum = (int) pos / resolution;
                values[tileNum] = val;
            }
        }
        double[] filteredVals = medianFilterValues(values, medianFilterWindow, threshold);
        writePositiveRegions(filteredVals, bedFileWriter, currentChromosome, faidx, resolution, peakNum);

    }

    private static double[] medianFilterValues(double[] values, int medianFilterWindow, double threshold) {
        double[] filteredValues = new double[values.length];
        int idx = 0;
        while (idx <= medianFilterWindow / 2) {
            filteredValues[idx] = values[idx] > threshold ? values[idx] : 0.0;
            idx++;
        }

        while (idx < filteredValues.length - medianFilterWindow / 2) {
            double[] filterWindow = Arrays.copyOfRange(values, idx - medianFilterWindow / 2, idx + medianFilterWindow / 2 + 1);
            Arrays.sort(filterWindow);
            double medianValue = filterWindow[medianFilterWindow / 2];
            if (medianValue > threshold) {
                filteredValues[idx] = values[idx];
            } else {
                filteredValues[idx] = 0;
            }
            idx++;
        }

        while (idx < filteredValues.length) {
            filteredValues[idx] = values[idx] > threshold ? values[idx] : 0.0;
            idx++;
        }
        return filteredValues;
    }

    private static int writePositiveRegions(double[] filteredVals, BufferedWriter bedFileWriter, String currentChromosome, FaidxFileHelper faidx, int resolution, int peakNum) throws IOException {
        boolean inPositivePeak = false;
        long peakStart = 0;
        int idx = 0;
        double peakMax = 0;

        while (idx < filteredVals.length) {
            long pos = idx * resolution;

            if (filteredVals[idx] > 0) {
                if (!inPositivePeak) {
                    peakStart = pos;
                    inPositivePeak = true;
                }
                peakMax = Math.max(peakMax, filteredVals[idx]);
            } else {
                if (inPositivePeak) {
                    long endPosition = pos - 1;
                    bedFileWriter.write(currentChromosome + "\t" + peakStart + "\t" + endPosition + "\t" + peakNum + "\t" + peakMax + "\n");
                    peakNum += 1;
                    inPositivePeak = false;
                    peakMax = 0;
                }
            }
            idx = idx + 1;
        }
        if (inPositivePeak) {
            long endPosition = faidx.getLengthForChromName(currentChromosome) - 1;
            if (endPosition < peakStart) return peakNum;
            bedFileWriter.write(currentChromosome + "\t" + peakStart + "\t" + endPosition + "\t" + peakNum + "\t" + 0 + "\n");
            peakNum += 1;
        }
        return peakNum;
    }

}
