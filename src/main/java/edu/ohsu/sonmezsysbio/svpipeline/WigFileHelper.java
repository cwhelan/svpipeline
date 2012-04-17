package edu.ohsu.sonmezsysbio.svpipeline;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

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

    public static void exportPositiveRegionsFromWig(String outputPrefix, BufferedReader averagedWigFileReader, BufferedWriter bedFileWriter, double threshold, FaidxFileHelper faidx) throws IOException {
        String trackName = outputPrefix + " peaks over " + threshold;
        bedFileWriter.write("track name = \"" + trackName + "\"\n");

        boolean inPositivePeak = false;
        double maxScoreInPeak = 0;
        double sumPeakScores = 0;
        long peakStart = 0;
        int peakNum = 1;

        String line;
        String currentChromosome = "";

        while ((line = averagedWigFileReader.readLine()) != null) {
            if (line.startsWith("track")) {
                continue;
            }
            
            if (line.startsWith("variableStep")) {
                if (inPositivePeak) {
                    long endPosition = faidx.getLengthForChromName(currentChromosome) - 1;
                    if (endPosition < peakStart) continue;
                    bedFileWriter.write(currentChromosome + "\t" + peakStart + "\t" + endPosition + "\t" + peakNum + "\t" + sumPeakScores + "\n");
                    peakNum += 1;
                }
                // reset counters
                inPositivePeak = false;
                maxScoreInPeak = 0;
                sumPeakScores = 0;
                peakStart = 0;

                currentChromosome = line.split(" ")[1].split("=")[1];
                continue;
            }
        
            String[] fields = line.split("\t");
            long pos = Long.valueOf(fields[0]);
            double val = Double.valueOf(fields[1]);

            if (val > threshold) {
                if (inPositivePeak) {
                    if (val > maxScoreInPeak) {
                        maxScoreInPeak = val;
                    }
                    sumPeakScores += val;
                } else {
                    peakStart = pos;
                    inPositivePeak = true;
                    maxScoreInPeak = val;
                    sumPeakScores = val;
                }
                
            } else {
                if (inPositivePeak) {
                    long endPosition = pos - 1;
                    bedFileWriter.write(currentChromosome + "\t" + peakStart + "\t" + endPosition + "\t" + peakNum + "\t" + sumPeakScores + "\n");
                    peakNum += 1;
                    inPositivePeak = false;
                }
            }
        }
    }
}
