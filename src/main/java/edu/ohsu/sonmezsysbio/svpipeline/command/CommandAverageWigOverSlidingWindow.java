package edu.ohsu.sonmezsysbio.svpipeline.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 6/6/11
 * Time: 3:42 PM
 */
@Parameters(separators = "=", commandDescription = "compute a transformed wig file averaged over a sliding window")
public class CommandAverageWigOverSlidingWindow implements SVPipelineCommand {
    public static final int WINDOW_SIZE_IN_LINES = 1000;
    @Parameter(names = {"--InFile"}, required = true)
    String inFile;

    @Parameter(names = {"--OutFile"}, required = true)
    String outFile;

    public void run() throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(new File(inFile)));
        BufferedWriter writer = new BufferedWriter(new FileWriter(new File(outFile)));

        try {
            String line;
            //LinkedList<Double> windowValues = new LinkedList<Double>();
            HashMap<Integer,Double> window = new HashMap<Integer, Double>();

            int numLines = 0;
            double windowTotal = 0;
            int lastPos = 0;

            while ((line = reader.readLine()) != null) {
                if (line.startsWith("track")) {
                    writer.write(line + "\n");
                    continue;
                }

                if (line.startsWith("variableStep")) {
                    writer.write(line + "\n");
                    continue;
                }

                String[] fields = line.split("\t");
                Integer pos = Integer.parseInt(fields[0]);
                if (pos - lastPos > 50) {
                    //System.err.println("hit a gap between " + lastPos + " and " + pos);
                }
                while (pos - lastPos > 50) {
                    lastPos = lastPos + 50;
                    //System.err.println("adding "+ lastPos + ", " + 0);
                    window.put(lastPos,0d);
                    if (window.keySet().size() > WINDOW_SIZE_IN_LINES) {
                        windowTotal = writeVal(writer, window, windowTotal, lastPos);
                    }
                }
                lastPos = pos;
                Double val = Double.parseDouble(fields[1]);
                //System.err.println("adding "+ pos + ", " + val);
                window.put(pos,val);
                windowTotal += val;


                if (window.keySet().size() > WINDOW_SIZE_IN_LINES) {
                    windowTotal = writeVal(writer, window, windowTotal, pos);
                }
                numLines++;
            }
        } finally {
            reader.close();
            writer.close();
        }
    }

    private double writeVal(BufferedWriter writer, HashMap<Integer, Double> window, double windowTotal, Integer pos) throws IOException {
        Double avg = windowTotal / WINDOW_SIZE_IN_LINES;

        if (! window.containsKey(pos -  (WINDOW_SIZE_IN_LINES / 2) * 50)) {
            System.err.println("Current position = " + pos + ", but did not have mid position " + (pos - (WINDOW_SIZE_IN_LINES / 2) * 50));
        }
        Double midVal = window.get(pos - (WINDOW_SIZE_IN_LINES / 2) * 50);

        Double modVal = midVal - avg;

        int positionToLeave = pos - WINDOW_SIZE_IN_LINES * 50;
        if (! window.containsKey(pos -  WINDOW_SIZE_IN_LINES * 50)) {
            System.err.println("Current position = " + pos + ", but did not have begin position " + (pos - WINDOW_SIZE_IN_LINES * 50));
            List positions = new ArrayList(window.keySet());
            Collections.sort(positions);
            System.err.println("lowest positions = " + positions.get(0) + ", " + positions.get(1));
        }

        writer.write(Integer.toString(pos - (WINDOW_SIZE_IN_LINES / 2) * 50)  + "\t" + modVal + "\n");


        Double leaving = window.get(positionToLeave);

        //Double leaving = windowValues.removeLast();
        windowTotal -= leaving;
        //System.err.println("removing " + positionToLeave + ", new total = " + windowTotal);
        window.remove(positionToLeave);
        return windowTotal;
    }


}
