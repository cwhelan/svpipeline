package edu.ohsu.sonmezsysbio.svpipeline.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.ohsu.sonmezsysbio.svpipeline.SVPipeline;
import edu.ohsu.sonmezsysbio.svpipeline.WigFileHelper;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hdfs.DistributedFileSystem;

import java.io.*;
import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 3/16/12
 * Time: 9:30 AM
 */
@Parameters(separators = "=", commandDescription = "Export Wig files and Bed file of deletions")
public class CommandExportWigAndBedFiles implements SVPipelineCommand {
    @Parameter(names = {"--inputHDFSDir"}, required = true)
    String inputHDFSDir;

    @Parameter(names = {"--outputPrefix"}, required = true)
    String outputPrefix;

    static class ReaderAndLine implements Comparable<ReaderAndLine> {
        protected BufferedReader reader;
        protected String chr;
        protected Long loc;
        protected Double nextValue;

        ReaderAndLine(BufferedReader reader, String chr, Long loc, Double nextValue) {
            this.reader = reader;
            this.chr = chr;
            this.loc = loc;            
            this.nextValue = nextValue;
        }

        public int compareTo(ReaderAndLine o) {
            if (chr.compareTo(o.chr) != 0) {
                return chr.compareTo(o.chr);
            } else {
                return loc.compareTo(o.loc);
            }
        }
    }

    public void run(Configuration conf) throws Exception {
        String pileupFileName = outputPrefix + "_piledup_deletion_scores.wig";
        String averagedFileName = outputPrefix + "_windowed_average_deletion_scores.wig";
        String bedFileName = outputPrefix + "_positive_score_regions.bed";
        
        File outputFile = new File(pileupFileName);
        if (! outputFile.createNewFile()) {
            System.err.println("Failed to create file " + outputFile);
            return;
        }

        System.err.println("Writing file " + pileupFileName);
        BufferedWriter outputFileWriter = new BufferedWriter(new FileWriter(outputFile));
        writePiledUpDeletionScores(conf, outputFileWriter, inputHDFSDir);
        outputFileWriter.close();

        System.err.println("Averaging scores over sliding window into " + averagedFileName);
        BufferedReader inFileReader = new BufferedReader(new FileReader(new File(pileupFileName)));
        BufferedWriter outFileWriter = new BufferedWriter(new FileWriter(new File(averagedFileName)));
        try {
            WigFileHelper.averageWigOverSlidingWindow(SVPipeline.RESOLUTION, SVPipeline.WINDOW_SIZE_IN_LINES, inFileReader, outFileWriter);
        } finally {
            inFileReader.close();
            outFileWriter.close();
        }

        System.err.println("Exporing regions with positive scores into " + bedFileName);
        BufferedReader averagedWigFileReader = new BufferedReader(new FileReader(new File(averagedFileName)));
        BufferedWriter bedFileWriter = new BufferedWriter(new FileWriter(new File(averagedFileName)));
        try {
            WigFileHelper.exportPositiveRegionsFromWig(outputPrefix, averagedWigFileReader, bedFileWriter);
        } finally {
            averagedWigFileReader.close();
            bedFileWriter.close();
        }

    }

    private void writePiledUpDeletionScores(Configuration conf, Writer outputFileWriter, String inputHDFSDir1) throws IOException {
        String currentChromosome = "";
        outputFileWriter.write("track type=wiggle_0 name=\"" + outputPrefix + " Deletion Scores\"\n");

        FileSystem dfs = DistributedFileSystem.get(conf);
        FileStatus[] stati = dfs.listStatus(new Path(inputHDFSDir1));
        if (stati == null) {
            System.err.println("Could not find input directory " + inputHDFSDir1);
            return;
        }
        LinkedList<ReaderAndLine> fileReaders = new LinkedList<ReaderAndLine>();
        for (FileStatus s : stati) {
            if (s.getPath().getName().startsWith("part")) {
                Path path = s.getPath();
                System.out.println(path);
                BufferedReader reader = new BufferedReader(new InputStreamReader(dfs.open(path)));
                String line = reader.readLine();
                String[] fields = line.split("\t");
                fileReaders.add(new ReaderAndLine(reader, fields[0], new Long(fields[1]), new Double(fields[2])));
            }
        }

        while (! fileReaders.isEmpty()) {
            ReaderAndLine minNextLine = null;
            for (ReaderAndLine readerAndLine : fileReaders) {
                if (minNextLine == null || minNextLine.compareTo(readerAndLine) > 0) {
                    minNextLine = readerAndLine;
                }
            }
            if (! currentChromosome.equals(minNextLine.chr)) {
                outputFileWriter.write("variableStep chrom=" + minNextLine.chr + " span=" + SVPipeline.RESOLUTION + "\n");
                currentChromosome = minNextLine.chr;
            }
            outputFileWriter.write(minNextLine.loc + "\t" + minNextLine.nextValue + "\n");
            String line = minNextLine.reader.readLine();
            if (line != null) {
                String[] fields = line.split("\t");
                minNextLine.chr = fields[0];
                minNextLine.loc = new Long(fields[1]);
                minNextLine.nextValue = new Double(fields[2]);
            } else {
                minNextLine.reader.close();
                fileReaders.remove(minNextLine);
            }
        }
    }
}
