package edu.ohsu.sonmezsysbio.cloudbreak.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.ohsu.sonmezsysbio.cloudbreak.Cloudbreak;
import edu.ohsu.sonmezsysbio.cloudbreak.file.DFSFacade;
import edu.ohsu.sonmezsysbio.cloudbreak.file.FaidxFileHelper;
import edu.ohsu.sonmezsysbio.cloudbreak.io.GMMResultsReaderAndLine;
import edu.ohsu.sonmezsysbio.cloudbreak.io.TextReaderAndLine;
import edu.ohsu.sonmezsysbio.cloudbreak.reducer.GMMScorerResults;
import edu.ohsu.sonmezsysbio.svpipeline.io.GenomicLocation;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hdfs.DistributedFileSystem;
import org.apache.hadoop.io.SequenceFile;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.PriorityQueue;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 3/16/12
 * Time: 9:30 AM
 */
@Parameters(separators = "=", commandDescription = "Export Wig files and Bed file of deletions")
public class CommandExportGMMResults implements CloudbreakCommand {
    @Parameter(names = {"--inputHDFSDir"}, required = true)
    String inputHDFSDir;

    @Parameter(names = {"--outputPrefix"}, required = true)
    String outputPrefix;

    @Parameter(names = {"--faidx"}, required = true)
    private String faidxFileName;

    @Parameter(names = {"--medianFilterWindow"})
    int medianFilterWindow = 1;

    @Parameter(names = {"--resolution"})
    int resolution = Cloudbreak.DEFAULT_RESOLUTION;

    public String getFaidxFileName() {
        return faidxFileName;
    }

    public void setFaidxFileName(String faidxFileName) {
        this.faidxFileName = faidxFileName;
    }

    public void run(Configuration conf) throws Exception {

        FaidxFileHelper faidx = new FaidxFileHelper(faidxFileName);

        String w0fileName = outputPrefix + "_w0.wig";
        String mu1fileName = outputPrefix + "_mu1.wig";
        String l1fileName = outputPrefix + "_l1.wig";
        String l2fileName = outputPrefix + "_l2.wig";
        String l1ffileName = outputPrefix + "_l1f.wig";
        String lrHetfileName = outputPrefix + "_lrHet.wig";
        String lrHomfileName = outputPrefix + "_lrHom.wig";

        String pileupBedFileName = outputPrefix + "_piledup_positive_score_regions.bed";

        // todo: refactor to hide all of these writers
        BufferedWriter w0outputFileWriter = createWriter(w0fileName);
        if (w0outputFileWriter == null) return;

        BufferedWriter mu1outputFileWriter = createWriter(mu1fileName);
        if (mu1outputFileWriter == null) return;

        BufferedWriter l1outputFileWriter = createWriter(l1fileName);
        if (l1outputFileWriter == null) return;

        BufferedWriter l2outputFileWriter = createWriter(l2fileName);
        if (l2outputFileWriter == null) return;

        BufferedWriter l1foutputFileWriter = createWriter(l1ffileName);
        if (l1foutputFileWriter == null) return;

        BufferedWriter lrHetOutputFileWriter = createWriter(lrHetfileName);
        if (lrHetOutputFileWriter == null) return;

        BufferedWriter lrHomOutputFileWriter = createWriter(lrHomfileName);
        if (lrHomOutputFileWriter == null) return;

        try {
            writeGMMResultWigFiles(conf, w0outputFileWriter, mu1outputFileWriter, l1outputFileWriter,
                    l2outputFileWriter, l1foutputFileWriter, lrHetOutputFileWriter, lrHomOutputFileWriter,
                    inputHDFSDir, faidx);
        } finally {
            w0outputFileWriter.close();
            mu1outputFileWriter.close();
            l1outputFileWriter.close();
            l2outputFileWriter.close();
            lrHetOutputFileWriter.close();
            lrHomOutputFileWriter.close();
        }


    }

    private BufferedWriter createWriter(String fileName) throws IOException {
        File w0outputFile = new File(fileName);
        if (! w0outputFile.createNewFile()) {
            System.err.println("Failed to create file " + w0outputFile);
            return null;
        }

        System.err.println("Writing file " + fileName);
        BufferedWriter outputFileWriter = new BufferedWriter(new FileWriter(w0outputFile));
        return outputFileWriter;
    }

    private void writeGMMResultWigFiles(Configuration conf, Writer w0outputFileWriter, BufferedWriter mu1outputFileWriter,
                                        BufferedWriter l1outputFileWriter, BufferedWriter l2outputFileWriter,
                                        BufferedWriter l1fOutputFileWriter, BufferedWriter lrHetOutputFileWriter,
                                        BufferedWriter lrHomOutputFileWriter, String inputHDFSDir1,
                                        FaidxFileHelper faix
    ) throws IOException {
        w0outputFileWriter.write("track type=wiggle_0 name=\"" + outputPrefix + " w0\"\n");
        mu1outputFileWriter.write("track type=wiggle_0 name=\"" + outputPrefix + " mu1\"\n");
        l1outputFileWriter.write("track type=wiggle_0 name=\"" + outputPrefix + " l1\"\n");
        l2outputFileWriter.write("track type=wiggle_0 name=\"" + outputPrefix + " l2\"\n");
        l1fOutputFileWriter.write("track type=wiggle_0 name=\"" + outputPrefix + " l1f\"\n");
        lrHetOutputFileWriter.write("track type=wiggle_0 name=\"" + outputPrefix + " lrHet\"\n");
        lrHomOutputFileWriter.write("track type=wiggle_0 name=\"" + outputPrefix + " lrHom\"\n");

        FileSystem dfs = DistributedFileSystem.get(conf);
        FileStatus[] stati = dfs.listStatus(new Path(inputHDFSDir1));
        if (stati == null) {
            System.err.println("Could not find input directory " + inputHDFSDir1);
            return;
        }

        List<Path> inputStreams = new ArrayList<Path>();
        for (FileStatus s : stati) {
            if (s.getPath().getName().startsWith("part")) {
                Path path = s.getPath();
                System.out.println(path);
                inputStreams.add(path);
//                inputStreams.add(dfs.open(path));
            }
        }

        mergeSortedInputStreams(new DFSFacade(dfs, conf), w0outputFileWriter, mu1outputFileWriter, l1outputFileWriter,
                l2outputFileWriter, l1fOutputFileWriter, lrHetOutputFileWriter, lrHomOutputFileWriter,
                faix, inputStreams);
    }

    public void mergeSortedInputStreams(DFSFacade dfsFacade, Writer w0outputFileWriter, BufferedWriter mu1outputFileWriter,
                                        BufferedWriter l1outputFileWriter, BufferedWriter l2outputFileWriter,
                                        BufferedWriter l1fOutputFileWriter, BufferedWriter lrHetOutputFileWriter,
                                        BufferedWriter lrHomOutputFileWriter, FaidxFileHelper faix,
                                        List<Path> paths) throws IOException {
        short currentChromosome = -1;
        PriorityQueue<GMMResultsReaderAndLine> fileReaders = new PriorityQueue<GMMResultsReaderAndLine>();
        for (Path path : paths) {
            SequenceFile.Reader reader = new SequenceFile.Reader(dfsFacade.dfs, path, dfsFacade.conf);
            GenomicLocation gl = new GenomicLocation();
            GMMScorerResults results = new GMMScorerResults();
            reader.next(gl, results);
            //System.err.println("Read " + gl);
            fileReaders.add(new GMMResultsReaderAndLine(reader, gl, results));
        }

        while (! fileReaders.isEmpty()) {
            GMMResultsReaderAndLine minNextLine = fileReaders.poll();
            if (currentChromosome != minNextLine.getGenomicLocation().chromosome) {
                writeChromHeader(w0outputFileWriter, faix, minNextLine);
                writeChromHeader(mu1outputFileWriter, faix, minNextLine);
                writeChromHeader(l1outputFileWriter, faix, minNextLine);
                writeChromHeader(l2outputFileWriter, faix, minNextLine);
                writeChromHeader(l1fOutputFileWriter, faix, minNextLine);
                writeChromHeader(lrHetOutputFileWriter, faix, minNextLine);
                writeChromHeader(lrHomOutputFileWriter, faix, minNextLine);

                currentChromosome = minNextLine.getGenomicLocation().chromosome;
            }

            w0outputFileWriter.write(minNextLine.getGenomicLocation().pos + "\t" + minNextLine.getNextValue().w0 + "\n");
            mu1outputFileWriter.write(minNextLine.getGenomicLocation().pos + "\t" + minNextLine.getNextValue().mu2 + "\n");
            l1outputFileWriter.write(minNextLine.getGenomicLocation().pos + "\t" + minNextLine.getNextValue().nodelOneComponentLikelihood + "\n");
            l2outputFileWriter.write(minNextLine.getGenomicLocation().pos + "\t" + minNextLine.getNextValue().twoComponentLikelihood + "\n");
            l1fOutputFileWriter.write(minNextLine.getGenomicLocation().pos + "\t" + minNextLine.getNextValue().oneFreeComponentLikelihood + "\n");
            lrHetOutputFileWriter.write(minNextLine.getGenomicLocation().pos + "\t" + minNextLine.getNextValue().lrHeterozygous + "\n");
            lrHomOutputFileWriter.write(minNextLine.getGenomicLocation().pos + "\t" + minNextLine.getNextValue().lrHomozygous + "\n");

            boolean gotLine;
            gotLine = readNextDataLine((GMMResultsReaderAndLine) minNextLine);

            if (gotLine) {
                fileReaders.add(minNextLine);
            } else {
                minNextLine.closeInput();
            }
        }
    }

    private void writeChromHeader(Writer w0outputFileWriter, FaidxFileHelper faix, GMMResultsReaderAndLine minNextLine) throws IOException {
        w0outputFileWriter.write("variableStep chrom=" + faix.getNameForChromKey(minNextLine.getGenomicLocation().chromosome) + " span=" + resolution + "\n");
    }

    private boolean readNextTextLine(FaidxFileHelper faix, TextReaderAndLine minNextLine) throws IOException {
        String line = minNextLine.getDataInput().readLine();
        boolean gotLine = false;
        if (line != null) {
            String[] fields = line.split("\t");
            minNextLine.setGenomicLocation(new GenomicLocation(faix.getKeyForChromName(fields[0]), new Integer(fields[1])));
            minNextLine.setNextValue(new Double(fields[2]));
            gotLine = true;
        }
        return gotLine;
    }

    private boolean readNextDataLine(GMMResultsReaderAndLine minNextLine) throws IOException {
        GenomicLocation gl = new GenomicLocation();
        GMMScorerResults results = new GMMScorerResults();
        boolean gotLine = minNextLine.getReader().next(gl, results);
        minNextLine.setGenomicLocation(gl);
        minNextLine.setNextValue(results);
        return gotLine;
    }

}
