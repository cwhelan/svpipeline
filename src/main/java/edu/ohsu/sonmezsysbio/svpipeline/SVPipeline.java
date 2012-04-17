package edu.ohsu.sonmezsysbio.svpipeline;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import edu.ohsu.sonmezsysbio.svpipeline.command.*;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;
import org.apache.hadoop.mapred.JobConf;
import org.apache.hadoop.util.Tool;
import org.apache.hadoop.util.ToolRunner;

public class SVPipeline extends Configured implements Tool
{

    public static final String ALIGNMENT_SEPARATOR = "\tSVP_ALIGNMENT\t";
    public static final String READ_SEPARATOR = "\tSVP_READ\t";
    public static final int RESOLUTION = 100;
    public static final int WINDOW_SIZE_IN_LINES = 1000;

    public static void main(String[] args) throws Exception {
        int res = ToolRunner.run(new Configuration(), new SVPipeline(), args);

        System.exit(res);
    }

    public int run(String[] args) throws Exception {
        JCommander jc = buildJCommander();

        String parsedCommand = null;

        try {
            jc.parse(args);

            parsedCommand = jc.getParsedCommand();

            if (parsedCommand == null) {
                jc.usage();
                return 1;
            }
            SVPipelineCommand command = (SVPipelineCommand) jc.getCommands().get(parsedCommand).getObjects().get(0);
            command.run(getConf());

            return 0;
        } catch (ParameterException pe) {
            System.err.println(pe.getMessage());
            return 1;
        } catch (Exception e) {
            e.printStackTrace();
            return 1;
        }
    }

    protected static JCommander buildJCommander() {
        JCommander jc = new JCommander(new CommanderMain());
        CommandNovoalignMatePair novoalignMatePair = new CommandNovoalignMatePair();
        jc.addCommand("novoalignMatePair", novoalignMatePair);
        CommandPileupDeletionScores pileupDeletionScores = new CommandPileupDeletionScores();
        jc.addCommand("pileupDeletionScores", pileupDeletionScores);

        CommandReadPairedEndFilesIntoHDFS readFiles = new CommandReadPairedEndFilesIntoHDFS();
        jc.addCommand("readPairedEndFilesIntoHDFS", readFiles);
        CommandNovoalignSingleEnds singleEnds  = new CommandNovoalignSingleEnds();
        jc.addCommand("alignSingleEnds", singleEnds);
        CommandPileupSingleEndDeletionScores pileupSingleEndDeletionScores = new CommandPileupSingleEndDeletionScores();
        jc.addCommand("pileupSingleEndDeletionScores", pileupSingleEndDeletionScores);

        CommandIncrementalUpdateSingleEndDeletionScores incrementalUpdateSingleEndDeletionScores = new CommandIncrementalUpdateSingleEndDeletionScores();
        jc.addCommand("incrementalUpdateSingleEndDeletionScores", incrementalUpdateSingleEndDeletionScores);

        jc.addCommand("splitDeflatedOutput", new CommandSplitDeflatedOutput());

        CommandAverageWigOverSlidingWindow averageWigOverSlidingWindow = new CommandAverageWigOverSlidingWindow();
        jc.addCommand("averageWigOverSlidingWindow", averageWigOverSlidingWindow);

        CommandExportWigAndBedFiles exportWigAndBedFiles = new CommandExportWigAndBedFiles();
        jc.addCommand("exportWigAndBedFiles", exportWigAndBedFiles);

        CommandDumpReadsWithScores dumpReadsWithScores = new CommandDumpReadsWithScores();
        jc.addCommand("dumpReadsWithScores", dumpReadsWithScores);

        CommandExtractPositiveRegionsFromWig commandExtractPositiveRegionsFromWig = new CommandExtractPositiveRegionsFromWig();
        jc.addCommand("extractPositiveRegionsFromWig", commandExtractPositiveRegionsFromWig);

        CommandDebugReadPairInfo commandDebugReadPairInfo = new CommandDebugReadPairInfo();
        jc.addCommand("debugReadPairInfo", commandDebugReadPairInfo);

        CommandSummarizeAlignments commandSummarizeAlignments = new CommandSummarizeAlignments();
        jc.addCommand("summarizeAlignments", commandSummarizeAlignments);

        jc.setProgramName("SVPipeline");
        return jc;
    }
}
