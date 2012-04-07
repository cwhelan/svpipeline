package edu.ohsu.sonmezsysbio.svpipeline.mapper;

import edu.ohsu.sonmezsysbio.svpipeline.NovoalignNativeRecord;
import edu.ohsu.sonmezsysbio.svpipeline.PairedAlignmentScorer;
import edu.ohsu.sonmezsysbio.svpipeline.ProbabilisticPairedAlignmentScorer;
import edu.ohsu.sonmezsysbio.svpipeline.SVPipeline;
import edu.ohsu.sonmezsysbio.svpipeline.io.GenomicLocation;
import edu.ohsu.sonmezsysbio.svpipeline.io.ReadPairInfo;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/6/12
 * Time: 1:03 PM
 */
public class SingleEndAlignmentsToReadPairInfoMapper extends MapReduceBase implements Mapper<LongWritable, Text, GenomicLocation, ReadPairInfo> {

    private boolean matePairs;
    private Integer maxInsertSize = 500000;
    private PairedAlignmentScorer scorer;
    private Map<String, Short> chromosomeKeys;
    private String faidxFileName;

    public void map(LongWritable key, Text value, OutputCollector<GenomicLocation, ReadPairInfo> output, Reporter reporter) throws IOException {
        String line = value.toString();
        int firstTabIndex = line.indexOf('\t');
        String lineValues = line.substring(firstTabIndex + 1);

        String[] readAligments = lineValues.split(SVPipeline.READ_SEPARATOR);
        String read1AlignmentsString = readAligments[0];
        String[] read1Alignments = read1AlignmentsString.split(SVPipeline.ALINGMENT_SEPARATOR);
        List<NovoalignNativeRecord> read1AlignmentRecords = NovoalignSingleEndMapperHelper.parseAlignmentsIntoRecords(read1Alignments);

        String read2AlignmentsString = readAligments[1];
        String[] read2Alignments = read2AlignmentsString.split(SVPipeline.ALINGMENT_SEPARATOR);
        List<NovoalignNativeRecord> read2AlignmentRecords = NovoalignSingleEndMapperHelper.parseAlignmentsIntoRecords(read2Alignments);

        emitReadPairInfoForAllPairs(read1AlignmentRecords, read2AlignmentRecords, output);
    }

    private void emitReadPairInfoForAllPairs(List<NovoalignNativeRecord> read1AlignmentRecords, List<NovoalignNativeRecord> read2AlignmentRecords, OutputCollector<GenomicLocation, ReadPairInfo> output) throws IOException {
        for (NovoalignNativeRecord record1 : read1AlignmentRecords) {
            for (NovoalignNativeRecord record2 : read2AlignmentRecords) {
                emitReadPairInfoForPair(record1, record2, output);
            }
        }
    }

    private void emitReadPairInfoForPair(NovoalignNativeRecord record1, NovoalignNativeRecord record2, OutputCollector<GenomicLocation, ReadPairInfo> output) throws IOException {

        // todo: not handling translocations for now
        if (! record1.getChromosomeName().equals(record2.getChromosomeName())) return;

        int endPosterior1 = record1.getPosteriorProb();
        int endPosterior2 = record2.getPosteriorProb();

        int insertSize;
        NovoalignNativeRecord leftRead = record1.getPosition() < record2.getPosition() ?
                record1 : record2;
        NovoalignNativeRecord rightRead = record1.getPosition() < record2.getPosition() ?
                record2 : record1;

        // todo: not handling inversions for now
        if (!scorer.validateMappingOrientations(record1, record2, matePairs)) return;

        insertSize = rightRead.getPosition() + rightRead.getSequence().length() - leftRead.getPosition();

        if (! scorer.validateInsertSize(insertSize, record1.getReadId(), maxInsertSize)) return;

        int genomeOffset = leftRead.getPosition() - leftRead.getPosition() % SVPipeline.RESOLUTION;

        insertSize = insertSize + leftRead.getPosition() % SVPipeline.RESOLUTION + SVPipeline.RESOLUTION - rightRead.getPosition() % SVPipeline.RESOLUTION;
        ReadPairInfo readPairInfo = new ReadPairInfo(insertSize, scorer.probabilityMappingIsCorrect(endPosterior1, endPosterior2));

        for (int i = 0; i <= insertSize; i = i + SVPipeline.RESOLUTION) {
            Short chromosome = chromosomeKeys.get(record1.getChromosomeName());
            if (chromosome == null) {
                throw new RuntimeException("Bad chromosome in record: " + record1.getChromosomeName());
            }
            GenomicLocation genomicLocation = new GenomicLocation(chromosome, genomeOffset + i);
            output.collect(genomicLocation, readPairInfo);
        }

    }


    public void configure(JobConf job) {
        super.configure(job);
        matePairs = Boolean.parseBoolean(job.get("pileupDeletionScore.isMatePairs"));
        scorer = new ProbabilisticPairedAlignmentScorer();

        faidxFileName = job.get("alignment.faidx");
        try {
            chromosomeKeys = readFaidxFile(faidxFileName);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private Map<String, Short> readFaidxFile(String faidxFileName) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(faidxFileName));
        return readFaidx(reader);
    }

    public Map<String, Short> readFaidx(BufferedReader bufferedReader) throws IOException {
        String line;
        short chrIdx = 0;
        Map<String, Short> chrTable = new HashMap<String, Short>();
        while ((line = bufferedReader.readLine()) != null) {
            String chromosomeName = line.split("\\s+")[0];
            chrTable.put(chromosomeName, chrIdx);
            chrIdx++;
        }
        return chrTable;
    }
}
