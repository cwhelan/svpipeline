package edu.ohsu.sonmezsysbio.svpipeline.mapper;

import edu.ohsu.sonmezsysbio.svpipeline.*;
import edu.ohsu.sonmezsysbio.svpipeline.io.GenomicLocation;
import edu.ohsu.sonmezsysbio.svpipeline.io.ReadPairInfo;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.*;

import java.io.IOException;
import java.util.List;

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
    private String faidxFileName;
    FaidxFileHelper faix;

    // for debugging, restrict output to a particular region
    private String chromosomeFilter;
    private Long startFilter;
    private Long endFilter;

    public FaidxFileHelper getFaix() {
        return faix;
    }

    public void setFaix(FaidxFileHelper faix) {
        this.faix = faix;
    }

    public String getChromosomeFilter() {
        return chromosomeFilter;
    }

    public void setChromosomeFilter(String chromosomeFilter) {
        this.chromosomeFilter = chromosomeFilter;
    }

    public Long getStartFilter() {
        return startFilter;
    }

    public void setStartFilter(Long startFilter) {
        this.startFilter = startFilter;
    }

    public Long getEndFilter() {
        return endFilter;
    }

    public void setEndFilter(Long endFilter) {
        this.endFilter = endFilter;
    }

    public boolean isMatePairs() {
        return matePairs;
    }

    public void setMatePairs(boolean matePairs) {
        this.matePairs = matePairs;
    }

    public Integer getMaxInsertSize() {
        return maxInsertSize;
    }

    public void setMaxInsertSize(Integer maxInsertSize) {
        this.maxInsertSize = maxInsertSize;
    }

    public PairedAlignmentScorer getScorer() {
        return scorer;
    }

    public void setScorer(PairedAlignmentScorer scorer) {
        this.scorer = scorer;
    }

    public void map(LongWritable key, Text value, OutputCollector<GenomicLocation, ReadPairInfo> output, Reporter reporter) throws IOException {
        String line = value.toString();
        int firstTabIndex = line.indexOf('\t');
        String lineValues = line.substring(firstTabIndex + 1);

        String[] readAligments = lineValues.split(SVPipeline.READ_SEPARATOR);
        String read1AlignmentsString = readAligments[0];
        String[] read1Alignments = read1AlignmentsString.split(SVPipeline.ALIGNMENT_SEPARATOR);
        List<NovoalignNativeRecord> read1AlignmentRecords = NovoalignSingleEndMapperHelper.parseAlignmentsIntoRecords(read1Alignments);

        String read2AlignmentsString = readAligments[1];
        String[] read2Alignments = read2AlignmentsString.split(SVPipeline.ALIGNMENT_SEPARATOR);
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


        int genomicWindow = insertSize +
                leftRead.getPosition() % SVPipeline.RESOLUTION +
                SVPipeline.RESOLUTION - rightRead.getPosition() % SVPipeline.RESOLUTION;

        ReadPairInfo readPairInfo = new ReadPairInfo(insertSize, scorer.probabilityMappingIsCorrect(endPosterior1, endPosterior2));

        for (int i = 0; i <= genomicWindow; i = i + SVPipeline.RESOLUTION) {
            Short chromosome = faix.getKeyForChromName(record1.getChromosomeName());
            if (chromosome == null) {
                throw new RuntimeException("Bad chromosome in record: " + record1.getChromosomeName());
            }

            int pos = genomeOffset + i;

            if (getChromosomeFilter() != null) {
                if (! record1.getChromosomeName().equals(getChromosomeFilter()) ||
                    pos < getStartFilter() || pos > getEndFilter()) {
                    return;
                }
            }

            GenomicLocation genomicLocation = new GenomicLocation(chromosome, pos);
            output.collect(genomicLocation, readPairInfo);

        }

    }


    public void configure(JobConf job) {
        super.configure(job);
        matePairs = Boolean.parseBoolean(job.get("pileupDeletionScore.isMatePairs"));
        scorer = new ProbabilisticPairedAlignmentScorer();

        faidxFileName = job.get("alignment.faidx");
        faix = new FaidxFileHelper(faidxFileName);

        if (job.get("alignments.filterchr") != null) {
            setChromosomeFilter(job.get("alignments.filterchr"));
            setStartFilter(Long.parseLong(job.get("alignments.filterstart")));
            setEndFilter(Long.parseLong(job.get("alignments.filterend")));
        }
    }
}
