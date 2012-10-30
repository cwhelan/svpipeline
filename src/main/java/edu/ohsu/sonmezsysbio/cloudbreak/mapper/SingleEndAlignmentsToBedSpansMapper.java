package edu.ohsu.sonmezsysbio.cloudbreak.mapper;

import edu.ohsu.sonmezsysbio.cloudbreak.*;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.JobConf;
import org.apache.hadoop.mapred.Mapper;
import org.apache.hadoop.mapred.OutputCollector;
import org.apache.hadoop.mapred.Reporter;

import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/23/11
 * Time: 10:12 AM
 */
public class SingleEndAlignmentsToBedSpansMapper extends CloudbreakMapReduceBase implements Mapper<Text, Text, Text, Text> {

    private boolean matePairs;
    private Integer maxInsertSize = 500000;

    // region to dump spans over
    private String chromosome;
    private int regionStart;
    private int regionEnd;
    private PairedAlignmentScorer scorer;
    private int minScore = -1;

    public Integer getMaxInsertSize() {
        return maxInsertSize;
    }

    public void setMaxInsertSize(Integer maxInsertSize) {
        this.maxInsertSize = maxInsertSize;
    }

    public boolean isMatePairs() {
        return matePairs;
    }

    public void setMatePairs(boolean matePairs) {
        this.matePairs = matePairs;
    }

    public String getChromosome() {
        return chromosome;
    }

    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    public int getRegionStart() {
        return regionStart;
    }

    public void setRegionStart(int regionStart) {
        this.regionStart = regionStart;
    }

    public int getRegionEnd() {
        return regionEnd;
    }

    public void setRegionEnd(int regionEnd) {
        this.regionEnd = regionEnd;
    }

    public int getMinScore() {
        return minScore;
    }

    public void setMinScore(int minScore) {
        this.minScore = minScore;
    }

    @Override
    public void configure(JobConf job) {
        super.configure(job);

        maxInsertSize = Integer.parseInt(job.get("pileupDeletionScore.maxInsertSize"));
        matePairs = Boolean.parseBoolean(job.get("pileupDeletionScore.isMatePairs"));
        parseRegion(job.get("pileupDeletionScore.region"));

        minScore = Integer.parseInt(job.get("pileupDeletionScore.minScore"));

        scorer = new ProbabilisticPairedAlignmentScorer();
    }

    protected void parseRegion(String region) {
        chromosome = region.split(":")[0];
        regionStart = Integer.parseInt(region.split(":")[1].split("-")[0]);
        regionEnd = Integer.parseInt(region.split(":")[1].split("-")[1]);
    }

    public void map(Text key, Text value, OutputCollector<Text, Text> output, Reporter reporter)
            throws IOException {

        String line = value.toString();
        //System.err.println("LINE: " + line);

        ReadPairAlignments readPairAlignments = parsePairAlignmentLine(line);
        emitDeletionScoresForAllPairs(readPairAlignments, output);
    }

    private void emitDeletionScoresForAllPairs(ReadPairAlignments readPairAlignments, OutputCollector<Text, Text> output) throws IOException {
        for (AlignmentRecord record1 : readPairAlignments.getRead1Alignments()) {
            if (getMinScore() != -1) {
                if (record1.getAlignmentScore() < getMinScore()) {
                    continue;
                }
            }
            for (AlignmentRecord record2 : readPairAlignments.getRead2Alignments()) {
                if (record2.getAlignmentScore() < getMinScore()) {
                    continue;
                }
                emitBedSpanForPair(record1, record2, readPairAlignments, output);
            }
        }
    }

    public void emitBedSpanForPair(AlignmentRecord record1, AlignmentRecord record2, ReadPairAlignments readPairAlignments, OutputCollector<Text, Text> output) throws IOException {

        // todo: not handling translocations for now
        if (! record1.getChromosomeName().equals(record2.getChromosomeName())) return;

        // todo: not handling inversions for now
        if (!scorer.validateMappingOrientations(record1, record2, isMatePairs())) return;

        int insertSize;

        AlignmentRecord leftRead = record1.getPosition() < record2.getPosition() ?
                record1 : record2;
        AlignmentRecord rightRead = record1.getPosition() < record2.getPosition() ?
                record2 : record1;

        if (! (record1.getChromosomeName().equals(chromosome) &&
               leftRead.getPosition() < regionEnd && rightRead.getPosition() > regionStart))
            return;

        insertSize = rightRead.getPosition() + rightRead.getSequenceLength() - leftRead.getPosition();

        if (! scorer.validateInsertSize(insertSize, record1.getReadId(), maxInsertSize)) return;

        double pMappingCorrect = alignmentReader.probabilityMappingIsCorrect(record1, record2);

        output.collect(new Text(leftRead.getReadId()),
                new Text(leftRead.getChromosomeName() + "\t" + leftRead.getPosition() + "\t" + rightRead.getPosition() + "\t" + leftRead.getReadId() + "\t" + insertSize + "\t" + pMappingCorrect));

    }
}
