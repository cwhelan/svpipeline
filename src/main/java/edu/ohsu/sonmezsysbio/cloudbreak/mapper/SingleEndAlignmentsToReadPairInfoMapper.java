package edu.ohsu.sonmezsysbio.cloudbreak.mapper;

import edu.ohsu.sonmezsysbio.cloudbreak.*;
import edu.ohsu.sonmezsysbio.cloudbreak.file.BigWigFileHelper;
import edu.ohsu.sonmezsysbio.cloudbreak.file.FaidxFileHelper;
import edu.ohsu.sonmezsysbio.cloudbreak.file.GFFFileHelper;
import edu.ohsu.sonmezsysbio.cloudbreak.file.ReadGroupInfoFileHelper;
import edu.ohsu.sonmezsysbio.cloudbreak.io.ReadPairInfo;
import edu.ohsu.sonmezsysbio.svpipeline.io.GenomicLocation;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.JobConf;
import org.apache.hadoop.mapred.Mapper;
import org.apache.hadoop.mapred.OutputCollector;
import org.apache.hadoop.mapred.Reporter;

import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/6/12
 * Time: 1:03 PM
 */
public class SingleEndAlignmentsToReadPairInfoMapper extends CloudbreakMapReduceBase implements Mapper<LongWritable, Text, GenomicLocation, ReadPairInfo> {

    private boolean matePairs;
    private Integer maxInsertSize = Cloudbreak.DEFAULT_MAX_INSERT_SIZE;
    private PairedAlignmentScorer scorer;
    private String faidxFileName;
    FaidxFileHelper faix;

    // for debugging, restrict output to a particular region
    private String chromosomeFilter;
    private Long startFilter;
    private Long endFilter;
    private GFFFileHelper exclusionRegions;
    private BigWigFileHelper mapabilityWeighting;
    private double targetIsize;
    private double targetIsizeSD;
    private Short readGroupId;

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

    public GFFFileHelper getExclusionRegions() {
        return exclusionRegions;
    }

    public void setExclusionRegions(GFFFileHelper exclusionRegions) {
        this.exclusionRegions = exclusionRegions;
    }

    public Short getReadGroupId() {
        return readGroupId;
    }

    public void setReadGroupId(Short readGroupId) {
        this.readGroupId = readGroupId;
    }

    public void map(LongWritable key, Text value, OutputCollector<GenomicLocation, ReadPairInfo> output, Reporter reporter) throws IOException {
        String line = value.toString();
        int firstTabIndex = line.indexOf('\t');
        String lineValues = line.substring(firstTabIndex + 1);

        String[] readAligments = lineValues.split(Cloudbreak.READ_SEPARATOR);
        String read1AlignmentsString = readAligments[0];
        String[] read1Alignments = read1AlignmentsString.split(Cloudbreak.ALIGNMENT_SEPARATOR);
        List<AlignmentRecord> read1AlignmentRecords = alignmentReader.parseAlignmentsIntoRecords(read1Alignments);

        String read2AlignmentsString = readAligments[1];
        String[] read2Alignments = read2AlignmentsString.split(Cloudbreak.ALIGNMENT_SEPARATOR);
        List<AlignmentRecord> read2AlignmentRecords = alignmentReader.parseAlignmentsIntoRecords(read2Alignments);

        Set<AlignmentRecord> recordsInExcludedAreas = new HashSet<AlignmentRecord>();
        try {
            if (exclusionRegions != null) {
                for (AlignmentRecord record : read1AlignmentRecords) {
                    if (exclusionRegions.doesLocationOverlap(record.getChromosomeName(), record.getPosition(), record.getPosition() + record.getSequence().length())) {
                        recordsInExcludedAreas.add(record);
                    }
                }

                for (AlignmentRecord record : read2AlignmentRecords) {
                    if (exclusionRegions.doesLocationOverlap(record.getChromosomeName(), record.getPosition(), record.getPosition() + record.getSequence().length())) {
                        recordsInExcludedAreas.add(record);
                    }
                }

            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

        try {
            emitReadPairInfoForAllPairs(read1AlignmentRecords, read2AlignmentRecords, output, recordsInExcludedAreas);
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }

    private void emitReadPairInfoForAllPairs(List<AlignmentRecord> read1AlignmentRecords, List<AlignmentRecord> read2AlignmentRecords, OutputCollector<GenomicLocation, ReadPairInfo> output, Set<AlignmentRecord> recordsInExcludedAreas) throws Exception {
        for (AlignmentRecord record1 : read1AlignmentRecords) {
            for (AlignmentRecord record2 : read2AlignmentRecords) {
                if (recordsInExcludedAreas.contains(record1) || recordsInExcludedAreas.contains(record2)) return;
                emitReadPairInfoForPair(record1, record2, output);
            }
        }
    }

    private void emitReadPairInfoForPair(AlignmentRecord record1, AlignmentRecord record2, OutputCollector<GenomicLocation, ReadPairInfo> output) throws IOException {

        // todo: not handling translocations for now
        if (! record1.getChromosomeName().equals(record2.getChromosomeName())) {
            return;
        }

        if (getChromosomeFilter() != null) {
            if (! record1.getChromosomeName().equals(getChromosomeFilter())) return;
        }

        double endPosterior1 = ((NovoalignNativeRecord) record1).getPosteriorProb();
        double endPosterior2 = ((NovoalignNativeRecord) record2).getPosteriorProb();

        int insertSize;
        AlignmentRecord leftRead = record1.getPosition() < record2.getPosition() ?
                record1 : record2;
        AlignmentRecord rightRead = record1.getPosition() < record2.getPosition() ?
                record2 : record1;

        // todo: not handling inversions for now
        if (!scorer.validateMappingOrientations(record1, record2, matePairs)) return;

        insertSize = rightRead.getPosition() + rightRead.getSequence().length() - leftRead.getPosition();

        if (! scorer.validateInsertSize(insertSize, record1.getReadId(), maxInsertSize)) return;

        int genomeOffset = leftRead.getPosition() - leftRead.getPosition() % resolution;


        int genomicWindow = insertSize +
                leftRead.getPosition() % resolution +
                (resolution - rightRead.getPosition() % resolution);


        double pMappingCorrect = scorer.probabilityMappingIsCorrect(endPosterior1, endPosterior2);

        if (mapabilityWeighting != null) {
            if (insertSize > targetIsize + 6 * targetIsizeSD) {
                String chrom = record1.getChromosomeName();
                int leftReadStart = leftRead.getPosition();
                int leftReadEnd = leftRead.getPosition() + leftRead.getSequence().length();
                double leftReadMapability = mapabilityWeighting.getMinValueForRegion(chrom, leftReadStart, leftReadEnd);
                //System.err.println("left read mapability from " + leftRead.getPosition() + " to " + (leftRead.getPosition() + leftRead.getSequence().length()) + " = " + leftReadMapability);

                int rightReadStart = rightRead.getPosition() - rightRead.getSequence().length();
                int rightReadEnd = rightRead.getPosition();
                double rightReadMapability = mapabilityWeighting.getMinValueForRegion(chrom, rightReadStart, rightReadEnd);
                //System.err.println("right read mapability from " + (rightRead.getPosition() - rightRead.getSequence().length()) + " to " + rightRead.getPosition() + " = " + rightReadMapability);

                //System.err.println("old pmc: " + pMappingCorrect);
                pMappingCorrect = pMappingCorrect + Math.log(leftReadMapability) + Math.log(rightReadMapability);
                //System.err.println("new pmc: " + pMappingCorrect);
            }
        }

        ReadPairInfo readPairInfo = new ReadPairInfo(insertSize, pMappingCorrect, readGroupId);

        for (int i = 0; i <= genomicWindow; i = i + resolution) {
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

            //System.err.println("Emitting insert size " + insertSize);
            GenomicLocation genomicLocation = new GenomicLocation(chromosome, pos);
            output.collect(genomicLocation, readPairInfo);

        }

    }


    public void configure(JobConf job) {
        super.configure(job);

        // todo: if we change to the non-deprecated API, need to update this as described here
        // todo: https://issues.apache.org/jira/browse/MAPREDUCE-2166
        String inputFile = getInputPath(job.get("map.input.file"));

        String readGroupInfoFile = job.get("read.group.info.file");
        ReadGroupInfoFileHelper readGroupInfoFileHelper = new ReadGroupInfoFileHelper();
        Map<Short, ReadGroupInfo> readGroupInfos = null;
        try {
            readGroupInfos = readGroupInfoFileHelper.readReadGroupsById(readGroupInfoFile);
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.err.println("Looking up the read group for input file: " + inputFile);
        boolean configuredReadGroup = false;
        for (Short readGroupInfoId : readGroupInfos.keySet()) {
            System.err.println("comparing to: " + readGroupInfos.get(readGroupInfoId).hdfsPath);
            if (inputFile.startsWith(readGroupInfos.get(readGroupInfoId).hdfsPath)) {
                System.err.println("got it!");
                this.readGroupId = readGroupInfoId;
                ReadGroupInfo readGroupInfo = readGroupInfos.get(readGroupInfoId);
                matePairs = readGroupInfo.matePair;
                targetIsize = readGroupInfo.isize;
                targetIsizeSD = readGroupInfo.isizeSD;
                configuredReadGroup = true;
                break;
            }
        }
        if (! configuredReadGroup) throw new RuntimeException("Unable to configure read group for " + inputFile);

        scorer = new ProbabilisticPairedAlignmentScorer();

        faidxFileName = job.get("alignment.faidx");
        faix = new FaidxFileHelper(faidxFileName);

        if (job.get("alignments.filterchr") != null) {
            setChromosomeFilter(job.get("alignments.filterchr"));
            setStartFilter(Long.parseLong(job.get("alignments.filterstart")));
            setEndFilter(Long.parseLong(job.get("alignments.filterend")));
        }

        if (job.get("alignment.exclusionRegions") != null) {
            String exclusionRegionsFileName = job.get("alignment.exclusionRegions");
            try {
                exclusionRegions = new GFFFileHelper(exclusionRegionsFileName);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        if (job.get("alignment.mapabilityWeighting") != null) {
            String mapabilityWeightingFileName = job.get("alignment.mapabilityWeighting");
            mapabilityWeighting = new BigWigFileHelper();
            try {
                mapabilityWeighting.open(mapabilityWeightingFileName);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        if (job.get("pileupDeletionScore.maxInsertSize") != null) {
            maxInsertSize = Integer.parseInt(job.get("pileupDeletionScore.maxInsertSize"));
        }
    }

    protected static String getInputPath(String mapInputProperty) {
        String path = new Path(mapInputProperty).toUri().getPath();
        return path;
    }
}
