package edu.ohsu.sonmezsysbio.cloudbreak;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 8/8/12
 * Time: 3:18 PM
 */
public class ReadPairAlignments {
    List<AlignmentRecord> read1Alignments;
    List<AlignmentRecord> read2Alignments;

    public ReadPairAlignments(List<AlignmentRecord> read1Alignments, List<AlignmentRecord> read2Alignments) {
        this.read1Alignments = read1Alignments;
        this.read2Alignments = read2Alignments;
    }

    public List<AlignmentRecord> getRead1Alignments() {
        return read1Alignments;
    }

    public void setRead1Alignments(List<AlignmentRecord> read1Alignments) {
        this.read1Alignments = read1Alignments;
    }

    public List<AlignmentRecord> getRead2Alignments() {
        return read2Alignments;
    }

    public void setRead2Alignments(List<AlignmentRecord> read2Alignments) {
        this.read2Alignments = read2Alignments;
    }
}
