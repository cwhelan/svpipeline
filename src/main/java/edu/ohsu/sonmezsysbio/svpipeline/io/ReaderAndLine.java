package edu.ohsu.sonmezsysbio.svpipeline.io;

import org.apache.hadoop.io.SequenceFile;

import java.io.IOException;

/**
* Created by IntelliJ IDEA.
* User: cwhelan
* Date: 4/10/12
* Time: 11:43 AM
*/
public interface ReaderAndLine extends Comparable<ReaderAndLine> {

    GenomicLocation getGenomicLocation();

    void setGenomicLocation(GenomicLocation genomicLocation);

    Double getNextValue();

    void setNextValue(Double nextValue);

    void closeInput() throws IOException;
}
