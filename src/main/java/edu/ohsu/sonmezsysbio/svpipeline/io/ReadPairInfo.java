package edu.ohsu.sonmezsysbio.svpipeline.io;


import org.apache.hadoop.io.Writable;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/6/12
 * Time: 12:58 PM
 */
public class ReadPairInfo implements Writable {
    public int insertSize;
    public double pMappingCorrect;

    public ReadPairInfo() {
    }

    public ReadPairInfo(int insertSize, double pMappingCorrect) {
        this.insertSize = insertSize;
        this.pMappingCorrect = pMappingCorrect;
    }

    public void write(DataOutput out) throws IOException {
        out.writeLong(insertSize);
        out.writeDouble(pMappingCorrect);
    }

    public void readFields(DataInput in) throws IOException {
        insertSize = in.readInt();
        pMappingCorrect = in.readDouble();
    }

    @Override
    public String toString() {
        return "ReadPairInfo{" +
                "insertSize=" + insertSize +
                ", pMappingCorrect=" + pMappingCorrect +
                '}';
    }
}
