package edu.ohsu.sonmezsysbio.svpipeline.io;

import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.WritableComparable;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/6/12
 * Time: 1:05 PM
 */
public class GenomicLocation implements WritableComparable<GenomicLocation> {
    public Text chromosome;
    public int pos;

    public GenomicLocation(Text chromosome, int pos) {
        this.chromosome = chromosome;
        this.pos = pos;
    }

    public void write(DataOutput out) throws IOException {
        chromosome.write(out);
        out.writeInt(pos);
    }

    public void readFields(DataInput in) throws IOException {
        chromosome.readFields(in);
        pos = in.readInt();
    }

    public int compareTo(GenomicLocation o) {
        if (chromosome.toString().compareTo(o.chromosome.toString()) != 0) {
            return chromosome.toString().compareTo(o.chromosome.toString());
        } else {
            if (pos < o.pos) return -1;
            if (pos > o.pos) return 1;
        }
        return 0;
    }
}
