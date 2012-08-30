package edu.ohsu.sonmezsysbio.cloudbreak.io;

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
public class GenomicLocationWithQuality implements WritableComparable<GenomicLocationWithQuality> {
    public short chromosome;
    public int pos;
    public double quality;

    public GenomicLocationWithQuality() {
    }

    public GenomicLocationWithQuality(short chromosome, int pos, double quality) {
        this.chromosome = chromosome;
        this.pos = pos;
        this.quality = quality;
    }

    public void write(DataOutput out) throws IOException {
        out.writeShort(chromosome);
        out.writeInt(pos);
        out.writeDouble(quality);
    }

    public void readFields(DataInput in) throws IOException {
        chromosome = in.readShort();
        pos = in.readInt();
        quality = in.readDouble();
    }

    public int compareTo(GenomicLocationWithQuality o) {
        if (chromosome < o.chromosome) {
            return -1;
        } else if (chromosome > o.chromosome) {
            return 1;
        } else {
            if (pos < o.pos) return -1;
            if (pos > o.pos) return 1;
            if (quality > o.quality) return 1;
            if (quality < o.quality) return -1;
        }
        return 0;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        GenomicLocationWithQuality that = (GenomicLocationWithQuality) o;

        if (chromosome != that.chromosome) return false;
        if (pos != that.pos) return false;
        if (quality != that.quality) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = (int) chromosome;
        result = 31 * result + pos;
        temp = quality != +0.0d ? Double.doubleToLongBits(quality) : 0L;
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }

    @Override
    public String toString() {
        return "GenomicLocation{" +
                "chromosome=" + chromosome +
                ", pos=" + pos +
                ", quality=" + quality +
                '}';
    }
}
