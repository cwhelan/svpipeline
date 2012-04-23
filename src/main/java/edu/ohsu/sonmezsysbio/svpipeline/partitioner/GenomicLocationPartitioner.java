package edu.ohsu.sonmezsysbio.svpipeline.partitioner;

import edu.ohsu.sonmezsysbio.svpipeline.io.GenomicLocation;
import edu.ohsu.sonmezsysbio.svpipeline.io.ReadPairInfo;
import org.apache.hadoop.mapred.MapReduceBase;
import org.apache.hadoop.mapred.Partitioner;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/22/12
 * Time: 10:01 PM
 */
public class GenomicLocationPartitioner extends MapReduceBase implements Partitioner<GenomicLocation, ReadPairInfo> {
    public int getPartition(GenomicLocation key, ReadPairInfo value, int numPartitions) {
        return key.pos % numPartitions;
    }
}
