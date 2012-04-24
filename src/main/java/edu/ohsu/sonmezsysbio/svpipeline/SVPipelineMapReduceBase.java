package edu.ohsu.sonmezsysbio.svpipeline;

import org.apache.hadoop.mapred.JobConf;
import org.apache.hadoop.mapred.MapReduceBase;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/23/12
 * Time: 10:05 PM
 */
public class SVPipelineMapReduceBase extends MapReduceBase {


    protected int resolution = SVPipeline.DEFAULT_RESOLUTION;

    public int getResolution() {
        return resolution;
    }

    public void setResolution(int resolution) {
        this.resolution = resolution;
    }

    @Override
    public void configure(JobConf job) {
        super.configure(job);
        if (job.get("svpipeline.resolution") != null) {
            resolution = Integer.parseInt(job.get("svpipeline.resolution"));
        }
    }
}
