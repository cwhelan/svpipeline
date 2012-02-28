package edu.ohsu.sonmezsysbio.svpipeline.command;

import org.apache.hadoop.conf.Configuration;

import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/20/11
 * Time: 1:55 PM
 */
public interface SVPipelineCommand {
    public void run(Configuration conf) throws Exception;
}
