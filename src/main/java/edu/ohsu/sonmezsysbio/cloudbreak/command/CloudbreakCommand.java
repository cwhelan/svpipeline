package edu.ohsu.sonmezsysbio.cloudbreak.command;

import org.apache.hadoop.conf.Configuration;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/20/11
 * Time: 1:55 PM
 */
public interface CloudbreakCommand {
    public void run(Configuration conf) throws Exception;
}
