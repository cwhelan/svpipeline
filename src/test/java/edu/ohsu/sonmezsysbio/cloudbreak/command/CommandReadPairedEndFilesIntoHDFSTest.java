package edu.ohsu.sonmezsysbio.cloudbreak.command;

import org.junit.Test;
import static junit.framework.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 7/22/12
 * Time: 1:48 PM
 */
public class CommandReadPairedEndFilesIntoHDFSTest {

    @Test
    public void testGetLongestCommonPrefix() throws Exception {
        String s1 = "@EAS139_44:4:1:2:5501";
        String s2 = "@EAS139_44:4:1:2:5502";
        assertEquals("@EAS139_44:4:1:2:550", CommandReadPairedEndFilesIntoHDFS.greatestCommonPrefix(s1, s2));
    }
}
