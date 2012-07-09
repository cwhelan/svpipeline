package edu.ohsu.sonmezsysbio.svpipeline.file;

import edu.ohsu.sonmezsysbio.svpipeline.ReadGroupInfo;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 7/9/12
 * Time: 10:03 AM
 */
public class ReadGroupInfoFileHelper {

    public Map<String,Short> readReadGroupIdsByName(String filename) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        try {
            return readReadGroupIdsByName(reader);
        } finally {
            reader.close();
        }
    }

    public Map<String,Short> readReadGroupIdsByName(BufferedReader reader) throws IOException {
         String line;
         short rgIdx = 0;
         Map<String,Short> rgsByName = new HashMap<String,Short>();
         while ((line = reader.readLine()) != null) {
             String rgName = line.split("\t")[0];
             rgsByName.put(rgName, rgIdx);
             rgIdx++;
         }
         return rgsByName;
    }

    public Map<Short, ReadGroupInfo> readReadGroupsById(String filename) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        try {
            return readReadGroupsById(reader);
        }finally {
            reader.close();
        }
    }

    public Map<Short, ReadGroupInfo> readReadGroupsById(BufferedReader reader) throws IOException {
        String line;
        Map<Short,ReadGroupInfo> rgsById = new HashMap<Short,ReadGroupInfo>();
        short rgIdx = 0;
        while ((line = reader.readLine()) != null) {
            String[] fields = line.split("\\t");
            ReadGroupInfo readGroupInfo = new ReadGroupInfo();
            readGroupInfo.readGroupName = fields[0];
            readGroupInfo.libraryName = fields[1];
            readGroupInfo.isize = Integer.parseInt(fields[2]);
            readGroupInfo.isizeSD = Integer.parseInt(fields[3]);
            readGroupInfo.matePair = Boolean.parseBoolean(fields[4]);
            rgsById.put(rgIdx, readGroupInfo);
            rgIdx++;
        }
        return rgsById;
    }
}
