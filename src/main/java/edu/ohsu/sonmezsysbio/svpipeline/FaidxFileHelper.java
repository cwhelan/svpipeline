package edu.ohsu.sonmezsysbio.svpipeline;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/9/12
 * Time: 4:02 PM
 */
public class FaidxFileHelper {

    Map<Short, String> chromNamesByKey;
    Map<String, Short> chromKeysByName;
    private String faidxFileName;

    public FaidxFileHelper(String faidxFileName) {
        this.faidxFileName = faidxFileName;
    }

    public String getNameForChromKey(Short key) throws IOException {
        if (chromNamesByKey == null) {
            chromNamesByKey = readChromNamesByKey(faidxFileName);
        }
        return chromNamesByKey.get(key);
    }

    public Short getKeyForChromName(String name) throws IOException {
        if (chromKeysByName == null) {
            chromKeysByName = readChromKeysByName(faidxFileName);
        }
        return chromKeysByName.get(name);
    }

    public Map<Short, String> readChromNamesByKey(String faidxFileName) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(faidxFileName));
        return readChromNamesByKey(reader);
    }

    public Map<Short, String> readChromNamesByKey(BufferedReader bufferedReader) throws IOException {
        String line;
        short chrIdx = 0;
        Map<Short, String> chrTable = new HashMap<Short, String>();
        while ((line = bufferedReader.readLine()) != null) {
            String chromosomeName = line.split("\\s+")[0];
            chrTable.put(chrIdx, chromosomeName);
            chrIdx++;
        }
        return chrTable;
    }

    public Map<String, Short> readChromKeysByName(String faidxFileName) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(faidxFileName));
        return readChromKeysByName(reader);
    }

    public Map<String, Short> readChromKeysByName(BufferedReader bufferedReader) throws IOException {
        String line;
        short chrIdx = 0;
        Map<String, Short> chrTable = new HashMap<String, Short>();
        while ((line = bufferedReader.readLine()) != null) {
            String chromosomeName = line.split("\\s+")[0];
            chrTable.put(chromosomeName, chrIdx);
            chrIdx++;
        }
        return chrTable;
    }

}
