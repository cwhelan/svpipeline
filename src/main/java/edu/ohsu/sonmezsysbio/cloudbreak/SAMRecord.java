package edu.ohsu.sonmezsysbio.cloudbreak;

import java.util.HashMap;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 5/23/11
 * Time: 10:44 AM
 */
public class SAMRecord {

    int flag;
    Map<String,String> tags = new HashMap<String, String>();
    String referenceName;
    String pairReferenceName;
    int position;
    int insertSize;

    public static SAMRecord parseSamRecord(String[] fields) {
        SAMRecord samRecord = new SAMRecord();
        samRecord.setFlag(Integer.parseInt(fields[1]));
        samRecord.setReferenceName(fields[2]);
        samRecord.setPosition(Integer.parseInt(fields[3]));
        if (! samRecord.isMapped()) {
            return samRecord;
        }

        samRecord.setPairReferenceName(fields[6]);
        samRecord.setInsertSize(Integer.parseInt(fields[8]));

        int optionalFieldsStart = 11;
        for (int i = optionalFieldsStart; i < fields.length; i++) {
            samRecord.addTag(fields[i]);
        }
        return samRecord;
    }

    public int getPosition() {
        return position;
    }

    public void setPosition(int position) {
        this.position = position;
    }

    public int getInsertSize() {
        return insertSize;
    }

    public void setInsertSize(int insertSize) {
        this.insertSize = insertSize;
    }

    public String getReferenceName() {
        return referenceName;
    }

    public void setReferenceName(String referenceName) {
        this.referenceName = referenceName;
    }

    public String getPairReferenceName() {
        return pairReferenceName;
    }

    public void setPairReferenceName(String pairReferenceName) {
        this.pairReferenceName = pairReferenceName;
    }

    public int getFlag() {
        return flag;
    }

    public void setFlag(int flag) {
        this.flag = flag;
    }

    public void addTag(String tag) {
        String[] tagFields = tag.split(":");
        tags.put(tagFields[0], tagFields[2]);
    }

    public boolean isMapped() {
        return ! ((flag & 0x4) == 0x4);
    }

    public boolean isMateMapped() {
        return ! ((flag & 0x8) == 0x8);
    }

    public int getPairPosterior() {
        return Integer.parseInt(tags.get("PQ"));
    }

    public int getEndPosterior() {
        return Integer.parseInt(tags.get("UQ"));
    }

    public boolean isPairMapped() {
        return isMapped() && isMateMapped();
    }

    public boolean isReverseComplemented() {
        return ! ((flag & 0x10) == 0x10);
    }

    public boolean isInterChromosomal() {
        return ! "=".equals(pairReferenceName) && ! referenceName.equals(pairReferenceName);
    }

    public boolean isProperPair() {
        return ((flag & 0x2) == 0x2);
    }


    public boolean isAlignmentOfFirstRead() {
        return ((flag & 0x40) == 0x40);
    }
}
