package edu.ohsu.sonmezsysbio.svpipeline;

import javax.xml.soap.Text;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 6/5/11
 * Time: 2:26 PM
 */
public class NovoalignNativeRecord {

    String referenceName;
    int position;
    String mappingStatus;
    int posteriorProb;
    boolean forward;
    String readId;

    public static NovoalignNativeRecord parseRecord(String[] fields) {
        NovoalignNativeRecord record = new NovoalignNativeRecord();
        record.setReadId(fields[0]);
        record.setMappingStatus(fields[4]);
        if (record.isMapped()) {
            String recordReferenceName = fields[7];
            // cut off the ">" that starts the chromosome name
            if (recordReferenceName.startsWith(">")) {
                recordReferenceName = recordReferenceName.substring(1);
            }
            record.setReferenceName(recordReferenceName);

            record.setPosition(Integer.parseInt(fields[8]));
            record.setPosteriorProb(Integer.parseInt(fields[6]));
            record.setForward("F".equals(fields[9]));
        }

        return record;

    }

    public boolean isMapped() {
        return "U".equals(getMappingStatus()) || "R".equals(getMappingStatus());
    }

    public String getChromosomeName() {
        return referenceName;
    }

    public void setReferenceName(String referenceName) {
        this.referenceName = referenceName;
    }

    public int getPosition() {
        return position;
    }

    public void setPosition(int position) {
        this.position = position;
    }

    public String getMappingStatus() {
        return mappingStatus;
    }

    public void setMappingStatus(String mappingStatus) {
        this.mappingStatus = mappingStatus;
    }

    public int getPosteriorProb() {
        return posteriorProb;
    }

    public void setPosteriorProb(int posteriorProb) {
        this.posteriorProb = posteriorProb;
    }

    public boolean isForward() {
        return forward;
    }

    public void setForward(boolean forward) {
        this.forward = forward;
    }

    public String getReadId() {
        return readId;
    }

    public void setReadId(String readId) {
        this.readId = readId;
    }
}
