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
    double posteriorProb;
    boolean forward;
    String readId;
    String sequence;

    public static NovoalignNativeRecord parseRecord(String[] fields) {
        NovoalignNativeRecord record = new NovoalignNativeRecord();
        record.setReadId(fields[0]);
        record.setSequence(fields[2]);
        record.setMappingStatus(fields[4]);
        if (record.isMapped()) {
            String recordReferenceName = fields[7];
            // cut off the ">" that starts the chromosome name
            if (recordReferenceName.startsWith(">")) {
                recordReferenceName = recordReferenceName.substring(1);
            }
            record.setReferenceName(recordReferenceName);

            record.setPosition(Integer.parseInt(fields[8]));
            record.setPosteriorProb(Double.parseDouble(fields[6]));
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

    public double getPosteriorProb() {
        return posteriorProb;
    }

    public void setPosteriorProb(double posteriorProb) {
        this.posteriorProb = posteriorProb;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        NovoalignNativeRecord that = (NovoalignNativeRecord) o;

        if (forward != that.forward) return false;
        if (position != that.position) return false;
        if (Double.compare(that.posteriorProb, posteriorProb) != 0) return false;
        if (mappingStatus != null ? !mappingStatus.equals(that.mappingStatus) : that.mappingStatus != null)
            return false;
        if (readId != null ? !readId.equals(that.readId) : that.readId != null) return false;
        if (referenceName != null ? !referenceName.equals(that.referenceName) : that.referenceName != null)
            return false;
        if (sequence != null ? !sequence.equals(that.sequence) : that.sequence != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = referenceName != null ? referenceName.hashCode() : 0;
        result = 31 * result + position;
        result = 31 * result + (mappingStatus != null ? mappingStatus.hashCode() : 0);
        temp = posteriorProb != +0.0d ? Double.doubleToLongBits(posteriorProb) : 0L;
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        result = 31 * result + (forward ? 1 : 0);
        result = 31 * result + (readId != null ? readId.hashCode() : 0);
        result = 31 * result + (sequence != null ? sequence.hashCode() : 0);
        return result;
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

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }
}
