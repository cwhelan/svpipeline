package edu.ohsu.sonmezsysbio.cloudbreak;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 8/8/12
 * Time: 1:46 PM
 */
public interface AlignmentRecord {
    boolean isMapped();

    String getChromosomeName();

    void setChromosomeName(String chromosomeName);

    int getPosition();

    void setPosition(int position);

    boolean isForward();

    void setForward(boolean forward);

    String getReadId();

    void setReadId(String readId);

    int getSequenceLength();

}
