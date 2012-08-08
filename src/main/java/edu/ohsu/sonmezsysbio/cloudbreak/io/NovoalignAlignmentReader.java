package edu.ohsu.sonmezsysbio.cloudbreak.io;

import edu.ohsu.sonmezsysbio.cloudbreak.AlignmentRecord;
import edu.ohsu.sonmezsysbio.cloudbreak.NovoalignNativeRecord;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 8/8/12
 * Time: 1:56 PM
 */
public class NovoalignAlignmentReader extends BaseAlignmentReader {
    public AlignmentRecord parseRecord(String[] fields) {
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
            record.setChromosomeName(recordReferenceName);

            record.setPosition(Integer.parseInt(fields[8]));
            record.setPosteriorProb(Double.parseDouble(fields[6]));
            record.setForward("F".equals(fields[9]));
        }

        return record;

    }

}
