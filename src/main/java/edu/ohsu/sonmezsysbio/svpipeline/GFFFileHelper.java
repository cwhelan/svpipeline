package edu.ohsu.sonmezsysbio.svpipeline;

import org.biojava3.genome.parsers.gff.FeatureList;
import org.biojava3.genome.parsers.gff.GFF3Reader;
import org.biojava3.genome.parsers.gff.Location;

import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 4/19/12
 * Time: 1:07 PM
 */
public class GFFFileHelper {

    FeatureList features;

    public GFFFileHelper() {

    }

    public GFFFileHelper(String filename) throws IOException {
        features = GFF3Reader.read(filename);
    }

    public boolean doesLocationOverlap(String chrom, int start, int end) throws Exception {
        Location query = new Location(start, end);
        return features.selectOverlapping(chrom, query, true).size() > 0;
    }
}
