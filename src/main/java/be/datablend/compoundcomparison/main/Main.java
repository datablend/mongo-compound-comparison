package be.datablend.compoundcomparison.main;

import be.datablend.compoundcomparison.query.Query;
import be.datablend.compoundcomparison.setup.CompoundDatabase;

import java.io.IOException;

/**
 * User: dsuvee
 * Date: 06/08/11
 */
public class Main {

    public static void main(String[] args) throws IOException {

        // Start by creating the compound database
        CompoundDatabase compounddatabase= new CompoundDatabase();
        compounddatabase.create();

        // Execute compound similarty queries
        Query query = new Query(compounddatabase);
        // Find compounds with 80% similarity for the 9000 compound in the input file
        query.findSimilarCompounds(9000,0.8);



    }

}
