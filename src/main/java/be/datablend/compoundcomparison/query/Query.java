package be.datablend.compoundcomparison.query;

import be.datablend.compoundcomparison.setup.CompoundDatabase;
import com.mongodb.*;
import fingerprinters.EncodingFingerprint;
import fingerprinters.features.FeatureMap;
import fingerprinters.features.IFeature;
import fingerprinters.topological.Encoding2DMolprint;
import io.reader.RandomAccessMDLReader;
import org.openscience.cdk.Molecule;

import static be.datablend.compoundcomparison.definition.Definition.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * User: dsuvee
 * Date: 04/08/11
 */
public class Query {

    private DBCollection compoundsCollection = null;
    private DBCollection fingerprintCountsCollection = null;

    public Query(CompoundDatabase compoundDatabase) {
         compoundsCollection = compoundDatabase.getCompoundsCollection();
         fingerprintCountsCollection = compoundDatabase.getFingerprintCountsCollection();
    }

    public void findSimilarCompounds(int compoundIndexInFile, double similarity) throws IOException {
        // Create the compound file reader and fingerprinter
        System.out.println("Started reading the compound input file ... ");
        RandomAccessMDLReader reader = new RandomAccessMDLReader(new File(getClass().getClassLoader().getResource("Compound_046200001_046225000.sdf").getFile()));
        EncodingFingerprint fingerprinter = new Encoding2DMolprint();
        System.out.println("Finished reading the compound input file ... \n");
        Molecule molecule = reader.getMol(compoundIndexInFile);
        FeatureMap fingerprints = new FeatureMap(fingerprinter.getFingerprint(molecule));

        // Retrieve the relevant properties
        String pubchemcid = (String)molecule.getProperties().get("PUBCHEM_COMPOUND_CID");
        String smiles = (String)molecule.getProperties().get("PUBCHEM_OPENEYE_CAN_SMILES");
        System.out.println("Trying to find " + (similarity * 100) + "% similar molecues for PubChem compound " + pubchemcid  + " ( " + smiles + " )");

        // Extract the fingerprints to find
        List<String> fingerprintstofind = new ArrayList<String>();
        for (IFeature feature : fingerprints.getKeySet()) {
            fingerprintstofind.add(feature.featureToString());
        }

        // Sort the fingeprints on total number of occurences
        fingerprintstofind = findSortedFingerprints(fingerprintstofind);

        // Execute the query using the build-in mongodb query functionalities
        long startnative = System.currentTimeMillis();
        executeNativeQuery(similarity, fingerprintstofind);
        long stopnative = System.currentTimeMillis();
        System.out.println("Total time for native query: " + (stopnative-startnative) + " ms\n");
        long startmr = System.currentTimeMillis();
        executeMapReduceQuery(similarity, fingerprintstofind);
        long stopmr = System.currentTimeMillis();
        System.out.println("Total time for map reduce query: " + (stopmr-startmr) + " ms");

    }

    private void executeNativeQuery(double similarity, List<String> fingerprintsToFind) {
        System.out.println("Executing using native query ...");

        // Calculate the essential numbers
        int maxnumberofcompoundfingerprints = (int) (fingerprintsToFind.size() / similarity);
        int minnumberofcompoundfingerprints = (int) (fingerprintsToFind.size() * similarity);
        int numberoffingerprintstoconsider = fingerprintsToFind.size() - minnumberofcompoundfingerprints;

        List<String> fingerprintsToConsider = fingerprintsToFind.subList(0,numberoffingerprintstoconsider+1);

        // Find all compounds that have one (or more) fingerprints in the list that we consider, but for which the total count of fingeprints is within the defined limits
        DBObject compoundquery = QueryBuilder.start(FINGERPRINTS_PROPERTY).in(fingerprintsToConsider).and(FINGERPRINTCOUNT_PROPERTY).lessThanEquals(maxnumberofcompoundfingerprints).and(FINGERPRINTCOUNT_PROPERTY).greaterThanEquals(minnumberofcompoundfingerprints).get();

        // Execute the query
        DBCursor compounds = compoundsCollection.find(compoundquery);

        // Let's process the found compounds locally
        while(compounds.hasNext()) {

            DBObject compound = compounds.next();
            // Retrieve all fingerprints
            BasicDBList fingerprints = ((BasicDBList) compound.get(FINGERPRINTS_PROPERTY));
            // Calculate the intersection on the total list of fingerprints
            fingerprints.retainAll(fingerprintsToFind);

            // If the remaining list of fingeprints contains at least the minimun number of featues
            if (fingerprints.size() >= minnumberofcompoundfingerprints) {

                // Retrieve the total count
                int totalcount = (Integer)compound.get(FINGERPRINTCOUNT_PROPERTY);
                // Calculate the tanimoto coefficient
                double tanimoto = (double) fingerprints.size() / (totalcount + fingerprintsToFind.size() - fingerprints.size());
                // Although we reduced the search space, we still need to check whether the tanimoto is really >= the required similarity
                if (tanimoto >= similarity) {
                    // Print it
                    System.out.println("Found compound with PubChemd ID " + compound.get(COMPOUNDCID_PROPERTY) + " ( " + compound.get(SMILES_PROPERTY) + " )  with similarity score of " + (int)(tanimoto * 100) + "%");
                }
            }
        }

    }

    private void executeMapReduceQuery(double similarity, List<String> fingerprintsToFind) {
        System.out.println("Executing using map reduce query ...");

        // Calculate the essential numbers
        int maxnumberofcompoundfingerprints = (int) (fingerprintsToFind.size() / similarity);
        int minnumberofcompoundfingerprints = (int) (fingerprintsToFind.size() * similarity);
        int numberoffingerprintstoconsider = fingerprintsToFind.size() - minnumberofcompoundfingerprints;


        List<String> fingerprintsToConsider = fingerprintsToFind.subList(0,numberoffingerprintstoconsider+1);

        // Find all compounds that have one (or more) fingerprints in the list that we consider, but for which the total count of fingeprints is within the defined limits
        DBObject compoundquery = QueryBuilder.start(FINGERPRINTS_PROPERTY).in(fingerprintsToConsider).and(FINGERPRINTCOUNT_PROPERTY).lessThanEquals(maxnumberofcompoundfingerprints).and(FINGERPRINTCOUNT_PROPERTY).greaterThanEquals(minnumberofcompoundfingerprints).get();

        // The map fuction
        String map = "function() {  " +
                        "var found = 0; " +
                        "var fingerprintslength = this.fingerprints.length; " +
                        "for (i = 0; i < fingerprintslength; i++) { " +
                            "if (fingerprintstofind[this.fingerprints[i]] === true) { found++; } " +
                        "} " +
                        "if (found >= minnumberofcompoundfingerprints) { emit (this.compound_cid, {found : found, total : this.fingerprint_count, smiles: this.smiles} ); } " +
                     "}";

        // Execute the map reduce function
        MapReduceCommand mr = new MapReduceCommand(compoundsCollection, map, "", null, MapReduceCommand.OutputType.INLINE, compoundquery);

        // Create a hashmap for the fingerprints to find (to speed up the javascript execution)
        Map<String,Boolean> tofind = new HashMap<String,Boolean>();
        for(String fingerprinttofind : fingerprintsToFind) {
            tofind.put(fingerprinttofind,true);
        }

        // Set the map reduce scope
        Map<String,Object> scope = new HashMap<String,Object>();
        scope.put("fingerprintstofind",tofind);
        scope.put("minnumberofcompoundfingerprints",minnumberofcompoundfingerprints);
        mr.setScope(scope);

        // Execute the map reduce
        MapReduceOutput out = compoundsCollection.mapReduce(mr);

        // Iterate the results
        for (DBObject result : out.results()) {
            String compound_cid = (String)result.get("_id");
            DBObject value = (DBObject)result.get("value");

            // Calculate the tanimoto coefficient
            double totalcount = (Double)value.get("total");
            double found = (Double)value.get("found");
            double tanimoto = (Double)value.get("found") / ((Double)value.get("total") + fingerprintsToFind.size() - (Double)value.get("found"));
            // Although we reduced the search space, we still need to check whether the tanimoto is really >= the required similarity
            if (tanimoto >= similarity) {
                // Print it
                System.out.println("Found compound with PubChemd ID " + compound_cid + " ( " + value.get("smiles") + " )  with similarity score of " + (int)(tanimoto * 100) + "%");
            }
        }

    }

    // Helper method to retrieve the fingeprints, but sorted on total count over the entire compound population
    private List<String> findSortedFingerprints(List<String> fingerprintsToFind) {
        System.out.println("Compound has " + fingerprintsToFind.size() + " unique fingerprints\n");

        List<String> sortedFingerprintsToFind = new ArrayList<String>();

        // Find all fingerprint count documents that have a fingerpint in the list of fingerprints to find
        DBObject fingerprintcountquery = QueryBuilder.start(FINGERPRINT_PROPERTY).in(fingerprintsToFind.toArray()).get();
        // Only retrieve the fingerprint string itself
        DBObject fingerprintcountselection = QueryBuilder.start(FINGERPRINT_PROPERTY).is(1).get();
        // Sort the result on count
        DBObject fingerprintcountsort = QueryBuilder.start(COUNT_PROPERTY).is(1).get();

        // Execute the query on the fingerprint counts collection
        DBCursor fingerprintcounts = fingerprintCountsCollection.find(fingerprintcountquery, fingerprintcountselection).sort(fingerprintcountsort);

        // Iterate and add them one by one (fingerprints with the smallest count will come first)
        while (fingerprintcounts.hasNext()) {
            DBObject fingerprintcount = fingerprintcounts.next();
            sortedFingerprintsToFind.add((String)fingerprintcount.get(FINGERPRINT_PROPERTY));
        }

        return sortedFingerprintsToFind;

    }

}
