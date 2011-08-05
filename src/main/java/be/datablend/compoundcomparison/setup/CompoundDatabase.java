package be.datablend.compoundcomparison.setup;

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
import java.net.UnknownHostException;
import java.util.ArrayList;

/**
 * User: dsuvee
 * Date: 04/08/11
 */
public class CompoundDatabase {

    private DBCollection compoundsCollection = null;
    private DBCollection fingerprintCountsCollection = null;

    public CompoundDatabase() throws UnknownHostException {
        // Init the mongo connection
        Mongo m = new Mongo("localhost", 27017);
		DB db = m.getDB(COMPOUNDS_DATABASE);

        compoundsCollection = db.getCollection(COMPOUNDS_COLLECTION);
        fingerprintCountsCollection = db.getCollection(FINGERPRINTCOUNTS_COLLECTION);
    }

    public void create() throws IOException {
        // Start by creating the required indexes
        createIndexes();
        // Import some compound data
        //importCompounds();
    }

    private void createIndexes() {
        System.out.println("Started setup of compounds database ... ");
        // Ensure indexes to speed up queries
		compoundsCollection.ensureIndex(new BasicDBObject().append(FINGERPRINTS_PROPERTY, 1));
        compoundsCollection.ensureIndex(new BasicDBObject().append(FINGERPRINTCOUNT_PROPERTY, 1));
        fingerprintCountsCollection.ensureIndex(new BasicDBObject().append(FINGERPRINT_PROPERTY, 1));
        System.out.println("Finished setup of compounds database ... ");
    }

    private void importCompounds() throws IOException {
        // Create the compound file reader and fingerprinter
        System.out.println("Started reading the compound input file ... ");
        RandomAccessMDLReader reader = new RandomAccessMDLReader(new File(getClass().getClassLoader().getResource("Compound_046200001_046225000.sdf").getFile()));
        EncodingFingerprint fingerprinter = new Encoding2DMolprint();
        System.out.println("Finished reading the compound input file ... ");

        System.out.println("Started import of " + reader.getSize() + " compounds  ... ");

        // Iterate the compounds one by one
        for (int i = 0; i < reader.getSize(); i++) {

            // Retrieve the molecule and the fingerprints for this molecule
            Molecule molecule = reader.getMol(i);
            FeatureMap fingerprints = new FeatureMap(fingerprinter.getFingerprint(molecule));

            // Retrieve some of the compound properties we want to use later on
            String compound_cid = (String)molecule.getProperty("PUBCHEM_COMPOUND_CID");
            String smiles = (String)molecule.getProperties().get("PUBCHEM_OPENEYE_CAN_SMILES");

            // Create the new document that will hold the compound
            BasicDBObject compound = new BasicDBObject();

            // Add the simple properties one by one to this compound object
            compound.put(COMPOUNDCID_PROPERTY,compound_cid);
            compound.put(SMILES_PROPERTY,smiles);
            compound.put(FINGERPRINTCOUNT_PROPERTY, fingerprints.getSize());

            // Create temporary holder for the fingerprints we need to save
            ArrayList<String> fingerprintstosave = new ArrayList<String>();

            // Iterate the fingerprints
            for (IFeature fingerprint : fingerprints.getKeySet()) {

                // First add it the string representation to the temporary holder
                fingerprintstosave.add(fingerprint.featureToString());

                // Create a +1 increment count for this fingerprint
                // This enables to increment the 'count' property of some document by 1 (or setting it on 1 in case it does not yet exist)
                BasicDBObject countplusone = new BasicDBObject();
                countplusone.put(COUNT_PROPERTY,1);
                BasicDBObject increment = new BasicDBObject();
                increment.put("$inc", countplusone);

                // Create the fingerprint document itself
                BasicDBObject the_fingerprint = new BasicDBObject();
                the_fingerprint.put(FINGERPRINT_PROPERTY, fingerprint.featureToString());

                // Perform an update on the fingerprints counts collection as an upsert (so creating the document itself if it woul not yet exist)
                fingerprintCountsCollection.update(the_fingerprint,increment,true,false);
            }

            // Add the full fingerprint list to the compound document
            compound.put(FINGERPRINTS_PROPERTY, fingerprintstosave.toArray(new String[]{}));

            // Compound document is created, save it.
            compoundsCollection.insert(compound);
        }

        System.out.println("Finished import of " + reader.getSize() + " compounds  ... \n");

    }

    public DBCollection getCompoundsCollection() {
        return compoundsCollection;
    }

    public DBCollection getFingerprintCountsCollection() {
        return fingerprintCountsCollection;
    }

}
