package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.*;

import java.io.PrintStream;
import java.io.FileNotFoundException;
import java.util.*;

import net.sf.samtools.SAMRecord;

/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 3, 2009
 *
 * This walker is designed to work as the first pass in a two-pass processing step.
 * It does a by-locus traversal operating only at sites that are not in dbSNP.
 * We assume that all mismatches we see are therefore errors.
 * This walker generates tables based on various user-specified covariates (such as read group, reported quality score, cycle, and dinuc)
 * Since there is a large amount of data one can then calculate an empirical probability of error
 *   given the particular covariates seen at this site, where p(error) = num mismatches / num observations
 * The output file is a CSV list of (several covariate values, num observations, num mismatches, empirical quality score)
 * The first lines of the output file give the name of the covariate classes that were used for this calculation.
 *
 * Note: ReadGroupCovariate and QualityScoreCovariate are required covariates and must be at the start of the list.
 * Note: This walker is designed to be used in conjunction with TableRecalibrationWalker.
 */

@WalkerName("CountCovariatesRefactored")
public class CovariateCounterWalker extends LocusWalker<Integer, PrintStream> {

    @Argument(fullName="list", shortName="ls", doc="List the available covariates and exit", required=false)
    private Boolean LIST_ONLY = false;
    @Argument(fullName="covariate", shortName="cov", doc="Covariates to be used in the recalibration. Each covariate is given as a separate cov parameter. ReadGroup and ReportedQuality are already added for you.", required=false)
    private String[] COVARIATES = null;
    @Argument(fullName="min_mapping_quality", shortName="minmap", required=false, doc="Only use reads with at least this mapping quality score")
    private int MIN_MAPPING_QUALITY = 1;
    @Argument(fullName = "use_original_quals", shortName="OQ", doc="If provided, we will use use the quals from the original qualities OQ attribute field instead of the quals in the regular QUALS field", required=false)
    private boolean USE_ORIGINAL_QUALS = false;
    @Argument(fullName = "platform", shortName="pl", doc="Which sequencing technology was used? This is important for the cycle covariate. Options are SLX, 454, and SOLID.", required=false)
    private String PLATFORM = "SLX";
    @Argument(fullName = "windowSizeNQS", shortName="nqs", doc="How big of a window should the MinimumNQSCovariate use for its calculation", required=false)
    private int WINDOW_SIZE = 3;
    @Argument(fullName="recal_file", shortName="rf", required=false, doc="Filename for the outputted covariates table recalibration file")
    private String RECAL_FILE = "output.recal_data.csv";
    @Argument(fullName="noPrintHeader", shortName="noHeader", required=false, doc="Don't print the usual header on the table recalibration file")
    private boolean NO_PRINT_HEADER = false;
    @Argument(fullName="validateOldRecalibrator", shortName="valor", required=false, doc="Depricated, no longer replicates previous behavior exactly due to clean up of structure of code. If yes will reorder the output to match the old recalibrator exactly but makes the file useless to the refactored TableRecalibrationWalker")
    private boolean validateOldRecalibrator = false;

    private RecalDataManager dataManager; // Holds the data HashMap, mostly used by TableRecalibrationWalker to create collapsed data hashmaps
    private ArrayList<Covariate> requestedCovariates; // A list to hold the covariate objects that were requested
    //private HashMap<SAMRecord, String> readGroupHashMap; // A hash map that hashes the read object itself into the read group name
                                                           // This is done for optimization purposes because pulling the read group out of the SAMRecord is expensive
    private long countedSites = 0; // Number of loci used in the calculations, used for reporting in the output file
    private long countedBases = 0; // Number of bases used in the calculations, used for reporting in the output file
    private long skippedSites = 0; // Number of loci skipped because it was a dbSNP site, used for reporting in the output file

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * Parse the -cov arguments and create a list of covariates to be used here
     * Based on the covariates' estimates for initial capacity allocate the data hashmap
     */
    public void initialize() {

        // Get a list of all available covariates
        final List<Class<? extends Covariate>> classes = PackageUtils.getClassesImplementingInterface(Covariate.class);
        //int estimatedCapacity = 1; // start at one because capacity is multiplicative for each covariate

        // Print and exit if that's what was requested
        if ( LIST_ONLY ) {
            out.println( "Available covariates:" );
            for( Class<?> covClass : classes ) {
                out.println( covClass.getSimpleName() );
            }
            out.println();
            
            System.exit( 0 ); // early exit here because user requested it
        }

        // Initialize the requested covariates by parsing the -cov argument
        requestedCovariates = new ArrayList<Covariate>();
        if( validateOldRecalibrator ) {
            requestedCovariates.add( new ReadGroupCovariate() );
            requestedCovariates.add( new CycleCovariate( PLATFORM ) ); // unfortunately the ordering here is different to match the output of the old recalibrator
            requestedCovariates.add( new QualityScoreCovariate() );
            requestedCovariates.add( new DinucCovariate() );
            //estimatedCapacity = 300 * 200 * 40 * 16;
        } else if( COVARIATES != null ) {
            int covNumber = 1;
            for( String requestedCovariateString : COVARIATES ) {
                boolean foundClass = false;
                for( Class<?> covClass : classes ) {
                    if( requestedCovariateString.equalsIgnoreCase( covClass.getSimpleName() ) ) { // -cov argument matches the class name for an implementing class
                        foundClass = true;
                        // Read Group Covariate and Quality Score Covariate are required covariates for the recalibration calculation and must begin the list
                        if( (covNumber == 1 && !requestedCovariateString.equalsIgnoreCase( "ReadGroupCovariate" )) ||
                            (covNumber == 2 && !requestedCovariateString.equalsIgnoreCase( "QualityScoreCovariate" )) ) {
                            throw new StingException("ReadGroupCovariate and QualityScoreCovariate are required covariates for the recalibration calculation and must begin the list" );
                        }
                        covNumber++;
                        try {
                            // Now that we've found a matching class, try to instantiate it
                            Covariate covariate = (Covariate)covClass.newInstance();
                            // some covariates need parameters (user supplied command line arguments) passed to them
                            if( covariate instanceof CycleCovariate ) { covariate = new CycleCovariate( PLATFORM ); }
                            else if( covariate instanceof PrimerRoundCovariate ) { covariate = new PrimerRoundCovariate( PLATFORM ); }
                            else if( covariate instanceof MinimumNQSCovariate ) { covariate = new MinimumNQSCovariate( WINDOW_SIZE ); }
                            requestedCovariates.add( covariate );
                            //estimatedCapacity *= covariate.estimatedNumberOfBins(); // update the estimated initial capacity
                        } catch ( InstantiationException e ) {
                            throw new StingException( String.format("Can not instantiate covariate class '%s': must be concrete class.", covClass.getSimpleName()) );
                        } catch ( IllegalAccessException e ) {
                            throw new StingException( String.format("Can not instantiate covariate class '%s': must have no-arg constructor.", covClass.getSimpleName()) );
                        }
                    }
                }

                if( !foundClass ) {
                    throw new StingException( "The requested covariate type (" + requestedCovariateString + ") isn't a valid covariate option. Use --list to see possible covariates." );
                }
            }
        } else { // Not validating old recalibrator and no covariates were specified by the user so add the default ones
            logger.info( "Note: Using default set of covariates because none were specified." );
            requestedCovariates.add( new ReadGroupCovariate() );
            requestedCovariates.add( new QualityScoreCovariate() );
            requestedCovariates.add( new CycleCovariate( PLATFORM ) );
            requestedCovariates.add( new DinucCovariate() );
            //estimatedCapacity = 300 * 40 * 200 * 16;
        }

        logger.info( "The covariates being used here: " );
        logger.info( requestedCovariates );

        //dataManager = new RecalDataManager( estimatedCapacity );
        dataManager = new RecalDataManager();
        //readGroupHashMap = new HashMap<SAMRecord, String>( 50000000, 0.97f );
    }


    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * For each read at this locus get the various covariate values and increment that location in the map based on
     *   whether or not the base matches the reference at this particular location
     * @param tracker The reference metadata tracker
     * @param ref The reference context
     * @param context The alignment context
     * @return Returns 1, but this value isn't used in the reduce step
     */
    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        final rodDbSNP dbsnp = rodDbSNP.getFirstRealSNP(tracker.getTrackData("dbsnp", null));

        // Only use data from non-dbsnp sites
        // Assume every mismatch at a non-dbsnp site is indicitive of poor quality
        if( dbsnp == null ) {
            final List<SAMRecord> reads = context.getReads();
            final List<Integer> offsets = context.getOffsets();
            SAMRecord read;
            int offset;
            String readGroup;
            byte[] quals;
            char[] bases;
            char refBase;
            char prevBase;
            
            final int numReads = reads.size();
            // For each read at this locus
            for( int iii = 0; iii < numReads; iii++ ) {
                read = reads.get(iii);
                
                // Only use data from reads with mapping quality above specified quality value
                if( read.getMappingQuality() >= MIN_MAPPING_QUALITY ) {
                
                	//readGroup = readGroupHashMap.get( read );
                    //if( readGroup == null ) { // read is not in the hashmap so add it
                    //    readGroup = read.getReadGroup().getReadGroupId();
                    //    readGroupHashMap.put( read, readGroup );
                    //}
                	
                	offset = offsets.get(iii); // offset is zero based so quals[offset] and bases[offset] is correct
                    
                	// skip first and last base because there is no dinuc, this is mainly done for speed so we don't have to check cases
                	if( offset > 0 && offset < read.getReadLength() - 1) {
	                	
                		quals = read.getBaseQualities();
                        // Check if we need to use the original quality scores instead
                        if ( USE_ORIGINAL_QUALS && read.getAttribute(RecalDataManager.ORIGINAL_QUAL_ATTRIBUTE_TAG) != null ) {
                            Object obj = read.getAttribute(RecalDataManager.ORIGINAL_QUAL_ATTRIBUTE_TAG);
                            if ( obj instanceof String )
                                quals = QualityUtils.fastqToPhred((String)obj);
                            else {
                                throw new RuntimeException(String.format("Value encoded by %s in %s isn't a string!", RecalDataManager.ORIGINAL_QUAL_ATTRIBUTE_TAG, read.getReadName()));
                            }
                        }

                		// skip if base quality is zero
                		if( quals[offset] > 0 ) {
                		    bases = read.getReadString().toCharArray();
		                    refBase = ref.getBase();
		                    prevBase = bases[offset-1];
		                    // Get the complement base strand if we are a negative strand read
		                    if( read.getReadNegativeStrandFlag() ) {
		                    	bases = BaseUtils.simpleComplement( read.getReadString() ).toCharArray(); // this is an expensive call
		                        refBase = BaseUtils.simpleComplement( refBase );
		                        prevBase = bases[offset+1];
		                    }
		                    
		                    // skip if this base or the previous one was an 'N' or etc.
		                    if( BaseUtils.isRegularBase(prevBase) && BaseUtils.isRegularBase(bases[offset]) ) {
		                    	readGroup = read.getReadGroup().getReadGroupId(); // this is an expensive call
		                    	updateDataFromRead( read, offset, readGroup, quals, bases, refBase );
		                    }
                		}
                	}
                }

            }
            countedSites++;

        } else { // We skipped over the dbSNP site
            skippedSites++;
        }
        
        return 1; // This value isn't actually used anywhere
    }

    /**
     * Major workhorse routine for this walker.
     * Loop through the list of requested covariates and pick out the value from the read, offset, and reference
     * Using the list of covariate values as a key, pick out the RecalDatum and increment,
     *   adding one to the number of observations and potentially one to the number of mismatches
     * Lots of things are passed as parameters to this method as a strategy for optimizing the covariate.getValue calls
     *   because pulling things out of the SAMRecord is an expensive operation.
     * @param read The read
     * @param offset The offset in the read for this locus
     * @param readGroup The read group the read is in
     * @param quals List of base quality scores
     * @param bases The bases which make up the read
     * @param refBase The reference base at this locus
     */
    private void updateDataFromRead(final SAMRecord read, final int offset, final String readGroup, final byte[] quals, final char[] bases, final char refBase) {

        List<Comparable> key = new ArrayList<Comparable>();
        
        // Loop through the list of requested covariates and pick out the value from the read, offset, and reference
        for( Covariate covariate : requestedCovariates ) {
        	key.add( covariate.getValue( read, offset, readGroup, quals, bases, refBase ) );
        }
       
    	// Using the list of covariate values as a key, pick out the RecalDatum from the data HashMap
        RecalDatum datum = dataManager.data.get( key );
        if( datum == null ) { // key doesn't exist yet in the map so make a new bucket and add it
            datum = new RecalDatum(); // initialized with zeros, will be incremented at end of method
            if( validateOldRecalibrator ) {
                dataManager.data.myPut( key, datum );
            } else {
                dataManager.data.put( key, datum );
            }
        }
        
        // Need the bases to determine whether or not we have a mismatch
        char base = bases[offset];
        
        // Add one to the number of observations and potentially one to the number of mismatches
        datum.increment( base, refBase );
        countedBases++;
    }


    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * Initialize the reudce step by creating a PrintStream from the filename specified as an argument to the walker.
     * @return returns A PrintStream created from the -rf filename
     */
    public PrintStream reduceInit() {
        try {
            return new PrintStream( RECAL_FILE );
        } catch ( FileNotFoundException e ) {
            throw new RuntimeException( "Couldn't open output file: ", e );
        }
    }

    /**
     * The Reduce method doesn't do anything for this walker.
     * @param value Result of the map. This value is immediately ignored.
     * @param recalTableStream The PrintStream used to output the CSV data
     * @return returns The PrintStream used to output the CSV data
     */
    public PrintStream reduce( Integer value, PrintStream recalTableStream ) {
        return recalTableStream; // nothing to do here
    }

    /**
     * Write out the full data hashmap to disk in CSV format
     * @param recalTableStream The PrintStream to write out to
     */
    public void onTraversalDone( PrintStream recalTableStream ) {
        out.print( "Writing raw recalibration data..." );
        outputToCSV( recalTableStream );
        out.println( "...done!" );
    
        recalTableStream.close();
    }

    /**
     * For each entry (key-value pair) in the data hashmap output the Covariate's values as well as the RecalDatum's data in CSV format
     * @param recalTableStream The PrintStream to write out to
     */
    private void outputToCSV( final PrintStream recalTableStream ) {

        if( validateOldRecalibrator ) {
            final boolean collapsePos = false;
            final boolean collapseDinuc = false;
            recalTableStream.printf("# collapsed_pos    %b%n", collapsePos);
            recalTableStream.printf("# collapsed_dinuc  %b%n", collapseDinuc);
            recalTableStream.printf("# counted_sites    %d%n", countedSites);
            recalTableStream.printf("# counted_bases    %d%n", countedBases);
            recalTableStream.printf("# skipped_sites    %d%n", skippedSites);
            recalTableStream.printf("# fraction_skipped 1 / %.0f bp%n", (double)countedSites / skippedSites);
            recalTableStream.println("rg,pos,Qrep,dn,nBases,nMismatches,Qemp");
            for(Pair<List<? extends Comparable>, RecalDatum> entry : dataManager.data.entrySetSorted4() ) {
                // For each Covariate in the key
                for( Comparable comp : entry.getFirst() ) {
                    // Output the Covariate's value
                    recalTableStream.print( comp + "," );
                }
                // Output the RecalDatum entry
                recalTableStream.println( entry.getSecond().outputToCSV() );
            }
        } else {
            if( !NO_PRINT_HEADER ) {
                recalTableStream.printf("# Counted Sites    %d%n", countedSites);
                recalTableStream.printf("# Counted Bases    %d%n", countedBases);
                recalTableStream.printf("# Skipped Sites    %d%n", skippedSites);
                recalTableStream.printf("# Fraction Skipped 1 / %.0f bp%n", (double)countedSites / skippedSites);
                for( Covariate cov : requestedCovariates ) {
                    // The "@!" is a code for TableRecalibrationWalker to recognize this line as a Covariate class name
                    recalTableStream.println( "@!" + cov.getClass().getSimpleName() );
                }
            }
            // For each entry in the data hashmap
            for( Map.Entry<List<? extends Comparable>, RecalDatum> entry : dataManager.data.entrySet() ) {
                // For each Covariate in the key
                for( Comparable comp : entry.getKey() ) {
                    // Output the Covariate's value
                    if( NO_PRINT_HEADER && comp instanceof String ) { continue; } // BUGBUG
                    recalTableStream.print( comp + "," );
                }
                // Output the RecalDatum entry
                recalTableStream.println( entry.getValue().outputToCSV() );
            }
        }
    }

}

