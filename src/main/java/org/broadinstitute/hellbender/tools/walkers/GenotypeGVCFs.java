package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantLocusWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBArgumentCollection;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBOptions;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandBiasBySample;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.CombineGVCFs;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotation;
import org.broadinstitute.hellbender.utils.variant.writers.IntervalFilteringVcfWriter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.utils.variant.writers.IntervalFilteringVcfWriter.Mode.STARTS_IN;

/**
 * Perform joint genotyping on one or more samples pre-called with HaplotypeCaller
 *
 * <p>
 * This tool is designed to perform joint genotyping on a single input, which may contain one or many samples. In any
 * case, the input samples must possess genotype likelihoods produced by HaplotypeCaller with `-ERC GVCF` or
 * `-ERC BP_RESOLUTION`.
 *
 *
 * <h3>Input</h3>
 * <p>
 * The GATK4 GenotypeGVCFs tool can take only one input track.  Options are 1) a single single-sample GVCF 2) a single
 * multi-sample GVCF created by CombineGVCFs or 3) a GenomicsDB workspace created by GenomicsDBImport.
 * A sample-level GVCF is produced by HaplotypeCaller with the `-ERC GVCF` setting.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A final VCF in which all samples have been jointly genotyped.
 * </p>
 *
 * <h3>Usage example</h3>
 *
 * <h4>Perform joint genotyping on a singular sample by providing a single-sample GVCF or on a cohort by providing a combined multi-sample GVCF</h4>
 * <pre>
 * gatk --java-options "-Xmx4g" GenotypeGVCFs \
 *   -R Homo_sapiens_assembly38.fasta \
 *   -V input.g.vcf.gz \
 *   -O output.vcf.gz
 * </pre>
 *
 * <h4>Perform joint genotyping on GenomicsDB workspace created with GenomicsDBImport</h4>
 * <pre>
 * gatk --java-options "-Xmx4g" GenotypeGVCFs \
 *   -R Homo_sapiens_assembly38.fasta \
 *   -V gendb://my_database \
 *   -O output.vcf.gz \
 *   --tmp-dir /path/to/large/tmp
 * </pre>
 *
 * <h3>Caveats</h3>
 * <ul>
 *   <li>Only GVCF files produced by HaplotypeCaller (or CombineGVCFs) can be used as input for this tool. Some other
 * programs produce files that they call GVCFs but those lack some important information (accurate genotype likelihoods
 * for every position) that GenotypeGVCFs requires for its operation.</li>
 *   <li>Cannot take multiple GVCF files in one command.</li>
 *   <li>The amount of temporary disk storage required by GenomicsDBImport may exceed what is available in the default location: `/tmp`. The command line argument `--tmp-dir` can be used to specify an alternate temperary storage location with sufficient space.</li>
 * </ul>
 *
 * <h3>Using the genotype priors in the input</h3>
 * <p>HaplotypeCaller assigns genotypes from posteriors when run with {@code --genotype-assignment-method USE_POSTERIOR_PROBABILITIES}
 * (as in {@code --dragen-mode}), using genotype priors such as the DRAGstr STR priors, and writes each sample's priors to the
 * GVCF as FORMAT PG (Illumina DRAGEN writes them as PRI). By default this tool ignores them and assigns genotypes from the
 * likelihoods alone. With {@code --use-input-genotype-priors} it recomputes each sample's posteriors from its likelihoods
 * and its own priors, and assigns GT, GQ, GP and PG from them. GP and PG in the output then come only from this
 * calculation: samples that keep their likelihood-based genotype have none, and PRI is not written.</p>
 * <ul>
 *   <li>A sample's priors are used only where they apply exactly: at a site starting where the sample's record started,
 *   whose alleles are all alleles of that record. At a site with an allele the sample's record did not have (e.g. another
 *   sample's allele, which that sample's {@code <NON_REF>} stands for) or a deletion spanning from an earlier record, the
 *   sample keeps its likelihood-based genotype. Reference blocks carry no priors.</li>
 *   <li>The site QUAL is computed from the posteriors only when every sample with informative likelihoods got posteriors;
 *   otherwise it is computed as usual from the likelihoods. Which alleles are kept, and the LowQual filter, are decided
 *   from the likelihood-based QUAL as usual, so the calls are not guaranteed to match HaplotypeCaller's. If the priors
 *   leave no sample with an alternate allele, the site is written only with {@code --include-non-variant-sites}.</li>
 *   <li>The priors of a single-sample GVCF are taken as written for each record, provided they fit it: one per genotype,
 *   and, where the record has GP, GP equal to PL plus the priors. Priors that do not fit (as when a sample is taken out
 *   of a combined GVCF, whose PL was remapped but whose PG and GP were not) are not used, with a warning. A reordering
 *   that leaves PL unchanged, or a record without GP, cannot be detected this way, so use GVCFs as written by the caller.</li>
 *   <li>For a cohort, combine the GVCFs with {@code CombineGVCFs --preserve-genotype-priors} or import them with
 *   {@code GenomicsDBImport --preserve-genotype-priors}, which record each sample's priors with the alleles they were
 *   written for (FORMAT {@value GATKVCFConstants#GENOTYPE_PRIOR_PROVENANCE_KEY}). Priors in multi-sample input, or in
 *   any GenomicsDB workspace, without that record are not used, with a warning.</li>
 *   <li>Records whose priors were dropped upstream (for example by ReblockGVCF, which removes them when it changes alleles)
 *   keep their likelihood-based genotypes.</li>
 * </ul>
 *
 * <h3>Special note on ploidy</h3>
 * <p>This tool is able to handle any ploidy (or mix of ploidies) intelligently; there is no need to specify ploidy
 * for non-diploid organisms.</p>
 *
 */
@CommandLineProgramProperties(summary = "Perform joint genotyping on a single-sample GVCF from HaplotypeCaller or a multi-sample GVCF from CombineGVCFs or GenomicsDBImport",
        oneLineSummary = "Perform joint genotyping on one or more samples pre-called with HaplotypeCaller",
        programGroup = ShortVariantDiscoveryProgramGroup.class)
@DocumentedFeature
public final class GenotypeGVCFs extends VariantLocusWalker {

    public static final String PHASED_HOM_VAR_STRING = "1|1";
    public static final String ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME = "only-output-calls-starting-in-intervals";
    public static final String ALL_SITES_LONG_NAME = "include-non-variant-sites";
    public static final String USE_INPUT_GENOTYPE_PRIORS_LONG_NAME = "use-input-genotype-priors";
    public static final String ALL_SITES_SHORT_NAME = "all-sites";
    public static final String KEEP_COMBINED_LONG_NAME = "keep-combined-raw-annotations";
    public static final String KEEP_COMBINED_SHORT_NAME = "keep-combined";
    public static final String FORCE_OUTPUT_INTERVALS_NAME = "force-output-intervals";
    public static final String SOMATIC_QUALITY_THRESHOLD_NAME = "somatic-quality-threshold";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written", optional=false)
    private GATKPath outputFile;

    @Argument(fullName=ALL_SITES_LONG_NAME, shortName=ALL_SITES_SHORT_NAME, doc="Include loci found to be non-variant after genotyping", optional=true)
    private boolean includeNonVariants = false;

    /**
     * Import all data between specified intervals.   Improves performance using large lists of intervals, as in exome
     * sequencing, especially if GVCF data only exists for specified intervals.  Use with
     * --{@value StandardArgumentDefinitions#VARIANT_OUTPUT_INTERVAL_FILTERING_MODE_LONG_NAME} if input GVCFs contain calls outside the specified intervals.
     */
    @Argument(fullName = GenomicsDBImport.MERGE_INPUT_INTERVALS_LONG_NAME,
            shortName = GenomicsDBImport.MERGE_INPUT_INTERVALS_LONG_NAME,
            doc = "Boolean flag to import all data in between intervals.")
    private boolean mergeInputIntervals = false;

    /**
     * "Genotype" somatic GVCFs, outputting genotypes according to confidently called alt alleles, which may lead to inconsistent ploidy
     * Note that the Mutect2 reference confidence mode is in BETA -- the likelihoods model and output format are subject to change in subsequent versions.
     */
    @Argument(fullName= CombineGVCFs.SOMATIC_INPUT_LONG_NAME, doc = "Finalize input GVCF according to somatic (i.e. Mutect2) TLODs (BETA feature)")
    protected boolean somaticInput = false;

    /**
     * Only variants with tumor LODs exceeding this threshold will be written to the VCF, regardless of filter status.
     * Set to less than or equal to tumor_lod. Increase argument value to reduce false positives in the callset.
     */
    @Argument(fullName=M2ArgumentCollection.EMISSION_LOD_LONG_NAME, shortName = M2ArgumentCollection.EMISSION_LOG_SHORT_NAME,
    doc = "LOD threshold to emit variant to VCF.")
    protected double tlodThreshold = 3.5;  //allow for some lower quality variants

    /**
     * Only variants with somatic quality (SQ) exceeding this threshold will be written to the VCF, regardless of filter status.
     * Use with --{@value CombineGVCFs#SOMATIC_INPUT_LONG_NAME}
     * Increase argument value to reduce false positives in the callset.
     */
    @Argument(fullName=SOMATIC_QUALITY_THRESHOLD_NAME, doc = "SQ threshold to emit variant to VCF.")
    protected double sqThreshold = 3.0;  //default to threshold used by DRAGEN https://support-docs.illumina.com/SW/DRAGEN_v40/Content/SW/DRAGEN/MitochondrialCalling.htm

    /**
     * Margin of error in allele fraction to consider a somatic variant homoplasmic, i.e. if there is less than a 0.1% reference allele fraction, those reads are likely errors
     */
    @Argument(fullName=CombineGVCFs.ALLELE_FRACTION_DELTA_LONG_NAME, doc = "Margin of error in allele fraction to consider a somatic variant homoplasmic")
    protected double afTolerance = 1e-3;  //based on Q30 as a "good" base quality score

    /**
     * If specified, keep all the combined raw annotations (e.g. AS_SB_TABLE) after genotyping.  This is applicable to Allele-Specific annotations. See {@link ReducibleAnnotation}
     */
    @Argument(fullName=KEEP_COMBINED_LONG_NAME, shortName = KEEP_COMBINED_SHORT_NAME, doc = "If specified, keep the combined raw annotations",
            mutex = {GenotypeGVCFsAnnotationArgumentCollection.KEEP_SPECIFIED_RAW_COMBINED_ANNOTATION_LONG_NAME})
    protected boolean keepCombined = false;

    @ArgumentCollection
    private GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    /**
     * Assign genotypes from posteriors computed with each sample's own genotype priors from the input (PG, or PRI from
     * Illumina DRAGEN), where they apply, and compute the site QUAL from those posteriors when every sample with
     * informative likelihoods has them. See "Using the genotype priors in the input" above.
     */
    @Argument(fullName = USE_INPUT_GENOTYPE_PRIORS_LONG_NAME, optional = true,
            doc = "Assign genotypes from posteriors using each sample's genotype priors (PG or PRI) from the input where they apply")
    private boolean useInputGenotypePriors = false;

    @ArgumentCollection
    private GenomicsDBArgumentCollection genomicsdbArgs = new GenomicsDBArgumentCollection();

    /**
     * This option can only be activated if intervals are specified.
     */
    @DeprecatedFeature
    @Advanced
    @Argument(fullName= ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME,
            doc="Restrict variant output to sites that start within provided intervals, equivalent to '--"+StandardArgumentDefinitions.VARIANT_OUTPUT_INTERVAL_FILTERING_MODE_LONG_NAME+" STARTS_IN'",
            optional=true,
            mutex = {StandardArgumentDefinitions.VARIANT_OUTPUT_INTERVAL_FILTERING_MODE_LONG_NAME})
    private boolean onlyOutputCallsStartingInIntervals = false;

    @Argument(fullName = FORCE_OUTPUT_INTERVALS_NAME,
            suppressFileExpansion = true, doc = "sites at which to output genotypes even if non-variant in samples", optional = true)
    protected final List<String> forceOutputIntervalStrings = new ArrayList<>();

    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set
     * when appropriate. Note that dbSNP is not used in any way for the genotyping calculations themselves.
     */
    @ArgumentCollection
    private final DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    // @ArgumentCollection deliberately omitted since this is passed to the annotation plugin
    final GenotypeGVCFsAnnotationArgumentCollection genotypeGVCFsAnnotationArgs = new GenotypeGVCFsAnnotationArgumentCollection();

    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    private ReferenceConfidenceVariantContextMerger merger;

    private VariantContextWriter vcfWriter;

    private OverlapDetector<GenomeLoc> forceOutputIntervals;

    private boolean forceOutputIntervalsPresent;

    private GenotypeGVCFsEngine gvcfEngine;

    /**
     * Get the largest interval per contig that contains the intervals specified on the command line.
     * @param getIntervals intervals to be transformed
     * @param sequenceDictionary used to validate intervals
     * @return a list of one interval per contig spanning the input intervals after processing and validation
     */
    @Override
    protected List<SimpleInterval> transformTraversalIntervals(final List<SimpleInterval> getIntervals, final SAMSequenceDictionary sequenceDictionary) {
        if (mergeInputIntervals) {
            return IntervalUtils.getSpanningIntervals(getIntervals, sequenceDictionary);
        } else {
            return getIntervals;
        }
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected GenomicsDBOptions getGenomicsDBOptions() {
        if (genomicsDBOptions == null) {
            //extract called genotypes so hom refs with no PLs aren't ambiguous
            genomicsdbArgs.callGenotypes = true;
            genomicsDBOptions = new GenomicsDBOptions(referenceArguments.getReferencePath(), genomicsdbArgs, genotypeArgs);
        }
        return genomicsDBOptions;
    }

    @Override
    public boolean useVariantAnnotations() { return true;}

    @Override
    public List<? extends CommandLinePluginDescriptor<?>> getPluginDescriptors() {
        GATKReadFilterPluginDescriptor readFilterDescriptor = new GATKReadFilterPluginDescriptor(getDefaultReadFilters());
        return useVariantAnnotations()?
                Arrays.asList(readFilterDescriptor, new GATKAnnotationPluginDescriptor(
                        genotypeGVCFsAnnotationArgs,
                        getDefaultVariantAnnotations(), getDefaultVariantAnnotationGroups())):
                Collections.singletonList(readFilterDescriptor);
    }

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Arrays.asList(StandardAnnotation.class);
    }

    @Override
    public void onTraversalStart() {

        if (somaticInput) {
            logger.warn("Note that the Mutect2 reference confidence mode is in BETA -- the likelihoods model and output format are subject to change in subsequent versions.");
        }

        if (onlyOutputCallsStartingInIntervals) {
            logger.warn("The --" + ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME + " option is deprecated. Please use '--" + StandardArgumentDefinitions.VARIANT_OUTPUT_INTERVAL_FILTERING_MODE_LONG_NAME + " STARTS_IN' for an equivalent filtering.");
            this.userOutputVariantIntervalFilteringMode = STARTS_IN;
        }

        forceOutputIntervalsPresent = !forceOutputIntervalStrings.isEmpty();

        if (useInputGenotypePriors) {
            validateInputGenotypePriorArguments();
        }

        if (includeNonVariants && forceOutputIntervalsPresent ) {
            throw new CommandLineException.BadArgumentValue(String.format("Force output (--%s) is incompatible with including " +
                    "non-variants (--%s and --%s).  Use the latter to force genotyping at all sites and the former to force genotyping only at given sites." +
                    "In both cases, variant sites are genotyped as usual.", FORCE_OUTPUT_INTERVALS_NAME, ALL_SITES_LONG_NAME, ALL_SITES_SHORT_NAME));
        }

        final GenomeLocSortedSet forceOutputLocs = IntervalUtils.loadIntervals(forceOutputIntervalStrings, IntervalSetRule.UNION,
                 IntervalMergingRule.ALL, 0, new GenomeLocParser(getSequenceDictionaryForDrivingVariants()));

        forceOutputIntervals = OverlapDetector.create(forceOutputLocs.toList());

        if (!(includeNonVariants || forceOutputIntervalsPresent)) {
            changeTraversalModeToByVariant();
        }

        final VCFHeader inputVCFHeader = getHeaderForVariants();

        final Collection<Annotation> variantAnnotations = makeVariantAnnotations();
        final Set<Annotation> annotationsToKeep = getAnnotationsToKeep();
        annotationEngine = new VariantAnnotatorEngine(variantAnnotations, dbsnp.dbsnp, Collections.emptyList(), false, keepCombined, annotationsToKeep);

        merger = new ReferenceConfidenceVariantContextMerger(annotationEngine, getHeaderForVariants(), somaticInput, false, true);

        //methods that cannot be called in engine bc its protected
        final Set<VCFHeaderLine> defaultToolVCFHeaderLines = getDefaultToolVCFHeaderLines();
        vcfWriter = createVCFWriter(outputFile);

        //create engine object
        final boolean   keepSB = variantAnnotations.stream().map(a -> a.getClass().getSimpleName()).collect(Collectors.toList()).contains(StrandBiasBySample.class.getSimpleName());
        gvcfEngine = new GenotypeGVCFsEngine(annotationEngine, genotypeArgs, includeNonVariants, inputVCFHeader, keepSB,
                useInputGenotypePriors, !IOUtils.isGenomicsDBPath(drivingVariantFile));

        //call initialize method in engine class that creates VCFWriter object and writes a header to it
        vcfWriter = gvcfEngine.setupVCFWriter(defaultToolVCFHeaderLines, keepCombined, dbsnp, vcfWriter);
    }

    private void validateInputGenotypePriorArguments() {
        if (somaticInput) {
            throw new CommandLineException.BadArgumentValue(String.format("--%s cannot be used with --%s: somatic GVCFs have no genotype priors",
                    USE_INPUT_GENOTYPE_PRIORS_LONG_NAME, CombineGVCFs.SOMATIC_INPUT_LONG_NAME));
        }
    }

    private Set<Annotation> getAnnotationsToKeep() {
        final GATKAnnotationPluginDescriptor pluginDescriptor = getCommandLineParser().getPluginDescriptor(GATKAnnotationPluginDescriptor.class);
        final List<String> annotationStringsToKeep = genotypeGVCFsAnnotationArgs.getKeepSpecifiedCombinedAnnotationNames();
        final Map<String, Annotation> resolvedInstancesMap = pluginDescriptor.getResolvedInstancesMap();
        return annotationStringsToKeep.stream()
                .peek(s -> {Annotation a = resolvedInstancesMap.get(s); if (a == null)
                    throw new UserException("Requested --" + GenotypeGVCFsAnnotationArgumentCollection.KEEP_SPECIFIED_RAW_COMBINED_ANNOTATION_LONG_NAME + ": " + s + " was not found in annotation list. Was it excluded with --" + StandardArgumentDefinitions.ANNOTATIONS_TO_EXCLUDE_LONG_NAME + " or not provided with --" + StandardArgumentDefinitions.ANNOTATION_LONG_NAME + "?"); })
                .map(resolvedInstancesMap::get)
                .collect(Collectors.toSet());
    }

    @Override
    public void apply(final Locatable loc, List<VariantContext> variants, ReadsContext reads, ReferenceContext ref, FeatureContext features) {

        final boolean inForceOutputIntervals = forceOutputIntervalsPresent && forceOutputIntervals.overlapsAny(loc);
        final boolean forceOutput = includeNonVariants || inForceOutputIntervals;
        final VariantContext regenotypedVC = gvcfEngine.callRegion(loc, variants, ref, features, merger, somaticInput, tlodThreshold, sqThreshold, afTolerance, forceOutput);

        if (regenotypedVC != null) {
            if ((forceOutput || !GATKVariantContextUtils.isSpanningDeletionOnly(regenotypedVC))) {
                vcfWriter.add(regenotypedVC);
            }
        }
    }

    @Override
    public void closeTool() {
        if (gvcfEngine != null) {
            gvcfEngine.logInputGenotypePriorSummary();
        }
        if ( vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
