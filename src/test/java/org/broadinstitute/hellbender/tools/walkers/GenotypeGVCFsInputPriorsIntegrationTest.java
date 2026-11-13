package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.CommandLineProgramTester;
import org.broadinstitute.hellbender.testutils.GenomicsDBTestUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.tools.walkers.annotator.RMSMappingQuality;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypePriorProvenance;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * GenotypeGVCFs --use-input-genotype-priors on GVCFs built here: sample A (diploid) and sample B (haploid), each
 * with genotype priors (PG) on its variant records. Every route to joint genotyping (a single-sample GVCF,
 * CombineGVCFs and GenomicsDBImport) is run with and without the priors, and the two are compared site by site:
 * <ul>
 *     <li>SITE_SHARED: A and B both have alt allele X. A's likelihoods favor hom-alt and B's favor alt, but with
 *     their priors A is het and B is ref. Both have applicable priors, so QUAL comes from the posteriors.</li>
 *     <li>SITE_DIFFERENT_ALTS: A has alt X, B alt Y. Each sample's <NON_REF> stood for the other's allele, whose
 *     prior is therefore unknown: both keep their likelihood-based genotypes and QUAL.</li>
 *     <li>SITE_A_ONLY: A has alt X, B a reference block. A gets posteriors; B has none, so QUAL stays likelihood-based.</li>
 * </ul>
 */
public class GenotypeGVCFsInputPriorsIntegrationTest extends CommandLineProgramTest {

    private static final String REFERENCE = b37_reference_20_21;
    private static final String CONTIG = "20";
    private static final int REGION_START = 9_999_901;
    private static final int REGION_END = 10_000_300;
    private static final int SITE_SHARED = 10_000_000;
    private static final int SITE_DIFFERENT_ALTS = 10_000_100;
    private static final int SITE_A_ONLY = 10_000_200;
    private static final String SAMPLE_A = "A";
    private static final String SAMPLE_B = "B";
    private static final String SAMPLE_C = "C";

    private static final CommandLineProgramTester GENOTYPE_GVCFS = GenotypeGVCFs.class::getSimpleName;
    private static final CommandLineProgramTester COMBINE_GVCFS = CombineGVCFs.class::getSimpleName;
    private static final CommandLineProgramTester GENOMICSDB_IMPORT = GenomicsDBImport.class::getSimpleName;
    private static final CommandLineProgramTester SELECT_VARIANTS = SelectVariants.class::getSimpleName;

    private SAMSequenceDictionary dictionary;
    private ReferenceDataSource reference;
    private File gvcfA;
    private File gvcfB;
    private File gvcfC;

    @Override
    public String getTestedToolName() {
        return GenotypeGVCFs.class.getSimpleName();
    }

    @BeforeClass
    public void writeGvcfs() {
        reference = ReferenceDataSource.of(IOUtils.getPath(REFERENCE));
        dictionary = reference.getSequenceDictionary();
        gvcfA = writeGvcf(SAMPLE_A, 2, Arrays.asList(SITE_SHARED, SITE_DIFFERENT_ALTS, SITE_A_ONLY), site -> 0);
        gvcfB = writeGvcf(SAMPLE_B, 1, Arrays.asList(SITE_SHARED, SITE_DIFFERENT_ALTS), site -> site == SITE_DIFFERENT_ALTS ? 1 : 0);
        gvcfC = writeGvcf(SAMPLE_C, 2, Collections.emptyList(), site -> 0);
    }

    private Allele referenceAllele(final int position) {
        return Allele.create(new byte[] {reference.queryAndPrefetch(CONTIG, position, position).getBases()[0]}, true);
    }

    /** One of the bases other than the reference: 0 picks the first, 1 the second. */
    private Allele altAllele(final int position, final int which) {
        final String ref = referenceAllele(position).getBaseString();
        final List<String> others = Arrays.stream("ACGT".split("")).filter(b -> !b.equals(ref)).collect(Collectors.toList());
        return Allele.create(others.get(which));
    }

    /**
     * Writes a GVCF with variant records at {@code variantSites} (alt chosen by {@code whichAlt}) and reference blocks elsewhere.
     * Diploid variants: PL 60,2,0 favors hom-alt; PG 0,35,45 makes het the most likely.
     * Haploid variants: PL 30,0 favors alt; PG 0,50 makes ref the most likely. At SITE_DIFFERENT_ALTS the haploid's
     * alt is supported strongly (PL 90,0), so that both samples' alts are kept.
     */
    private File writeGvcf(final String sample, final int ploidy, final List<Integer> variantSites, final Function<Integer, Integer> whichAlt) {
        final List<VariantContext> records = new ArrayList<>();
        int blockStart = REGION_START;
        for (final int site : variantSites) {
            records.add(referenceBlock(sample, ploidy, blockStart, site - 1));
            records.add(variant(sample, ploidy, site, altAllele(site, whichAlt.apply(site))));
            blockStart = site + 1;
        }
        records.add(referenceBlock(sample, ploidy, blockStart, REGION_END));
        return writeGvcf(sample, records);
    }

    private File writeGvcf(final String sample, final List<VariantContext> records) {
        final File gvcf = createTempFile("sample" + sample, ".g.vcf.gz");
        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        lines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));  // as HaplotypeCaller declares it; GenomicsDB drops FORMAT DP without it
        for (final String key : Arrays.asList(VCFConstants.GENOTYPE_KEY, VCFConstants.GENOTYPE_ALLELE_DEPTHS, VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY, VCFConstants.GENOTYPE_PL_KEY)) {
            lines.add(VCFStandardHeaderLines.getFormatLine(key));
        }
        lines.add(GATKVCFHeaderLines.getGenotypePriorsFormatLine());
        lines.add(GATKVCFHeaderLines.getGenotypePosteriorsFormatLine());
        final VCFHeader header = new VCFHeader(lines, Collections.singletonList(sample));
        header.setSequenceDictionary(dictionary);
        try (final VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(gvcf.toPath(), dictionary, false, Options.INDEX_ON_THE_FLY)) {
            writer.writeHeader(header);
            records.forEach(writer::add);
        }
        return gvcf;
    }

    /** Sets PL, PG and, as HaplotypeCaller writes it, GP: PL plus PG, scaled so the best genotype is 0. */
    private static GenotypeBuilder withPriors(final GenotypeBuilder builder, final int[] pl, final double[] pg) {
        final double[] gp = IntStream.range(0, pl.length).mapToDouble(i -> pl[i] + pg[i]).toArray();
        final double best = Arrays.stream(gp).min().getAsDouble();
        return builder.PL(pl).attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, pg)
                .attribute(VCFConstants.GENOTYPE_POSTERIORS_KEY, Arrays.stream(gp).map(x -> x - best).toArray());
    }

    private VariantContext referenceBlock(final String sample, final int ploidy, final int start, final int end) {
        return referenceBlock(sample, ploidy, start, end, false);
    }

    /** A reference block, with reads or, if {@code noCoverage}, without (and so without informative likelihoods). */
    private VariantContext referenceBlock(final String sample, final int ploidy, final int start, final int end, final boolean noCoverage) {
        final Allele ref = referenceAllele(start);
        final Genotype genotype = new GenotypeBuilder(sample, Collections.nCopies(ploidy, ref)).DP(noCoverage ? 0 : 20).GQ(noCoverage ? 0 : 60)
                .PL(noCoverage ? new int[ploidy + 1] : ploidy == 2 ? new int[] {0, 60, 600} : new int[] {0, 60}).make();
        return new VariantContextBuilder("test", CONTIG, start, end, Arrays.asList(ref, Allele.NON_REF_ALLELE))
                .attribute(VCFConstants.END_KEY, end).genotypes(genotype).make();
    }

    private VariantContext variant(final String sample, final int ploidy, final int site, final Allele alt) {
        final Allele ref = referenceAllele(site);
        final GenotypeBuilder genotype = ploidy == 2
                ? withPriors(new GenotypeBuilder(sample, Arrays.asList(alt, alt)).AD(new int[] {8, 12, 0}).DP(20).GQ(2),
                    new int[] {60, 2, 0, 90, 90, 90}, new double[] {0, 35, 45, 40, 80, 50})
                : withPriors(new GenotypeBuilder(sample, Collections.singletonList(alt)).AD(new int[] {6, 14, 0}).DP(20).GQ(30),
                    site == SITE_DIFFERENT_ALTS ? new int[] {90, 0, 120} : new int[] {30, 0, 60}, new double[] {0, 50, 60});
        return new VariantContextBuilder("test", CONTIG, site, site, Arrays.asList(ref, alt, Allele.NON_REF_ALLELE))
                .attribute(VCFConstants.DEPTH_KEY, 20).genotypes(genotype.make()).make();
    }

    private File genotype(final String input, final boolean useInputPriors, final String... extraArgs) {
        final File output = createTempFile("genotyped", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(new File(REFERENCE))
                .add("V", input)
                .addFlag(RMSMappingQuality.RMS_MAPPING_QUALITY_OLD_BEHAVIOR_OVERRIDE_ARGUMENT)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false)
                .addInterval(new SimpleInterval(CONTIG, REGION_START, REGION_END))
                .addOutput(output);
        if (useInputPriors) {
            args.addFlag(GenotypeGVCFs.USE_INPUT_GENOTYPE_PRIORS_LONG_NAME);
        }
        Arrays.stream(extraArgs).forEach(args::addRaw);
        GENOTYPE_GVCFS.runCommandLine(args);
        return output;
    }

    private File combine(final boolean preservePriors, final File... gvcfs) {
        return combine(preservePriors, Arrays.asList(gvcfs));
    }

    private File combine(final boolean preservePriors, final List<File> gvcfs, final String... extraArgs) {
        final File output = createTempFile("combined", ".g.vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder().addReference(new File(REFERENCE)).addOutput(output)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);
        gvcfs.forEach(args::addVCF);
        if (preservePriors) {
            args.addFlag(CombineGVCFs.PRESERVE_GENOTYPE_PRIORS_LONG_NAME);
        }
        Arrays.stream(extraArgs).forEach(args::addRaw);
        COMBINE_GVCFS.runCommandLine(args);
        return output;
    }

    private String importToGenomicsDB(final boolean preservePriors, final File... gvcfs) {
        final File workspace = new File(createTempDir("inputPriorsGenomicsDB"), "workspace");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add(GenomicsDBImport.WORKSPACE_ARG_LONG_NAME, workspace.getAbsolutePath());
        Arrays.stream(gvcfs).forEach(args::addVCF);
        runImport(args, preservePriors);
        return GenomicsDBTestUtils.makeGenomicsDBUri(workspace);
    }

    /** Adds samples to the workspace made by {@link #importToGenomicsDB}, given by its URI. */
    private void importIncrementally(final String workspaceUri, final boolean preservePriors, final File... gvcfs) {
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add(GenomicsDBImport.INCREMENTAL_WORKSPACE_ARG_LONG_NAME, workspaceUri.replaceFirst("^gendb://", ""));
        Arrays.stream(gvcfs).forEach(args::addVCF);
        runImport(args, preservePriors);
    }

    private void runImport(final ArgumentsBuilder args, final boolean preservePriors) {
        args.addInterval(new SimpleInterval(CONTIG, REGION_START, REGION_END))
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);
        if (preservePriors) {
            args.addFlag(GenomicsDBImport.PRESERVE_GENOTYPE_PRIORS_LONG_NAME);
        }
        GENOMICSDB_IMPORT.runCommandLine(args);
    }

    /** The records of a genotyped VCF by start, which must include all three test sites. */
    private static Map<Integer, VariantContext> bySite(final File vcf) {
        final Map<Integer, VariantContext> sites = VariantContextTestUtils.readEntireVCFIntoMemory(vcf.getAbsolutePath()).getRight().stream()
                .collect(Collectors.toMap(VariantContext::getStart, vc -> vc));
        Assert.assertTrue(sites.keySet().containsAll(Arrays.asList(SITE_SHARED, SITE_DIFFERENT_ALTS, SITE_A_ONLY)), vcf + " has sites " + sites.keySet());
        return sites;
    }

    /** The records of a combined GVCF by start. */
    private static Map<Integer, VariantContext> recordsBySite(final File gvcf) {
        return VariantContextTestUtils.readEntireVCFIntoMemory(gvcf.getAbsolutePath()).getRight().stream()
                .collect(Collectors.toMap(VariantContext::getStart, vc -> vc, (a, b) -> a));
    }

    /** The genotype as allele indices, e.g. "0/1", or "1" for a haploid. */
    private static String genotypeString(final VariantContext vc, final String sample) {
        return vc.getGenotype(sample).getAlleles().stream()
                .map(a -> a.isNoCall() ? VCFConstants.EMPTY_ALLELE : Integer.toString(vc.getAlleleIndex(a)))
                .collect(Collectors.joining("/"));
    }

    /**
     * Asserts that every sample at a site has the same genotype (as alleles, since routes may order alternate alleles
     * differently) and GQ in both outputs, and the site the same alleles and QUAL.
     */
    private static void assertSameCalls(final VariantContext actual, final VariantContext expected) {
        Assert.assertEquals(new LinkedHashSet<>(actual.getAlleles()), new LinkedHashSet<>(expected.getAlleles()), "alleles at " + expected.getStart());
        Assert.assertEquals(actual.getPhredScaledQual(), expected.getPhredScaledQual(), 1e-6, "QUAL at " + expected.getStart());
        for (final Genotype genotype : expected.getGenotypes()) {
            final Genotype other = actual.getGenotype(genotype.getSampleName());
            Assert.assertEquals(sortedAlleles(other), sortedAlleles(genotype), genotype.getSampleName() + " at " + expected.getStart());
            Assert.assertEquals(other.getGQ(), genotype.getGQ(), genotype.getSampleName() + " GQ at " + expected.getStart());
        }
    }

    private static List<String> sortedAlleles(final Genotype genotype) {
        return genotype.getAlleles().stream().map(Allele::getDisplayString).sorted().collect(Collectors.toList());
    }

    private static void assertNoTransportFields(final File vcf) {
        final Pair<VCFHeader, List<VariantContext>> output = VariantContextTestUtils.readEntireVCFIntoMemory(vcf.getAbsolutePath());
        Assert.assertNull(output.getLeft().getFormatHeaderLine(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY));
        for (final VariantContext vc : output.getRight()) {
            for (final Genotype genotype : vc.getGenotypes()) {
                Assert.assertFalse(genotype.hasExtendedAttribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY), genotype + " at " + vc.getStart());
            }
        }
        Assert.assertFalse(output.getRight().isEmpty(), vcf + " has no records");
    }

    @Test
    public void testSingleSampleGvcf() {
        final File withPriorsVcf = genotype(gvcfA.getAbsolutePath(), true);
        final Map<Integer, VariantContext> withPriors = bySite(withPriorsVcf);
        final Map<Integer, VariantContext> withoutPriors = bySite(genotype(gvcfA.getAbsolutePath(), false));
        Assert.assertEquals(withPriors.keySet(), new LinkedHashSet<>(Arrays.asList(SITE_SHARED, SITE_DIFFERENT_ALTS, SITE_A_ONLY)));
        for (final int site : withPriors.keySet()) {
            Assert.assertEquals(genotypeString(withoutPriors.get(site), SAMPLE_A), "1/1", "likelihoods favor hom-alt at " + site);
            final VariantContext vc = withPriors.get(site);
            final Genotype a = vc.getGenotype(SAMPLE_A);
            Assert.assertEquals(genotypeString(vc, SAMPLE_A), "0/1", "with the priors het is most likely at " + site);
            Assert.assertEquals(a.getGQ(), 8);
            Assert.assertEquals(a.getExtendedAttribute(VCFConstants.GENOTYPE_POSTERIORS_KEY).toString(), "23,0,8");
            Assert.assertEquals(a.getExtendedAttribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY).toString(), "0,35,45");
            Assert.assertEquals(a.getPL(), new int[] {60, 2, 0});
            // QUAL from posteriors, as HaplotypeCaller computes it: the sample's hom-ref GP
            Assert.assertEquals(vc.getPhredScaledQual(), 23.0, 1e-6);
            Assert.assertNotEquals(withoutPriors.get(site).getPhredScaledQual(), 23.0, 1e-6);
        }
        assertNoTransportFields(withPriorsVcf);
    }

    private void assertJointCalls(final Map<Integer, VariantContext> withPriors, final Map<Integer, VariantContext> withoutPriors) {
        Assert.assertEquals(withPriors.keySet(), withoutPriors.keySet());
        Assert.assertTrue(withPriors.keySet().containsAll(Arrays.asList(SITE_SHARED, SITE_DIFFERENT_ALTS, SITE_A_ONLY)), withPriors.keySet().toString());

        final VariantContext shared = withPriors.get(SITE_SHARED);
        Assert.assertEquals(genotypeString(withoutPriors.get(SITE_SHARED), SAMPLE_A), "1/1");
        Assert.assertEquals(genotypeString(withoutPriors.get(SITE_SHARED), SAMPLE_B), "1");
        Assert.assertEquals(genotypeString(shared, SAMPLE_A), "0/1");
        Assert.assertEquals(genotypeString(shared, SAMPLE_B), "0");
        Assert.assertEquals(shared.getGenotype(SAMPLE_B).getExtendedAttribute(VCFConstants.GENOTYPE_POSTERIORS_KEY).toString(), "0,20");
        // both samples have applicable priors: QUAL from posteriors, the sum of the samples' hom-ref GP (23 + 0)
        Assert.assertEquals(shared.getPhredScaledQual(), 23.0, 1e-6);

        // the other sample's allele has no known prior, so nothing changes
        assertSameCalls(withPriors.get(SITE_DIFFERENT_ALTS), withoutPriors.get(SITE_DIFFERENT_ALTS));
        Assert.assertFalse(withPriors.get(SITE_DIFFERENT_ALTS).getGenotype(SAMPLE_A).hasExtendedAttribute(VCFConstants.GENOTYPE_POSTERIORS_KEY));

        // A gets posteriors; B (a reference block, without priors) keeps its genotype and the QUAL stays likelihood-based
        final VariantContext aOnly = withPriors.get(SITE_A_ONLY);
        Assert.assertEquals(genotypeString(withoutPriors.get(SITE_A_ONLY), SAMPLE_A), "1/1");
        Assert.assertEquals(genotypeString(aOnly, SAMPLE_A), "0/1");
        Assert.assertEquals(genotypeString(aOnly, SAMPLE_B), genotypeString(withoutPriors.get(SITE_A_ONLY), SAMPLE_B));
        Assert.assertEquals(aOnly.getPhredScaledQual(), withoutPriors.get(SITE_A_ONLY).getPhredScaledQual(), 1e-6);
    }

    @Test
    public void testCombineGVCFs() {
        final File combined = combine(true, gvcfA, gvcfB);
        final File withPriors = genotype(combined.getAbsolutePath(), true);
        final File withoutPriors = genotype(combined.getAbsolutePath(), false);
        assertJointCalls(bySite(withPriors), bySite(withoutPriors));
        assertNoTransportFields(withPriors);
        assertNoTransportFields(withoutPriors);  // the provenance only carries priors, whether or not they are used
    }

    @Test
    public void testGenomicsDB() {
        final String workspace = importToGenomicsDB(true, gvcfA, gvcfB);
        final File withPriors = genotype(workspace, true);
        final File withoutPriors = genotype(workspace, false);
        assertJointCalls(bySite(withPriors), bySite(withoutPriors));
        assertNoTransportFields(withPriors);
        assertNoTransportFields(withoutPriors);
    }

    @Test
    public void testRoutesAgree() {
        final Map<Integer, VariantContext> combined = bySite(genotype(combine(true, gvcfA, gvcfB).getAbsolutePath(), true));
        final Map<Integer, VariantContext> genomicsDB = bySite(genotype(importToGenomicsDB(true, gvcfA, gvcfB), true));
        Assert.assertEquals(genomicsDB.keySet(), combined.keySet());
        for (final int site : combined.keySet()) {
            assertSameCalls(genomicsDB.get(site), combined.get(site));
        }
    }

    @Test
    public void testProvenanceSurvivesRepeatedCombining() {
        final File first = combine(true, gvcfA, gvcfB);
        final File second = combine(true, first, gvcfC);
        final Map<Integer, VariantContext> firstSites = recordsBySite(first);
        final Map<Integer, VariantContext> secondSites = recordsBySite(second);
        final Object provenance = firstSites.get(SITE_SHARED).getGenotype(SAMPLE_A).getExtendedAttribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY);
        Assert.assertEquals(provenance, new GenotypePriorProvenance(SITE_SHARED, 2,
                Arrays.asList(referenceAllele(SITE_SHARED), altAllele(SITE_SHARED, 0), Allele.NON_REF_ALLELE), new double[] {0, 35, 45, 40, 80, 50}).encode());
        Assert.assertEquals(secondSites.get(SITE_SHARED).getGenotype(SAMPLE_A).getExtendedAttribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY), provenance);

        // genotyping the twice-combined GVCF gives A and B the same calls as the once-combined one
        final Map<Integer, VariantContext> fromFirst = bySite(genotype(first.getAbsolutePath(), true));
        final Map<Integer, VariantContext> fromSecond = bySite(genotype(second.getAbsolutePath(), true));
        for (final String sample : Arrays.asList(SAMPLE_A, SAMPLE_B)) {
            Assert.assertEquals(genotypeString(fromSecond.get(SITE_SHARED), sample), genotypeString(fromFirst.get(SITE_SHARED), sample));
        }
    }

    @Test
    public void testPriorsWithoutProvenanceAreNotUsed() {
        // combined or imported without --preserve-genotype-priors, PG is copied unchanged through allele merging and cannot be trusted
        final File combined = combine(false, gvcfA, gvcfB);
        final String workspace = importToGenomicsDB(false, gvcfA, gvcfB);
        // nor in a single-sample workspace, since GenomicsDB records are never taken as written by the caller
        final String singleSampleWorkspace = importToGenomicsDB(false, gvcfA);
        // nor when a multi-sample GVCF without provenance is combined again with it
        final File recombined = combine(true, combined, gvcfC);
        for (final String input : Arrays.asList(combined.getAbsolutePath(), workspace, singleSampleWorkspace, recombined.getAbsolutePath())) {
            final Map<Integer, VariantContext> withPriors = bySite(genotype(input, true));
            final Map<Integer, VariantContext> withoutPriors = bySite(genotype(input, false));
            Assert.assertEquals(withPriors.keySet(), withoutPriors.keySet());
            for (final int site : withoutPriors.keySet()) {
                assertSameCalls(withPriors.get(site), withoutPriors.get(site));
            }
        }
    }

    @Test
    public void testDefaultOutputIgnoresPriors() {
        // without the option, the priors in the input change nothing
        final File gvcfWithoutPriors = createTempFile("noPriors", ".g.vcf.gz");
        final List<VariantContext> stripped = new ArrayList<>();
        final Pair<VCFHeader, List<VariantContext>> input = VariantContextTestUtils.readEntireVCFIntoMemory(gvcfA.getAbsolutePath());
        for (final VariantContext vc : input.getRight()) {
            stripped.add(new VariantContextBuilder(vc).genotypes(vc.getGenotypes().stream()
                    .map(g -> new GenotypeBuilder(g).noAttributes().make()).collect(Collectors.toList())).make());
        }
        try (final VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(gvcfWithoutPriors.toPath(), dictionary, false, Options.INDEX_ON_THE_FLY)) {
            writer.writeHeader(input.getLeft());
            stripped.forEach(writer::add);
        }
        final Map<Integer, VariantContext> fromPriors = bySite(genotype(gvcfA.getAbsolutePath(), false));
        final Map<Integer, VariantContext> fromNoPriors = bySite(genotype(gvcfWithoutPriors.getAbsolutePath(), false));
        Assert.assertEquals(fromPriors.keySet(), fromNoPriors.keySet());
        for (final int site : fromNoPriors.keySet()) {
            assertSameCalls(fromPriors.get(site), fromNoPriors.get(site));
        }
    }

    @Test
    public void testNonVariantSitesAndWeakVariants() {
        // a variant too weak to call (PL 0,10,40) is genotyped down to the reference allele; with all sites output, it
        // and the reference blocks are written, without the priors' transport field
        final File gvcf = createTempFile("weak", ".g.vcf.gz");
        final Pair<VCFHeader, List<VariantContext>> input = VariantContextTestUtils.readEntireVCFIntoMemory(gvcfA.getAbsolutePath());
        try (final VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(gvcf.toPath(), dictionary, false, Options.INDEX_ON_THE_FLY)) {
            writer.writeHeader(input.getLeft());
            for (final VariantContext vc : input.getRight()) {
                writer.add(vc.getStart() != SITE_DIFFERENT_ALTS ? vc : new VariantContextBuilder(vc).genotypes(withPriors(new GenotypeBuilder(vc.getGenotype(SAMPLE_A)),
                        new int[] {0, 10, 40, 30, 50, 90}, new double[] {0, 35, 45, 40, 80, 50}).make()).make());
            }
        }
        for (final String allSites : Arrays.asList("--" + GenotypeGVCFs.ALL_SITES_LONG_NAME,
                "--" + GenotypeGVCFs.FORCE_OUTPUT_INTERVALS_NAME + " " + CONTIG + ":" + SITE_DIFFERENT_ALTS)) {
            final File withPriors = genotype(gvcf.getAbsolutePath(), true, allSites);
            assertNoTransportFields(withPriors);
            final Map<Integer, VariantContext> withPriorsSites = recordsBySite(withPriors);
            final Map<Integer, VariantContext> withoutPriorsSites = recordsBySite(genotype(gvcf.getAbsolutePath(), false, allSites));
            Assert.assertEquals(withPriorsSites.keySet(), withoutPriorsSites.keySet(), allSites);
            Assert.assertTrue(withPriorsSites.containsKey(SITE_DIFFERENT_ALTS), allSites);
            Assert.assertEquals(genotypeString(withPriorsSites.get(SITE_SHARED), SAMPLE_A), "0/1", allSites);
        }
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testRejectsSomaticInput() {
        genotype(gvcfA.getAbsolutePath(), true, "--" + CombineGVCFs.SOMATIC_INPUT_LONG_NAME);
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testGenomicsDBImportRejectsHtslibReader() {
        final File workspace = new File(createTempDir("inputPriorsBypass"), "workspace");
        GENOMICSDB_IMPORT.runCommandLine(new ArgumentsBuilder()
                .add(GenomicsDBImport.WORKSPACE_ARG_LONG_NAME, workspace.getAbsolutePath())
                .addInterval(new SimpleInterval(CONTIG, REGION_START, REGION_END))
                .addVCF(gvcfA)
                .addFlag(GenomicsDBImport.BYPASS_FEATURE_READER)
                .addFlag(GenomicsDBImport.PRESERVE_GENOTYPE_PRIORS_LONG_NAME));
    }

    /** gvcfA with sample A's PG at SITE_SHARED replaced. */
    private File gvcfAWithPriorsAtSharedSite(final double[] priors) {
        final File gvcf = createTempFile("modifiedPriors", ".g.vcf.gz");
        final Pair<VCFHeader, List<VariantContext>> input = VariantContextTestUtils.readEntireVCFIntoMemory(gvcfA.getAbsolutePath());
        try (final VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(gvcf.toPath(), dictionary, false, Options.INDEX_ON_THE_FLY)) {
            writer.writeHeader(input.getLeft());
            for (final VariantContext vc : input.getRight()) {
                writer.add(vc.getStart() != SITE_SHARED ? vc : new VariantContextBuilder(vc).genotypes(new GenotypeBuilder(vc.getGenotype(SAMPLE_A))
                        .attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, priors).make()).make());
            }
        }
        return gvcf;
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testRejectsMalformedPriors() {
        genotype(gvcfAWithPriorsAtSharedSite(new double[] {0, -35, 45, 40, 80, 50}).getAbsolutePath(), true);
    }

    @Test
    public void testDoesNotUsePriorsThatDoNotFitTheirRecord() {
        // PG with three values on a record with three alleles, which needs six: the sample keeps its likelihood genotype there
        final Map<Integer, VariantContext> sites = bySite(genotype(gvcfAWithPriorsAtSharedSite(new double[] {0, 35, 45}).getAbsolutePath(), true));
        Assert.assertEquals(genotypeString(sites.get(SITE_SHARED), SAMPLE_A), "1/1");
        Assert.assertEquals(genotypeString(sites.get(SITE_A_ONLY), SAMPLE_A), "0/1");
    }

    @Test
    public void testSampleTakenOutOfCombinedGvcf() {
        // combining without --preserve-genotype-priors adds B's allele to A's record at SITE_DIFFERENT_ALTS and remaps
        // A's PL, but leaves its PG and GP; taken out again, A's priors there no longer fit the record and are not used
        final File combined = combine(false, gvcfA, gvcfB);
        final File extracted = createTempFile("extractedA", ".g.vcf.gz");
        SELECT_VARIANTS.runCommandLine(new ArgumentsBuilder().addReference(new File(REFERENCE)).addVCF(combined).addOutput(extracted)
                .add(StandardArgumentDefinitions.SAMPLE_NAME_LONG_NAME, SAMPLE_A)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false));
        Assert.assertEquals(recordsBySite(extracted).get(SITE_DIFFERENT_ALTS).getNAlleles(), 4, "the extracted record keeps the combined alleles");

        final Map<Integer, VariantContext> fromExtracted = bySite(genotype(extracted.getAbsolutePath(), true));
        final Map<Integer, VariantContext> fromOriginal = bySite(genotype(gvcfA.getAbsolutePath(), true));
        final Map<Integer, VariantContext> withoutPriors = bySite(genotype(extracted.getAbsolutePath(), false));
        assertSameCalls(fromExtracted.get(SITE_DIFFERENT_ALTS), withoutPriors.get(SITE_DIFFERENT_ALTS));
        // where the combined record kept A's own alleles, its priors still fit and are used
        assertSameCalls(fromExtracted.get(SITE_SHARED), fromOriginal.get(SITE_SHARED));
        Assert.assertEquals(genotypeString(fromExtracted.get(SITE_SHARED), SAMPLE_A), "0/1");
    }

    @Test
    public void testSpanningDeletionWithUninformativeSample() {
        // D has a deletion and, inside it, a record with the spanning deletion and another allele, each
        // with priors; E has no reads there. QUAL at the spanning-deletion site is computed from D's posteriors alone.
        final int deletionStart = SITE_A_ONLY + 50;
        final Allele deletionRef = Allele.create(reference.queryAndPrefetch(CONTIG, deletionStart, deletionStart + 1).getBases(), true);
        final Allele deletion = Allele.create(new byte[] {deletionRef.getBases()[0]});
        final Allele inside = referenceAllele(deletionStart + 1);
        final Allele insideAlt = altAllele(deletionStart + 1, 0);
        final VariantContext deletionRecord = new VariantContextBuilder("test", CONTIG, deletionStart, deletionStart + 1,
                Arrays.asList(deletionRef, deletion, Allele.NON_REF_ALLELE)).attribute(VCFConstants.DEPTH_KEY, 20)
                .genotypes(withPriors(new GenotypeBuilder("D", Arrays.asList(deletionRef, deletion)).AD(new int[] {10, 10, 0}).DP(20).GQ(99),
                        new int[] {300, 0, 300, 330, 330, 660}, new double[] {0, 30, 40, 30, 50, 50}).make()).make();
        // genotypes of T, *, C, <NON_REF>: D is */C
        final VariantContext spanningRecord = new VariantContextBuilder("test", CONTIG, deletionStart + 1, deletionStart + 1,
                Arrays.asList(inside, Allele.SPAN_DEL, insideAlt, Allele.NON_REF_ALLELE)).attribute(VCFConstants.DEPTH_KEY, 20)
                .genotypes(withPriors(new GenotypeBuilder("D", Arrays.asList(Allele.SPAN_DEL, insideAlt)).AD(new int[] {0, 10, 10, 0}).DP(20).GQ(99),
                        new int[] {600, 300, 300, 300, 0, 300, 600, 600, 600, 600}, new double[] {0, 30, 40, 30, 40, 50, 50, 60, 60, 60}).make()).make();
        final File gvcfD = writeGvcf("D", Arrays.asList(referenceBlock("D", 2, REGION_START, deletionStart - 1), deletionRecord,
                spanningRecord, referenceBlock("D", 2, deletionStart + 2, REGION_END)));
        final File gvcfE = writeGvcf("E", Arrays.asList(referenceBlock("E", 2, REGION_START, deletionStart - 1),
                referenceBlock("E", 2, deletionStart, REGION_END, true)));

        final File combined = combine(true, gvcfD, gvcfE);
        final Map<Integer, VariantContext> withPriors = recordsBySite(genotype(combined.getAbsolutePath(), true));
        final VariantContext atSpanningSite = withPriors.get(deletionStart + 1);
        Assert.assertNotNull(atSpanningSite, "sites: " + withPriors.keySet());
        Assert.assertTrue(atSpanningSite.getAlleles().contains(Allele.SPAN_DEL), atSpanningSite.toString());
        Assert.assertTrue(atSpanningSite.getGenotype("D").hasExtendedAttribute(VCFConstants.GENOTYPE_POSTERIORS_KEY), "D is genotyped from its priors");
        Assert.assertFalse(atSpanningSite.getGenotype("E").hasExtendedAttribute(VCFConstants.GENOTYPE_POSTERIORS_KEY));
    }

    @Test
    public void testSiteWithReferenceExtendedByAnotherSamplesDeletion() {
        // W has a weak deletion starting at SITE_SHARED, so the combined record's reference allele is extended and A's
        // alleles with it; the deletion is not kept, and A's priors map onto its extended alleles
        final Allele deletionRef = Allele.create(reference.queryAndPrefetch(CONTIG, SITE_SHARED, SITE_SHARED + 1).getBases(), true);
        final Allele deletion = Allele.create(new byte[] {deletionRef.getBases()[0]});
        final VariantContext weakDeletion = new VariantContextBuilder("test", CONTIG, SITE_SHARED, SITE_SHARED + 1,
                Arrays.asList(deletionRef, deletion, Allele.NON_REF_ALLELE)).attribute(VCFConstants.DEPTH_KEY, 20)
                .genotypes(withPriors(new GenotypeBuilder("W", Arrays.asList(deletionRef, deletionRef)).AD(new int[] {18, 2, 0}).DP(20).GQ(20),
                        new int[] {0, 40, 400, 60, 400, 600}, new double[] {0, 30, 40, 30, 50, 50}).make()).make();
        final File gvcfW = writeGvcf("W", Arrays.asList(referenceBlock("W", 2, REGION_START, SITE_SHARED - 1), weakDeletion,
                referenceBlock("W", 2, SITE_SHARED + 2, REGION_END)));

        final File combined = combine(true, gvcfA, gvcfW);
        Assert.assertEquals(recordsBySite(combined).get(SITE_SHARED).getReference(), deletionRef, "the combined reference allele is extended");
        final VariantContext joint = bySite(genotype(combined.getAbsolutePath(), true)).get(SITE_SHARED);
        final VariantContext alone = bySite(genotype(gvcfA.getAbsolutePath(), true)).get(SITE_SHARED);
        Assert.assertEquals(sortedAlleles(joint.getGenotype(SAMPLE_A)), sortedAlleles(alone.getGenotype(SAMPLE_A)));
        Assert.assertEquals(joint.getGenotype(SAMPLE_A).getExtendedAttribute(VCFConstants.GENOTYPE_POSTERIORS_KEY).toString(), "23,0,8");
    }

    @Test
    public void testIncrementalGenomicsDBImport() {
        final String workspace = importToGenomicsDB(true, gvcfA);
        importIncrementally(workspace, true, gvcfB);
        assertJointCalls(bySite(genotype(workspace, true)), bySite(genotype(workspace, false)));
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testIncrementalImportWithPriorsRequiresWorkspaceWithPriors() {
        importIncrementally(importToGenomicsDB(false, gvcfA), true, gvcfB);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testIncrementalImportIntoWorkspaceWithPriorsRequiresThem() {
        importIncrementally(importToGenomicsDB(true, gvcfA), false, gvcfB);
    }

    @Test
    public void testGenomicsDBImportFromSampleNameMap() throws IOException {
        final File sampleMap = createTempFile("sampleMap", ".txt");
        Files.writeString(sampleMap.toPath(), SAMPLE_A + "\t" + gvcfA.getAbsolutePath() + "\n" + SAMPLE_B + "\t" + gvcfB.getAbsolutePath() + "\n");
        final File workspace = new File(createTempDir("inputPriorsSampleMap"), "workspace");
        runImport(new ArgumentsBuilder().add(GenomicsDBImport.WORKSPACE_ARG_LONG_NAME, workspace.getAbsolutePath())
                .add(GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME, sampleMap.getAbsolutePath()), true);
        final String uri = GenomicsDBTestUtils.makeGenomicsDBUri(workspace);
        assertJointCalls(bySite(genotype(uri, true)), bySite(genotype(uri, false)));
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testCombineGVCFsRejectsSomaticInput() {
        combine(true, Collections.singletonList(gvcfA), "--" + CombineGVCFs.SOMATIC_INPUT_LONG_NAME);
    }
}
