package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public class GenotypePriorProvenanceRecorderUnitTest extends GATKBaseTest {

    private static final Allele REF_A = Allele.create("A", true);
    private static final Allele C = Allele.create("C");
    private static final Allele G = Allele.create("G");
    private static final List<Allele> ALLELES = Arrays.asList(REF_A, C, Allele.NON_REF_ALLELE);

    /** As HaplotypeCaller writes it: GP is PL plus PG, scaled so the best genotype is 0. */
    private static Genotype written(final String sample) {
        return new GenotypeBuilder(sample, Arrays.asList(REF_A, C)).PL(new int[] {30, 0, 40, 33, 43, 76})
                .attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, "0,34.77,37.78,39.03,73.8,42.04")
                .attribute(VCFConstants.GENOTYPE_POSTERIORS_KEY, "0,4.77,47.78,42.03,86.8,88.04").make();
    }

    private static VariantContext record(final List<Allele> alleles, final Genotype... genotypes) {
        return new VariantContextBuilder("test", "chr1", 100, 100, alleles).genotypes(genotypes).make();
    }

    private static GenotypePriorProvenanceRecorder recorder(final boolean trustSingleSampleRecords) {
        return new GenotypePriorProvenanceRecorder(trustSingleSampleRecords, "Do the right thing.",
                LogManager.getLogger(GenotypePriorProvenanceRecorderUnitTest.class));
    }

    @Test
    public void testRecordsProvenanceOfSingleSampleRecord() {
        final Genotype genotype = written("s");
        final VariantContext vc = record(ALLELES, genotype);
        final GenotypePriorProvenanceRecorder recorder = recorder(true);
        final Genotype annotated = recorder.apply(vc, true).getGenotype("s");
        Assert.assertEquals(annotated.getExtendedAttribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY),
                "1|100|2|A/C/<NON_REF>|0.0/34.77/37.78/39.03/73.8/42.04");
        Assert.assertEquals(annotated.getExtendedAttribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY), genotype.getExtendedAttribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY));
        Assert.assertFalse(genotype.hasExtendedAttribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY), "the input genotype must not change");

        // a record that already has provenance keeps it, even if its PG has since gone stale
        final VariantContext stale = record(Arrays.asList(REF_A, C, G, Allele.NON_REF_ALLELE), new GenotypeBuilder(annotated).make());
        Assert.assertSame(recorder.apply(stale, true), stale);

        // and a record without priors is returned as it is
        final VariantContext noPriors = record(ALLELES, new GenotypeBuilder("s", Arrays.asList(REF_A, C)).PL(new int[] {30, 0, 40, 33, 43, 76}).make());
        Assert.assertSame(recorder.apply(noPriors, true), noPriors);

        Assert.assertEquals(recorder.getRecordsWithoutProvenance(), 0);
        Assert.assertEquals(recorder.getRecordsWithPriorsNotFittingThem(), 0);
    }

    @Test
    public void testDoesNotTrustMultiSampleOrUntrustedRecords() {
        final GenotypePriorProvenanceRecorder recorder = recorder(true);
        final VariantContext multiSample = record(ALLELES, written("s"), written("t"));
        Assert.assertSame(recorder.apply(multiSample, true), multiSample);
        Assert.assertSame(recorder.apply(multiSample, false), multiSample);
        Assert.assertEquals(recorder.getRecordsWithoutProvenance(), 1, "a record given again at a later position is not counted again");

        // e.g. GenomicsDB, which merges samples' records, even when one sample is returned
        final GenotypePriorProvenanceRecorder untrusting = recorder(false);
        final VariantContext singleSample = record(ALLELES, written("s"));
        Assert.assertSame(untrusting.apply(singleSample, true), singleSample);
        Assert.assertEquals(untrusting.getRecordsWithoutProvenance(), 1);
        untrusting.logSummary();
    }

    @Test
    public void testDoesNotTrustPriorsThatDoNotFitTheRecord() {
        final GenotypePriorProvenanceRecorder recorder = recorder(true);

        // taken out of a combined GVCF whose merge added an allele: PL was remapped, PG and GP were not
        final Genotype added = new GenotypeBuilder(written("s")).PL(new int[] {30, 90, 90, 0, 90, 40, 33, 90, 43, 76}).make();
        final VariantContext withAddedAllele = record(Arrays.asList(REF_A, G, C, Allele.NON_REF_ALLELE), added);
        Assert.assertSame(recorder.apply(withAddedAllele, true), withAddedAllele);

        // ... or reordered the alleles: the length fits, but GP is no longer PL plus PG
        final Genotype reordered = new GenotypeBuilder(written("s")).PL(new int[] {30, 33, 76, 0, 43, 40}).make();
        final VariantContext withReorderedAlleles = record(Arrays.asList(REF_A, Allele.NON_REF_ALLELE, C), reordered);
        Assert.assertSame(recorder.apply(withReorderedAlleles, true), withReorderedAlleles);

        Assert.assertEquals(recorder.getRecordsWithPriorsNotFittingThem(), 2);
        Assert.assertEquals(recorder.getRecordsWithoutProvenance(), 0);
        recorder.logSummary();
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testRejectsMalformedPriors() {
        final Genotype negative = new GenotypeBuilder(written("s")).attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, "0,-3,37.78,39.03,73.8,42.04").make();
        recorder(true).apply(record(ALLELES, negative), true);
    }
}
