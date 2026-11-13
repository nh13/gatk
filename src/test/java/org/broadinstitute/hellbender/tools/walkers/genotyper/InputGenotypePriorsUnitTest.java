package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class InputGenotypePriorsUnitTest extends GATKBaseTest {

    private static final Allele A = Allele.create("A", true);
    private static final Allele C = Allele.create("C");
    private static final Allele G = Allele.create("G");
    private static final List<Allele> RECORD = Arrays.asList(A, C, Allele.NON_REF_ALLELE);
    private static final List<Allele> SITE = Arrays.asList(A, C);

    // record priors for A/A, A/C, C/C, A/<NON_REF>, C/<NON_REF>, <NON_REF>/<NON_REF>
    private static final String PROVENANCE = new GenotypePriorProvenance(100, 2, RECORD, new double[] {0, 35, 45, 40, 80, 50}).encode();

    /** A genotype already assigned from its likelihoods: PL 60,2,0 favors C/C, but with priors 0,35,45 A/C has the highest posterior. */
    private static Genotype likelihoodCalled(final String sample, final String provenance) {
        final GenotypeBuilder builder = new GenotypeBuilder(sample, Arrays.asList(C, C)).PL(new int[] {60, 2, 0}).GQ(2)
                .attribute(GATKVCFConstants.DRAGEN_GENOTYPE_PRIOR_KEY, "0,35,45,40,80,50")
                .attribute(VCFConstants.DEPTH_KEY, 20);
        if (provenance != null) {
            builder.attribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY, provenance);
        }
        return builder.make();
    }

    @Test
    public void testAssignsGenotypeFromPosteriors() {
        final Genotype input = likelihoodCalled("s1", PROVENANCE);
        final InputGenotypePriors priors = new InputGenotypePriors();
        final InputGenotypePriors.Result result = priors.assign("chr1", 100, SITE, GenotypesContext.create(input));

        final Genotype assigned = result.getGenotypes().get("s1");
        // posteriors in Phred, before scaling: 60+0, 2+35, 0+45 = 60, 37, 45; scaled so the best is 0: 23, 0, 8
        Assert.assertEquals(assigned.getAlleles(), Arrays.asList(A, C));
        Assert.assertEquals(VariantContextGetters.getAttributeAsDoubleArray(assigned, VCFConstants.GENOTYPE_POSTERIORS_KEY, () -> null, Double.NaN),
                new double[] {23, 0, 8}, 1e-9);
        Assert.assertEquals(VariantContextGetters.getAttributeAsDoubleArray(assigned, GATKVCFConstants.GENOTYPE_PRIOR_KEY, () -> null, Double.NaN),
                new double[] {0, 35, 45}, 1e-9);
        // GQ: -10 * log10(10^-2.3 + 10^-0.8) = 7.86
        Assert.assertEquals(assigned.getGQ(), 8);
        Assert.assertEquals(assigned.getPL(), new int[] {60, 2, 0});
        Assert.assertEquals(assigned.getExtendedAttribute(VCFConstants.DEPTH_KEY), 20);
        Assert.assertFalse(assigned.hasExtendedAttribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY));
        Assert.assertFalse(assigned.hasExtendedAttribute(GATKVCFConstants.DRAGEN_GENOTYPE_PRIOR_KEY));
        Assert.assertTrue(result.allInformativeSamplesHavePosteriors());

        // the input genotype is unchanged
        Assert.assertEquals(input.getAlleles(), Arrays.asList(C, C));
        Assert.assertTrue(input.hasExtendedAttribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY));

        Assert.assertEquals(priors.getGenotypeCounts(), Collections.singletonMap(InputGenotypePriors.Outcome.POSTERIOR, 1L));
        Assert.assertEquals(priors.getSitesWithPosteriorQual(), 1);
        Assert.assertEquals(priors.getSitesWithLikelihoodQual(), 0);
    }

    @Test
    public void testKeepsLikelihoodGenotypeWithoutApplicablePriors() {
        final InputGenotypePriors priors = new InputGenotypePriors();
        final Genotype noPriors = likelihoodCalled("noPriors", null);
        final Genotype newAllele = likelihoodCalled("newAllele", PROVENANCE);
        final InputGenotypePriors.Result result = priors.assign("chr1", 100, SITE, GenotypesContext.create(noPriors));
        Assert.assertFalse(result.allInformativeSamplesHavePosteriors());
        assertLikelihoodGenotypeKept(result.getGenotypes().get("noPriors"));

        // G was not among the record's alleles: its prior is unknown, although <NON_REF> stood for it
        final Genotype atNewAlleleSite = new GenotypeBuilder(newAllele).PL(new int[] {60, 2, 0, 70, 80, 90}).make();
        final InputGenotypePriors.Result newAlleleResult = priors.assign("chr1", 100, Arrays.asList(A, C, G), GenotypesContext.create(atNewAlleleSite));
        Assert.assertFalse(newAlleleResult.allInformativeSamplesHavePosteriors());
        assertLikelihoodGenotypeKept(newAlleleResult.getGenotypes().get("newAllele"));

        // a spanning deletion's priors were written for the record where it starts
        final InputGenotypePriors.Result spanningResult = priors.assign("chr1", 101, SITE, GenotypesContext.create(newAllele));
        Assert.assertFalse(spanningResult.allInformativeSamplesHavePosteriors());

        // a site that keeps <NON_REF> is left to the likelihoods
        final Genotype withNonRef = new GenotypeBuilder(newAllele).PL(new int[] {60, 2, 0, 70, 80, 90}).make();
        Assert.assertFalse(priors.assign("chr1", 100, RECORD, GenotypesContext.create(withNonRef)).allInformativeSamplesHavePosteriors());

        // priors written for a diploid record do not apply to a haploid genotype
        final Genotype haploid = new GenotypeBuilder("haploid", Collections.singletonList(C)).PL(new int[] {60, 0})
                .attribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY, PROVENANCE).make();
        final InputGenotypePriors.Result haploidResult = priors.assign("chr1", 100, SITE, GenotypesContext.create(haploid));
        Assert.assertFalse(haploidResult.allInformativeSamplesHavePosteriors());
        Assert.assertEquals(haploidResult.getGenotypes().get("haploid").getAlleles(), Collections.singletonList(C));

        // nor to a site whose reference allele differs
        final InputGenotypePriors.Result otherReference = priors.assign("chr1", 100, Arrays.asList(Allele.create("G", true), C),
                GenotypesContext.create(likelihoodCalled("otherReference", PROVENANCE)));
        Assert.assertFalse(otherReference.allInformativeSamplesHavePosteriors());
        assertLikelihoodGenotypeKept(otherReference.getGenotypes().get("otherReference"));

        final Map<InputGenotypePriors.Outcome, Long> counts = priors.getGenotypeCounts();
        Assert.assertEquals(counts.get(InputGenotypePriors.Outcome.NO_PRIORS), Long.valueOf(1));
        Assert.assertEquals(counts.get(InputGenotypePriors.Outcome.NEW_ALLELE), Long.valueOf(1));
        Assert.assertEquals(counts.get(InputGenotypePriors.Outcome.DIFFERENT_START), Long.valueOf(1));
        Assert.assertEquals(counts.get(InputGenotypePriors.Outcome.NON_REF_AT_SITE), Long.valueOf(1));
        Assert.assertEquals(counts.get(InputGenotypePriors.Outcome.DIFFERENT_PLOIDY), Long.valueOf(1));
        Assert.assertEquals(counts.get(InputGenotypePriors.Outcome.DIFFERENT_REFERENCE), Long.valueOf(1));
        Assert.assertFalse(counts.containsKey(InputGenotypePriors.Outcome.POSTERIOR));
        Assert.assertEquals(priors.getSitesWithLikelihoodQual(), 6);
        Assert.assertEquals(priors.getSitesWherePriorsRemovedAlternateCalls(), 0);
        priors.logSummary(LogManager.getLogger(InputGenotypePriorsUnitTest.class));
    }

    @Test
    public void testEveryInapplicabilityHasAnOutcome() {
        for (final GenotypePriorProvenance.Inapplicability inapplicability : GenotypePriorProvenance.Inapplicability.values()) {
            Assert.assertEquals(InputGenotypePriors.Outcome.of(inapplicability).name(), inapplicability.name());
        }
    }

    @Test
    public void testCountsSitesWherePriorsRemoveEveryAlternateCall() {
        // PL 10,0,50 calls A/C; with priors 0,35,45 the posteriors 10,35,95 make it A/A
        final Genotype het = new GenotypeBuilder("het", Arrays.asList(A, C)).PL(new int[] {10, 0, 50})
                .attribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY, PROVENANCE).make();
        final InputGenotypePriors priors = new InputGenotypePriors();
        final InputGenotypePriors.Result result = priors.assign("chr1", 100, SITE, GenotypesContext.create(het));
        Assert.assertEquals(result.getGenotypes().get("het").getAlleles(), Arrays.asList(A, A));
        Assert.assertEquals(priors.getSitesWherePriorsRemovedAlternateCalls(), 1);

        // a site still called variant is not counted
        priors.assign("chr1", 100, SITE, GenotypesContext.create(likelihoodCalled("s1", PROVENANCE)));
        Assert.assertEquals(priors.getSitesWherePriorsRemovedAlternateCalls(), 1);
    }

    private static void assertLikelihoodGenotypeKept(final Genotype genotype) {
        Assert.assertEquals(genotype.getAlleles(), Arrays.asList(C, C));
        Assert.assertEquals(genotype.getGQ(), 2);
        Assert.assertFalse(genotype.hasExtendedAttribute(VCFConstants.GENOTYPE_POSTERIORS_KEY));
        Assert.assertFalse(genotype.hasExtendedAttribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY));
        Assert.assertFalse(genotype.hasExtendedAttribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY));
        Assert.assertFalse(genotype.hasExtendedAttribute(GATKVCFConstants.DRAGEN_GENOTYPE_PRIOR_KEY));
    }

    @Test
    public void testSiteQualRequiresPosteriorsForEveryInformativeSample() {
        final Genotype withPriors = likelihoodCalled("withPriors", PROVENANCE);
        final Genotype withoutPriors = likelihoodCalled("withoutPriors", null);
        final Genotype noData = new GenotypeBuilder("noData", Arrays.asList(Allele.NO_CALL, Allele.NO_CALL)).PL(new int[] {0, 0, 0}).make();
        final Genotype noLikelihoods = new GenotypeBuilder("noLikelihoods", Arrays.asList(A, A)).GQ(20).make();
        final InputGenotypePriors priors = new InputGenotypePriors();

        // a sample with informative likelihoods but no priors keeps the site QUAL likelihood-based ...
        final InputGenotypePriors.Result mixed = priors.assign("chr1", 100, SITE, GenotypesContext.create(withPriors, withoutPriors));
        Assert.assertFalse(mixed.allInformativeSamplesHavePosteriors());
        Assert.assertEquals(mixed.getGenotypes().get("withPriors").getAlleles(), Arrays.asList(A, C), "the sample with priors still gets posteriors");

        // ... while samples without informative likelihoods do not count
        final InputGenotypePriors.Result withUninformative = priors.assign("chr1", 100, SITE, GenotypesContext.create(withPriors, noData, noLikelihoods));
        Assert.assertTrue(withUninformative.allInformativeSamplesHavePosteriors());
        Assert.assertEquals(withUninformative.getGenotypes().get("noData").getAlleles(), noData.getAlleles());
        Assert.assertEquals(withUninformative.getGenotypes().get("noLikelihoods").getAlleles(), noLikelihoods.getAlleles());

        // and with no informative sample at all there are no posteriors to use
        Assert.assertFalse(priors.assign("chr1", 100, SITE, GenotypesContext.create(noData)).allInformativeSamplesHavePosteriors());
    }

    @Test
    public void testHaploidPriors() {
        final String haploidProvenance = new GenotypePriorProvenance(100, 1, RECORD, new double[] {0, 50, 60}).encode();
        // PL 40,0 favors C; prior 50 on C makes A the more likely: posteriors 40, 50 -> 0, 10
        final Genotype haploid = new GenotypeBuilder("h", Collections.singletonList(C)).PL(new int[] {40, 0})
                .attribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY, haploidProvenance).make();
        final Genotype assigned = new InputGenotypePriors().assign("chrX", 100, SITE, GenotypesContext.create(haploid)).getGenotypes().get("h");
        Assert.assertEquals(assigned.getAlleles(), Collections.singletonList(A));
        Assert.assertEquals(VariantContextGetters.getAttributeAsDoubleArray(assigned, VCFConstants.GENOTYPE_POSTERIORS_KEY, () -> null, Double.NaN),
                new double[] {0, 10}, 1e-9);
        Assert.assertEquals(assigned.getGQ(), 10);
    }
}
