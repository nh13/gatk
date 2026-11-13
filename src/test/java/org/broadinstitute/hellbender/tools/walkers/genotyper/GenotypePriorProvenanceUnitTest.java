package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.IntStream;

public class GenotypePriorProvenanceUnitTest extends GATKBaseTest {

    private static final Allele REF_A = Allele.create("A", true);
    private static final Allele C = Allele.create("C");
    private static final Allele G = Allele.create("G");
    private static final Allele T = Allele.create("T");
    private static final Allele NON_REF = Allele.NON_REF_ALLELE;

    /** Distinct priors 0, 1, 2, ... so that every genotype's prior identifies the record genotype it came from. */
    private static double[] indexPriors(final int ploidy, final int alleleCount) {
        return IntStream.range(0, genotypeOrder(ploidy, alleleCount).size()).mapToDouble(i -> i).toArray();
    }

    /**
     * The VCF order of genotypes, written independently of GenotypeIndexCalculator: each genotype is its allele indices
     * in non-decreasing order, and genotypes are ordered by their last index, then the one before it, and so on.
     */
    private static List<int[]> genotypeOrder(final int ploidy, final int alleleCount) {
        final List<int[]> genotypes = new ArrayList<>();
        enumerate(new int[ploidy], 0, 0, alleleCount, genotypes);
        genotypes.sort((a, b) -> {
            for (int i = ploidy - 1; i >= 0; i--) {
                if (a[i] != b[i]) return Integer.compare(a[i], b[i]);
            }
            return 0;
        });
        return genotypes;
    }

    private static void enumerate(final int[] current, final int position, final int minimum, final int alleleCount, final List<int[]> out) {
        if (position == current.length) {
            out.add(current.clone());
            return;
        }
        for (int a = minimum; a < alleleCount; a++) {
            current[position] = a;
            enumerate(current, position + 1, a, alleleCount, out);
        }
    }

    /** The expected log10 priors at a site, by mapping each site genotype's alleles to record alleles through the oracle order. */
    private static double[] expectedLog10Priors(final int ploidy, final int recordAlleleCount, final double[] phredPriors, final int[] siteToRecordAllele) {
        final List<int[]> recordOrder = genotypeOrder(ploidy, recordAlleleCount);
        final List<int[]> siteOrder = genotypeOrder(ploidy, siteToRecordAllele.length);
        final double[] result = new double[siteOrder.size()];
        for (int i = 0; i < siteOrder.size(); i++) {
            final int[] mapped = Arrays.stream(siteOrder.get(i)).map(a -> siteToRecordAllele[a]).sorted().toArray();
            final int recordIndex = IntStream.range(0, recordOrder.size()).filter(j -> Arrays.equals(recordOrder.get(j), mapped)).findFirst().getAsInt();
            result[i] = -0.1 * phredPriors[recordIndex];
        }
        return result;
    }

    @DataProvider
    public Object[][] applicableMappings() {
        return new Object[][] {
                // ploidy, record alleles, site alleles, site-to-record allele map
                {2, Arrays.asList(REF_A, C, G, NON_REF), Arrays.asList(REF_A, C, G, NON_REF), new int[] {0, 1, 2, 3}},  // identity
                {2, Arrays.asList(REF_A, C, G, NON_REF), Arrays.asList(REF_A, G, C), new int[] {0, 2, 1}},               // permutation, NON_REF removed
                {2, Arrays.asList(REF_A, C, G, NON_REF), Arrays.asList(REF_A, G), new int[] {0, 2}},                     // restriction
                {1, Arrays.asList(REF_A, C, G, NON_REF), Arrays.asList(REF_A, G, C), new int[] {0, 2, 1}},               // haploid
                {3, Arrays.asList(REF_A, C, G, NON_REF), Arrays.asList(REF_A, G, C, NON_REF), new int[] {0, 2, 1, 3}},   // triploid
                {4, Arrays.asList(REF_A, C, G, T), Arrays.asList(REF_A, T, C), new int[] {0, 3, 1}},                     // tetraploid
                {2, Arrays.asList(REF_A, Allele.SPAN_DEL, C, NON_REF), Arrays.asList(REF_A, C, Allele.SPAN_DEL), new int[] {0, 2, 1}}, // spanning deletion recorded
        };
    }

    @Test(dataProvider = "applicableMappings")
    public void testMapToApplicableSite(final int ploidy, final List<Allele> recordAlleles, final List<Allele> siteAlleles, final int[] siteToRecordAllele) {
        final double[] phred = indexPriors(ploidy, recordAlleles.size());
        final GenotypePriorProvenance provenance = new GenotypePriorProvenance(100, ploidy, recordAlleles, phred);
        final GenotypePriorProvenance.Mapping mapping = provenance.mapTo(100, siteAlleles, ploidy);
        Assert.assertTrue(mapping.isApplicable());
        Assert.assertEquals(mapping.getLog10Priors(), expectedLog10Priors(ploidy, recordAlleles.size(), phred, siteToRecordAllele), 1e-12);
    }

    @Test
    public void testMapToSiteWithExtendedReference() {
        // the record's C is CTG at a site whose reference allele was extended to ATG by another sample's deletion
        final List<Allele> recordAlleles = Arrays.asList(REF_A, C, NON_REF);
        final double[] phred = indexPriors(2, 3);
        final GenotypePriorProvenance provenance = new GenotypePriorProvenance(100, 2, recordAlleles, phred);
        final GenotypePriorProvenance.Mapping mapping = provenance.mapTo(100, Arrays.asList(Allele.create("ATG", true), Allele.create("CTG")), 2);
        Assert.assertTrue(mapping.isApplicable());
        Assert.assertEquals(mapping.getLog10Priors(), expectedLog10Priors(2, 3, phred, new int[] {0, 1}), 1e-12);

        // the other sample's deletion (A) is not one of this record's alleles
        final GenotypePriorProvenance.Mapping withDeletion = provenance.mapTo(100, Arrays.asList(Allele.create("ATG", true), Allele.create("CTG"), Allele.create("A")), 2);
        Assert.assertEquals(withDeletion.getInapplicability(), GenotypePriorProvenance.Inapplicability.NEW_ALLELE);
    }

    @DataProvider
    public Object[][] inapplicableMappings() {
        final List<Allele> record = Arrays.asList(REF_A, C, NON_REF);
        return new Object[][] {
                {record, 101, Arrays.asList(REF_A, C), 2, GenotypePriorProvenance.Inapplicability.DIFFERENT_START},
                {record, 100, Arrays.asList(REF_A, C), 1, GenotypePriorProvenance.Inapplicability.DIFFERENT_PLOIDY},
                {record, 100, Arrays.asList(Allele.create("G", true), C), 2, GenotypePriorProvenance.Inapplicability.DIFFERENT_REFERENCE},
                {Arrays.asList(Allele.create("AT", true), C, NON_REF), 100, Arrays.asList(REF_A, C), 2, GenotypePriorProvenance.Inapplicability.DIFFERENT_REFERENCE},
                {record, 100, Arrays.asList(REF_A, C, G), 2, GenotypePriorProvenance.Inapplicability.NEW_ALLELE},   // G stood behind <NON_REF>
                {record, 100, Arrays.asList(REF_A, G), 2, GenotypePriorProvenance.Inapplicability.NEW_ALLELE},
                {record, 100, Arrays.asList(REF_A, C, Allele.SPAN_DEL), 2, GenotypePriorProvenance.Inapplicability.NEW_ALLELE},
        };
    }

    @Test(dataProvider = "inapplicableMappings")
    public void testMapToInapplicableSite(final List<Allele> recordAlleles, final int siteStart, final List<Allele> siteAlleles, final int sitePloidy,
                                         final GenotypePriorProvenance.Inapplicability expected) {
        final GenotypePriorProvenance provenance = new GenotypePriorProvenance(100, 2, recordAlleles, indexPriors(2, recordAlleles.size()));
        final GenotypePriorProvenance.Mapping mapping = provenance.mapTo(siteStart, siteAlleles, sitePloidy);
        Assert.assertFalse(mapping.isApplicable());
        Assert.assertEquals(mapping.getInapplicability(), expected);
    }

    @Test
    public void testEncodeDecodeRoundTrip() {
        final List<Allele> alleles = Arrays.asList(Allele.create("ACGT", true), Allele.create("A"), Allele.SPAN_DEL,
                Allele.create("<DEL:ME:ALU>"), Allele.create("A[chr2:123["), NON_REF);
        final double[] phred = new double[GenotypeIndexCalculator.genotypeCount(2, alleles.size())];
        for (int i = 0; i < phred.length; i++) {
            phred[i] = i * 3.77 + 1e-5;
        }
        final GenotypePriorProvenance provenance = new GenotypePriorProvenance(123456, 2, alleles, phred);
        final String encoded = provenance.encode();
        Assert.assertTrue(encoded.startsWith(GenotypePriorProvenance.VERSION + "|123456|2|ACGT/A/*/"), encoded);
        for (final String forbidden : Arrays.asList(":", ",", ";", "=", "\t", " ", "[")) {
            Assert.assertFalse(encoded.contains(forbidden), "encoded value contains '" + forbidden + "': " + encoded);
        }
        final GenotypePriorProvenance decoded = GenotypePriorProvenance.decode(encoded);
        Assert.assertEquals(decoded, provenance);
        Assert.assertEquals(decoded.getAlleles(), alleles);
        Assert.assertEquals(decoded.getPhredPriors(), phred);
    }

    @DataProvider
    public Object[][] malformedValues() {
        return new Object[][] {
                {"2|100|2|A/C|0/1/2"},           // unsupported version
                {"1|100|2|A/C|0/1"},             // too few priors for 2 alleles at ploidy 2
                {"1|100|2|A/C|0/1/2/3"},         // too many priors
                {"1|100|2|A/C|0/x/2"},           // not a number
                {"1|100|2|A/C|0/-1/2"},          // negative
                {"1|100|2|A/C|0/NaN/2"},         // not finite
                {"1|100|2|A/C|0/Infinity/2"},
                {"1|0|2|A/C|0/1/2"},             // start must be positive
                {"1|100|0|A/C|0"},               // ploidy must be positive
                {"1|100|2|A/A|0/1/2"},           // duplicate alleles
                {"1|100|2|A//C|0/1/2/3/4/5"},    // empty allele
                {"1|100|2|A/%3|0/1/2"},          // truncated percent-encoding
                {"1|100|2|A/C:G|0/1/2"},         // unencoded reserved character
                {"1|100|2|A/C"},                 // missing field
                {"1|100|2|A/C|0/1/2|extra"},     // extra field
        };
    }

    @Test(dataProvider = "malformedValues", expectedExceptions = IllegalArgumentException.class)
    public void testDecodeRejectsMalformedValues(final String value) {
        GenotypePriorProvenance.decode(value);
    }

    @Test
    public void testFromGenotypeRejectsMalformedValueWithLocation() {
        final Genotype genotype = new GenotypeBuilder("sample1", Arrays.asList(REF_A, C))
                .attribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY, "1|100|2|A/C|0/1").make();
        try {
            GenotypePriorProvenance.fromGenotype(genotype, "chr1", 100);
            Assert.fail("expected a malformed provenance to be rejected");
        } catch (final UserException.BadInput e) {
            Assert.assertTrue(e.getMessage().contains("sample1") && e.getMessage().contains("chr1:100"), e.getMessage());
        }
        Assert.assertNull(GenotypePriorProvenance.fromGenotype(new GenotypeBuilder("sample1", Arrays.asList(REF_A, C)).make(), "chr1", 100));
        Assert.assertNull(GenotypePriorProvenance.fromGenotype(new GenotypeBuilder("sample1", Arrays.asList(REF_A, C))
                .attribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY, ".").make(), "chr1", 100));
    }

    private static VariantContext record(final List<Allele> alleles, final Genotype... genotypes) {
        return new VariantContextBuilder("test", "chr1", 100, 100, alleles).genotypes(genotypes).make();
    }

    private static GenotypeBuilder diploid(final String sample) {
        return new GenotypeBuilder(sample, Arrays.asList(REF_A, C)).PL(new int[] {30, 0, 40, 33, 43, 76});
    }

    @DataProvider
    public Object[][] recordPriors() {
        final double[] priors = {0, 34.77, 37.78, 39.03, 73.8, 42.04};
        return new Object[][] {
                // genotype, expected priors (null for none)
                {diploid("s").attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, "0,34.77,37.78,39.03,73.8,42.04").make(), priors},  // as read from VCF text
                {diploid("s").attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, priors).make(), priors},
                {diploid("s").attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, Arrays.asList("0", "34.77", "37.78", "39.03", "73.8", "42.04")).make(), priors},
                {diploid("s").attribute(GATKVCFConstants.DRAGEN_GENOTYPE_PRIOR_KEY, "0,34.77,37.78,39.03,73.8,42.04").make(), priors},
                {diploid("s").attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, priors)
                        .attribute(GATKVCFConstants.DRAGEN_GENOTYPE_PRIOR_KEY, "0,34.77,37.78,39.03,73.8,42.04").make(), priors},
                {diploid("s").make(), null},
                {diploid("s").attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, ".,.,.,.,.,.").make(), null},
                {diploid("s").attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, ".").make(), null},
                {new GenotypeBuilder("s", Collections.emptyList()).attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, priors).make(), null},  // no ploidy
        };
    }

    @Test(dataProvider = "recordPriors")
    public void testReadCallerPriors(final Genotype genotype, final double[] expectedPriors) {
        final VariantContext vc = record(Arrays.asList(REF_A, C, NON_REF), genotype);
        final double[] priors = GenotypePriorProvenance.readCallerPriors(vc, genotype);
        if (expectedPriors == null) {
            Assert.assertNull(priors);
        } else {
            Assert.assertEquals(priors, expectedPriors);
        }
    }

    @DataProvider
    public Object[][] malformedRecordPriors() {
        return new Object[][] {
                {diploid("s").attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, "0,.,37.78,39.03,73.8,42.04").make()}, // partly missing
                {diploid("s").attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, "0,x,37.78,39.03,73.8,42.04").make()}, // not a number
                {diploid("s").attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, "0,-3,37.78,39.03,73.8,42.04").make()}, // negative
                {diploid("s").attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, "0,Infinity,37.78,39.03,73.8,42.04").make()}, // not finite
                {diploid("s").attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, "0,34.77,37.78,39.03,73.8,42.04")
                        .attribute(GATKVCFConstants.DRAGEN_GENOTYPE_PRIOR_KEY, "0,34.77,37.78,39.03,73.8,40").make()},  // PG and PRI disagree
        };
    }

    @Test(dataProvider = "malformedRecordPriors")
    public void testReadCallerPriorsRejectsMalformedPriors(final Genotype genotype) {
        try {
            GenotypePriorProvenance.readCallerPriors(record(Arrays.asList(REF_A, C, NON_REF), genotype), genotype);
            Assert.fail("expected malformed priors to be rejected");
        } catch (final UserException.BadInput e) {
            Assert.assertTrue(e.getMessage().contains("s at chr1:100"), e.getMessage());
            Assert.assertTrue(e.getMessage().contains("--use-input-genotype-priors"), "the error should say which options read priors: " + e.getMessage());
        }
    }

    /** A genotype as HaplotypeCaller writes it: GP is PL plus PG, scaled so the best genotype is 0. */
    private static GenotypeBuilder withPosteriors(final GenotypeBuilder builder, final int[] pl, final double[] pg) {
        final double[] gp = IntStream.range(0, pl.length).mapToDouble(i -> pl[i] + pg[i]).toArray();
        final double best = Arrays.stream(gp).min().getAsDouble();
        return builder.PL(pl).attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, pg)
                .attribute(VCFConstants.GENOTYPE_POSTERIORS_KEY, Arrays.stream(gp).map(x -> Math.round((x - best) * 100) / 100.0).toArray());
    }

    @Test
    public void testWhyPriorsDoNotFitRecord() {
        final List<Allele> alleles = Arrays.asList(REF_A, C, NON_REF);
        final int[] pl = {30, 0, 40, 33, 43, 76};
        final double[] pg = {0, 34.77, 37.78, 39.03, 73.8, 42.04};

        // as the caller wrote it, with or without GP
        final Genotype written = withPosteriors(new GenotypeBuilder("s", Arrays.asList(REF_A, C)), pl, pg).make();
        Assert.assertNull(GenotypePriorProvenance.whyPriorsDoNotFitRecord(record(alleles, written), written, pg));
        final Genotype withoutPosteriors = diploid("s").attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, pg).make();
        Assert.assertNull(GenotypePriorProvenance.whyPriorsDoNotFitRecord(record(alleles, withoutPosteriors), withoutPosteriors, pg));

        // too few priors for the record's alleles, as when a combined GVCF added an allele
        final List<Allele> merged = Arrays.asList(REF_A, G, C, NON_REF);
        final Genotype extracted = new GenotypeBuilder(written).PL(new int[] {30, 90, 90, 0, 90, 40, 33, 90, 43, 76}).make();
        Assert.assertNotNull(GenotypePriorProvenance.whyPriorsDoNotFitRecord(record(merged, extracted), extracted, pg));

        // PL remapped to reordered alleles while PG and GP stayed: right length, but GP is no longer PL plus PG
        final Genotype reordered = new GenotypeBuilder(written).PL(new int[] {30, 33, 76, 0, 43, 40}).make();
        final String misfit = GenotypePriorProvenance.whyPriorsDoNotFitRecord(record(Arrays.asList(REF_A, NON_REF, C), reordered), reordered, pg);
        Assert.assertNotNull(misfit);
        Assert.assertTrue(misfit.contains("GP"), misfit);

        // within the tolerance for rounding
        final double[] nudged = ((double[]) written.getExtendedAttribute(VCFConstants.GENOTYPE_POSTERIORS_KEY)).clone();
        nudged[2] += GenotypePriorProvenance.MAX_POSTERIOR_DISCREPANCY * 0.9;
        final Genotype rounded = new GenotypeBuilder(written).attribute(VCFConstants.GENOTYPE_POSTERIORS_KEY, nudged).make();
        Assert.assertNull(GenotypePriorProvenance.whyPriorsDoNotFitRecord(record(alleles, rounded), rounded, pg));
        nudged[2] += GenotypePriorProvenance.MAX_POSTERIOR_DISCREPANCY * 0.2;
        final Genotype off = new GenotypeBuilder(written).attribute(VCFConstants.GENOTYPE_POSTERIORS_KEY, nudged).make();
        Assert.assertNotNull(GenotypePriorProvenance.whyPriorsDoNotFitRecord(record(alleles, off), off, pg));
    }

    @Test
    public void testHasPriorsWithoutProvenance() {
        Assert.assertTrue(GenotypePriorProvenance.hasPriorsWithoutProvenance(diploid("s").attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, "0,1,2,3,4,5").make()));
        Assert.assertTrue(GenotypePriorProvenance.hasPriorsWithoutProvenance(diploid("s").attribute(GATKVCFConstants.DRAGEN_GENOTYPE_PRIOR_KEY, "0,1,2,3,4,5").make()));
        Assert.assertFalse(GenotypePriorProvenance.hasPriorsWithoutProvenance(diploid("s").attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, "0,1,2,3,4,5")
                .attribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY, "1|100|2|A/C/<NON_REF>|0.0/1.0/2.0/3.0/4.0/5.0").make()));
        Assert.assertFalse(GenotypePriorProvenance.hasPriorsWithoutProvenance(diploid("s").make()));
        Assert.assertFalse(GenotypePriorProvenance.hasPriorsWithoutProvenance(diploid("s").attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, ".").make()));
    }

    @Test
    public void testPercentEncoding() {
        Assert.assertEquals(GenotypePriorProvenance.percentEncode("<NON_REF>"), "<NON_REF>");
        Assert.assertEquals(GenotypePriorProvenance.percentEncode("A:B|C/D %\u00e9 "), "A%3AB%7CC%2FD%20%25%C3%A9%20");
        Assert.assertEquals(GenotypePriorProvenance.percentDecode("A%3AB%7CC%2FD%20%25%C3%A9%20"), "A:B|C/D %\u00e9 ");
        Assert.assertEquals(GenotypePriorProvenance.percentDecode("%3a"), ":");
    }

    @DataProvider
    public Object[][] malformedPercentEncodings() {
        return new Object[][] {{"%+F"}, {"%-1"}, {"%\uFF11\uFF21"}, {"%G0"}, {"%3"}, {"A%"}};
    }

    @Test(dataProvider = "malformedPercentEncodings", expectedExceptions = IllegalArgumentException.class)
    public void testPercentDecodeRejectsMalformedEscapes(final String value) {
        GenotypePriorProvenance.percentDecode(value);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testConstructorRejectsWrongNumberOfPriors() {
        new GenotypePriorProvenance(100, 2, Arrays.asList(REF_A, C), new double[] {0, 1});
    }

    @Test
    public void testPriorsAreCopied() {
        final double[] phred = {0, 1, 2};
        final GenotypePriorProvenance provenance = new GenotypePriorProvenance(100, 2, Arrays.asList(REF_A, C), phred);
        phred[1] = 99;
        provenance.getPhredPriors()[2] = 99;
        Assert.assertEquals(provenance.getPhredPriors(), new double[] {0, 1, 2});
    }
}
