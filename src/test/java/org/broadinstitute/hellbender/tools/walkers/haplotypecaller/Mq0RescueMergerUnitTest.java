package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link Mq0RescueMerger}, which merges HaplotypeCaller's MAPQ-0 rescue records into its regular output.
 */
public class Mq0RescueMergerUnitTest extends GATKBaseTest {

    private static final SAMFileHeader HEADER = ArtificialReadUtils.createArtificialSamHeader(2, 1, 1_000_000);
    private static final String CONTIG_1 = HEADER.getSequence(0).getSequenceName();
    private static final String CONTIG_2 = HEADER.getSequence(1).getSequenceName();
    private static final int DISTANCE = 10;
    private static final double MIN_QUAL = 3.0;

    private static Mq0RescueMerger merger() {
        return new Mq0RescueMerger(HEADER.getSequenceDictionary(), DISTANCE, MIN_QUAL);
    }

    private static VariantContext refBlock(final String contig, final int pos) {
        final Allele ref = Allele.create("A", true);
        return new VariantContextBuilder("test", contig, pos, pos, List.of(ref, Allele.NON_REF_ALLELE))
                .genotypes(new GenotypeBuilder("sample", List.of(ref, ref)).make()).make();
    }

    /** a record with a real alternate allele: a SNP, or a deletion if {@code length} > 1, het (or hom-ref) */
    private static VariantContext variant(final String contig, final int pos, final int length, final double qual, final boolean homRef) {
        final Allele ref = Allele.create("A".repeat(length), true);
        final Allele alt = Allele.create(length > 1 ? "A" : "C");
        return new VariantContextBuilder("test", contig, pos, pos + length - 1, List.of(ref, alt, Allele.NON_REF_ALLELE))
                .log10PError(qual / -10.0)
                .genotypes(new GenotypeBuilder("sample", List.of(ref, homRef ? ref : alt)).make()).make();
    }

    /** a record with only a reference allele, as --output-mode EMIT_ALL_ACTIVE_SITES writes them */
    private static VariantContext monomorphic(final String contig, final int pos, final int length) {
        final Allele ref = Allele.create("A".repeat(length), true);
        return new VariantContextBuilder("test", contig, pos, pos + length - 1, List.of(ref))
                .genotypes(new GenotypeBuilder("sample", List.of(ref, ref)).make()).make();
    }

    private static VariantContext rescue(final String contig, final int pos, final int length) {
        return new VariantContextBuilder(variant(contig, pos, length, 50, false))
                .filter(GATKVCFConstants.MQ0_RESCUE_FILTER_NAME).attribute(GATKVCFConstants.MQ0_RESCUE_KEY, true).make();
    }

    /** the regular records of a GVCF region: one reference record per base, with the calls in place of some */
    private static List<VariantContext> gvcfRegion(final String contig, final int from, final int to, final VariantContext... calls) {
        final List<VariantContext> records = new ArrayList<>();
        for (int pos = from; pos <= to; pos++) {
            final int p = pos;
            final VariantContext call = List.of(calls).stream().filter(vc -> vc.getStart() == p).findFirst().orElse(null);
            records.add(call != null ? call : refBlock(contig, pos));
        }
        return records;
    }

    private static boolean isRescue(final VariantContext vc) {
        return vc.hasAttribute(GATKVCFConstants.MQ0_RESCUE_KEY);
    }

    private static String describe(final List<VariantContext> records) {
        return records.stream().map(vc -> vc.getContig() + ":" + vc.getStart() + (isRescue(vc) ? "R" :
                Mq0RescueMerger.isReferenceConfidence(vc) ? "" : "V")).collect(Collectors.joining(" "));
    }

    private static List<Integer> rescueStarts(final List<VariantContext> records) {
        return records.stream().filter(Mq0RescueMergerUnitTest::isRescue).map(VariantContext::getStart).collect(Collectors.toList());
    }

    private static void assertSortedAndNotOverlappingRescue(final List<VariantContext> records) {
        for (int i = 1; i < records.size(); i++) {
            final VariantContext a = records.get(i - 1);
            final VariantContext b = records.get(i);
            Assert.assertTrue(HEADER.getSequenceIndex(a.getContig()) < HEADER.getSequenceIndex(b.getContig())
                    || (a.getContig().equals(b.getContig()) && a.getStart() < b.getStart()), "out of order: " + describe(records));
        }
        for (final VariantContext r : records) {
            if (isRescue(r)) {
                Assert.assertTrue(records.stream().noneMatch(vc -> vc != r && vc.overlaps(r)), r + " overlaps another record");
            }
        }
    }

    private static List<VariantContext> run(final Mq0RescueMerger merger, final List<Object[]> regions) {
        final List<VariantContext> out = new ArrayList<>();
        for (final Object[] region : regions) {
            @SuppressWarnings("unchecked") final List<VariantContext> regular = (List<VariantContext>) region[1];
            @SuppressWarnings("unchecked") final List<VariantContext> rescue = (List<VariantContext>) region[2];
            out.addAll(merger.add((SimpleInterval) region[0], regular, rescue));
        }
        out.addAll(merger.finish());
        assertSortedAndNotOverlappingRescue(out);
        return out;
    }

    private static Object[] region(final String contig, final int from, final int to, final List<VariantContext> regular, final VariantContext... rescue) {
        return new Object[]{new SimpleInterval(contig, from, to), regular, List.of(rescue)};
    }

    @Test
    public void testReplacesReferenceBlocks() {
        // a SNP takes the place of one reference record, a deletion of all those under it
        final VariantContext snp = rescue(CONTIG_1, 103, 1);
        final VariantContext deletion = rescue(CONTIG_1, 3000, 4);
        final List<VariantContext> out = run(merger(), List.of(
                region(CONTIG_1, 100, 110, gvcfRegion(CONTIG_1, 100, 110), snp),
                region(CONTIG_1, 2990, 3010, gvcfRegion(CONTIG_1, 2990, 3010), deletion)));
        Assert.assertEquals(out.size(), 11 + 21 - 4 + 1);
        Assert.assertSame(out.get(3), snp);
        Assert.assertEquals(out.get(4).getStart(), 104);
        Assert.assertSame(out.get(21), deletion);
        Assert.assertEquals(out.get(20).getStart(), 2999);
        Assert.assertEquals(out.get(22).getStart(), 3004);
    }

    @Test
    public void testRescueCoversLastPendingRecord() {
        // the rescue record covers the last records still pending, and is decided by finish()
        final VariantContext atEnd = rescue(CONTIG_1, 110, 1);
        List<VariantContext> out = run(merger(), List.<Object[]>of(region(CONTIG_1, 100, 110, gvcfRegion(CONTIG_1, 100, 110), atEnd)));
        Assert.assertEquals(out.size(), 11);
        Assert.assertSame(out.get(10), atEnd);
        final VariantContext deletionAtEnd = rescue(CONTIG_1, 108, 3);
        out = run(merger(), List.<Object[]>of(region(CONTIG_1, 100, 110, gvcfRegion(CONTIG_1, 100, 110), deletionAtEnd)));
        Assert.assertEquals(out.size(), 9);
        Assert.assertSame(out.get(8), deletionAtEnd);

        // the same, decided by add() once a later region on the same contig comes in, with nothing pending after it
        final Mq0RescueMerger merger = merger();
        final List<VariantContext> first = merger.add(new SimpleInterval(CONTIG_1, 100, 110), gvcfRegion(CONTIG_1, 100, 110), List.of(atEnd));
        Assert.assertEquals(first.size(), 10);
        final List<VariantContext> second = merger.add(new SimpleInterval(CONTIG_1, 200, 210), Collections.emptyList(), Collections.emptyList());
        Assert.assertEquals(second, List.of(atEnd));
        Assert.assertTrue(merger.finish().isEmpty());
    }

    @Test
    public void testRegularCallWins() {
        // a regular call wins over a rescue record it overlaps, whether at the same position or as a deletion over it
        final VariantContext call = variant(CONTIG_1, 103, 1, 50, false);
        final VariantContext deletion = variant(CONTIG_1, 205, 5, 50, false);
        final Mq0RescueMerger merger = merger();
        final List<VariantContext> out = run(merger, List.of(
                region(CONTIG_1, 100, 110, gvcfRegion(CONTIG_1, 100, 110, call), rescue(CONTIG_1, 103, 1)),
                region(CONTIG_1, 200, 210, gvcfRegion(CONTIG_1, 200, 210, deletion), rescue(CONTIG_1, 208, 1))));
        Assert.assertEquals(describe(out), describe(List.of(gvcfRegion(CONTIG_1, 100, 110, call), gvcfRegion(CONTIG_1, 200, 210, deletion)).stream()
                .flatMap(List::stream).collect(Collectors.toList())));
        Assert.assertSame(out.get(3), call);
        Assert.assertEquals(merger.droppedOverlap(), 2);
        Assert.assertEquals(merger.merged(), 0);
    }

    @Test
    public void testNonCallRecordsAreReplaced() {
        // A regular record that is not a call (a GVCF hom-ref record, or one with QUAL below the minimum) does not keep
        // a rescue record out, as it would not be in a VCF: the rescue record takes its place, the same allele or not.
        // One that starts before the rescue record may already be written, so it does keep the rescue record out.
        final VariantContext homRef = variant(CONTIG_1, 103, 1, 50, true);
        final VariantContext lowQual = variant(CONTIG_1, 106, 1, 2, false);
        final VariantContext lowQualDeletion = variant(CONTIG_1, 204, 4, 2, false);
        final Mq0RescueMerger merger = merger();
        final List<VariantContext> out = run(merger, List.of(
                region(CONTIG_1, 100, 110, gvcfRegion(CONTIG_1, 100, 110, homRef, lowQual), rescue(CONTIG_1, 103, 1), rescue(CONTIG_1, 105, 3)),
                region(CONTIG_1, 200, 210, gvcfRegion(CONTIG_1, 200, 210, lowQualDeletion), rescue(CONTIG_1, 206, 1))));
        Assert.assertEquals(rescueStarts(out), List.of(103, 105));
        Assert.assertTrue(out.stream().noneMatch(vc -> vc == homRef || vc == lowQual));
        Assert.assertTrue(out.contains(lowQualDeletion));
        Assert.assertEquals(merger.replaced(), 2);
        Assert.assertEquals(merger.droppedOverlap(), 1);
    }

    @Test
    public void testMonomorphicRecords() {
        // --output-mode EMIT_ALL_ACTIVE_SITES writes records with only a reference allele, which are not calls: one
        // under a rescue record is replaced, one starting before a rescue record that it overlaps keeps it out
        final VariantContext under = monomorphic(CONTIG_1, 120, 1);
        final VariantContext before = monomorphic(CONTIG_1, 140, 5);
        final Mq0RescueMerger merger = merger();
        final List<VariantContext> out = run(merger, List.<Object[]>of(
                region(CONTIG_1, 100, 200, List.of(under, before), rescue(CONTIG_1, 120, 1), rescue(CONTIG_1, 142, 1))));
        Assert.assertEquals(describe(out), CONTIG_1 + ":120R " + CONTIG_1 + ":140V");
        Assert.assertEquals(merger.droppedOverlap(), 1);
    }

    @Test
    public void testDistance() {
        // a regular call exactly the rescue distance away drops a rescue record, one a base further does not; a hom-ref
        // or low-quality regular record nearby does not count as a call
        final Mq0RescueMerger merger = merger();
        final List<VariantContext> out = run(merger, List.of(
                region(CONTIG_1, 100, 200, gvcfRegion(CONTIG_1, 100, 200, variant(CONTIG_1, 112, 1, 50, false),
                                variant(CONTIG_1, 152, 1, 2, false), variant(CONTIG_1, 162, 1, 50, true)),
                        rescue(CONTIG_1, 102, 1), rescue(CONTIG_1, 123, 1), rescue(CONTIG_1, 150, 1), rescue(CONTIG_1, 160, 1)),
                region(CONTIG_1, 201, 300, gvcfRegion(CONTIG_1, 201, 300), rescue(CONTIG_1, 290, 1))));
        Assert.assertEquals(rescueStarts(out), List.of(123, 150, 160, 290));
        Assert.assertEquals(merger.droppedDistance(), 1);

        final Mq0RescueMerger edge = merger();
        final List<VariantContext> edgeOut = run(edge, List.<Object[]>of(region(CONTIG_1, 100, 200, gvcfRegion(CONTIG_1, 100, 200, variant(CONTIG_1, 150, 1, 50, false)),
                rescue(CONTIG_1, 139, 1), rescue(CONTIG_1, 140, 1))));
        // 11 bp away is kept, 10 bp is dropped
        Assert.assertEquals(rescueStarts(edgeOut), List.of(139));
        Assert.assertEquals(edge.droppedDistance(), 1);
        Assert.assertEquals(edge.merged(), 1);
    }

    @Test
    public void testDistanceToLaterRegion() {
        // the regular call that drops a rescue record can come in a later region
        final Mq0RescueMerger merger = merger();
        final List<VariantContext> out = run(merger, List.of(
                region(CONTIG_1, 100, 200, gvcfRegion(CONTIG_1, 100, 200), rescue(CONTIG_1, 195, 1), rescue(CONTIG_1, 200, 1)),
                region(CONTIG_1, 201, 300, gvcfRegion(CONTIG_1, 201, 300, variant(CONTIG_1, 206, 1, 50, false)))));
        Assert.assertEquals(rescueStarts(out), List.of(195));
        Assert.assertEquals(merger.droppedDistance(), 1);
    }

    @Test
    public void testOverlappingRescueRecords() {
        // a rescue record overlapping an earlier one is dropped, from the same region or the next
        final Mq0RescueMerger merger = merger();
        final List<VariantContext> out = run(merger, List.of(
                region(CONTIG_1, 100, 110, gvcfRegion(CONTIG_1, 100, 110), rescue(CONTIG_1, 105, 3), rescue(CONTIG_1, 106, 1), rescue(CONTIG_1, 109, 4)),
                region(CONTIG_1, 111, 120, gvcfRegion(CONTIG_1, 111, 120), rescue(CONTIG_1, 111, 1), rescue(CONTIG_1, 115, 1))));
        Assert.assertEquals(merger.droppedOverlap(), 2);
        Assert.assertEquals(rescueStarts(out), List.of(105, 109, 115));
        // the deletion at 109-112 replaced the reference records of both regions under it
        Assert.assertEquals(out.size(), 21 - 3 + 1 - 4 + 1 - 1 + 1);
    }

    @Test
    public void testHoldsRecordsUntilDecided() {
        // the records from a held rescue record on are returned only once the regular records within the rescue
        // distance of it are known
        final Mq0RescueMerger merger = merger();
        final VariantContext snp = rescue(CONTIG_1, 105, 1);
        final List<VariantContext> first = merger.add(new SimpleInterval(CONTIG_1, 100, 110), gvcfRegion(CONTIG_1, 100, 110), List.of(snp));
        Assert.assertEquals(describe(first), describe(gvcfRegion(CONTIG_1, 100, 104)));
        final List<VariantContext> second = merger.add(new SimpleInterval(CONTIG_1, 111, 114), gvcfRegion(CONTIG_1, 111, 114), Collections.emptyList());
        Assert.assertTrue(second.isEmpty());
        final List<VariantContext> third = merger.add(new SimpleInterval(CONTIG_1, 115, 120), gvcfRegion(CONTIG_1, 115, 120), Collections.emptyList());
        Assert.assertSame(third.get(0), snp);
        Assert.assertEquals(third.size(), 120 - 105 + 1);
        Assert.assertTrue(merger.finish().isEmpty());
    }

    @Test
    public void testContigChange() {
        // a rescue record held near the end of one contig is decided, and written before the next contig's records
        final List<VariantContext> out = run(merger(), List.of(
                region(CONTIG_1, 100, 110, gvcfRegion(CONTIG_1, 100, 110), rescue(CONTIG_1, 108, 1)),
                region(CONTIG_2, 100, 110, gvcfRegion(CONTIG_2, 100, 110, variant(CONTIG_2, 108, 1, 50, false)), rescue(CONTIG_2, 104, 1))));
        Assert.assertEquals(describe(out.subList(7, 12)), CONTIG_1 + ":107 " + CONTIG_1 + ":108R " + CONTIG_1 + ":109 " + CONTIG_1 + ":110 " + CONTIG_2 + ":100");
        // the call on the second contig is within 10 bp of its rescue record, but not of the first contig's
        Assert.assertEquals(out.stream().filter(Mq0RescueMergerUnitTest::isRescue).count(), 1);
    }

    @Test
    public void testVcfMode() {
        // without reference records, a rescue record goes between the regular calls
        final VariantContext a = variant(CONTIG_1, 100, 1, 50, false);
        final VariantContext b = variant(CONTIG_1, 300, 1, 50, false);
        final VariantContext r = rescue(CONTIG_1, 200, 1);
        final List<VariantContext> out = run(merger(), List.of(region(CONTIG_1, 50, 250, List.of(a), r), region(CONTIG_1, 251, 400, List.of(b))));
        Assert.assertEquals(out, List.of(a, r, b));
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testRescueOutsideRegion() {
        merger().add(new SimpleInterval(CONTIG_1, 100, 110), Collections.emptyList(), List.of(rescue(CONTIG_1, 111, 1)));
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testRegionsOutOfOrder() {
        final Mq0RescueMerger merger = merger();
        merger.add(new SimpleInterval(CONTIG_1, 200, 210), Collections.emptyList(), Collections.emptyList());
        merger.add(new SimpleInterval(CONTIG_1, 100, 110), Collections.emptyList(), Collections.emptyList());
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testRegularRecordsOutOfOrder() {
        merger().add(new SimpleInterval(CONTIG_1, 100, 110), List.of(variant(CONTIG_1, 105, 1, 50, false), variant(CONTIG_1, 101, 1, 50, false)),
                Collections.emptyList());
    }
}
