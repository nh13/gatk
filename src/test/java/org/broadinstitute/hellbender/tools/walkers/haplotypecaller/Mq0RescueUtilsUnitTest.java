package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Unit tests for the MAPQ-0 rescue's allele counting, trigger, event selection and call checks in {@link Mq0RescueUtils}.
 */
public class Mq0RescueUtilsUnitTest extends GATKBaseTest {

    //                                       1    6    11   16   21   26   31   36
    private static final String REF_BASES = "AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT";
    private static final SAMFileHeader HEADER = ArtificialReadUtils.createArtificialSamHeader(2, 1, 1000);
    private static final String CONTIG = HEADER.getSequence(0).getSequenceName();
    private static final SimpleInterval SPAN = new SimpleInterval(CONTIG, 11, 30);
    private static final int MIN_BASE_QUALITY = 10;

    private static ReferenceContext referenceContext() {
        final SimpleInterval window = new SimpleInterval(CONTIG, 1, REF_BASES.length());
        return new ReferenceContext(new ReferenceMemorySource(new ReferenceBases(REF_BASES.getBytes(), window), HEADER.getSequenceDictionary()), window);
    }

    private static GATKRead read(final String name, final int start, final String bases, final String cigar) {
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 30);
        return ArtificialReadUtils.createArtificialRead(HEADER, name, 0, start, bases.getBytes(), quals, cigar);
    }

    private static List<GATKRead> copies(final int n, final String name, final int start, final String bases, final String cigar) {
        final List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            reads.add(read(name + i, start, bases, cigar));
        }
        return reads;
    }

    private static Mq0RescueUtils.AlleleCounts counts(final List<GATKRead> reads) {
        return new Mq0RescueUtils.AlleleCounts(reads, SPAN, referenceContext(), MIN_BASE_QUALITY);
    }

    private static Mq0RescueUtils.ReadAllele allele(final int position, final String ref, final String alt) {
        return new Mq0RescueUtils.ReadAllele(position, ref, alt);
    }

    @Test
    public void testAlleleCounts() {
        // a G>A at 13
        final GATKRead snp = read("snp", 11, "GGAGGTTTTT", "10M");
        // the same G>A, but at a base quality below the threshold
        final GATKRead lowQual = read("lowQual", 11, "GGAGGTTTTT", "10M");
        final byte[] quals = lowQual.getBaseQualities();
        quals[2] = MIN_BASE_QUALITY - 1;
        lowQual.setBaseQualities(quals);
        // a deletion of 18-19 in the T run, and one of 16-17: both left-align to a deletion of 16-17 anchored at 15
        final GATKRead deletion = read("deletion", 11, "GGGGGTTTAA", "7M2D3M");
        final GATKRead leftDeletion = read("leftDeletion", 11, "GGGGGTTTAA", "5M2D5M");
        // an extra T in the T run, which left-aligns to an insertion after 15
        final GATKRead insertion = read("insertion", 11, "GGGGGTTTTTTAA", "8M1I4M");
        // soft-clipped, with an A inserted after 15
        final GATKRead softClipped = read("softClipped", 11, "TTGGGGGATTT", "2S5M1I3M");
        // outside the span
        final GATKRead outside = read("outside", 1, "AAAAA", "5M");

        final Mq0RescueUtils.AlleleCounts counts = counts(List.of(snp, lowQual, deletion, leftDeletion, insertion, softClipped, outside));
        Assert.assertEquals(counts.count(allele(13, "G", "A")), 1);
        Assert.assertEquals(counts.count(allele(15, "GTT", "G")), 2);
        Assert.assertEquals(counts.count(allele(17, "TTT", "T")), 0);
        Assert.assertEquals(counts.count(allele(15, "G", "GT")), 1);
        Assert.assertEquals(counts.count(allele(15, "G", "GA")), 1);
        // a base below the quality threshold is not counted, a deleted base is
        final int[] expectedDepth = {6, 6, 5, 6, 6, 6, 6, 6, 5, 5, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0};
        for (int i = 0; i < expectedDepth.length; i++) {
            Assert.assertEquals(counts.depth(SPAN.getStart() + i), expectedDepth[i], "depth at " + (SPAN.getStart() + i));
        }

        // the deletion's 2 reads are enough at 15; at 13, 1 of 5 reads reaches a fraction of 0.2 exactly
        Assert.assertEquals(counts.altPositions(2, 0.1), List.of(15));
        Assert.assertEquals(counts.altPositions(1, 0.2), List.of(13, 15));
        Assert.assertEquals(counts.altPositions(1, 0.25), List.of(15));
        Assert.assertEquals(counts.altPositions(2, 0.34), Collections.emptyList());

        // two different indels at the same position are counted apart: one read each is not 2 reads with the same indel
        Assert.assertEquals(counts(List.of(insertion, softClipped)).altPositions(2, 0.0), Collections.emptyList());
        Assert.assertEquals(counts(List.of(insertion, softClipped)).altPositions(1, 0.0), List.of(15));
    }

    private static Event event(final int start, final String ref, final String alt) {
        return new Event(CONTIG, start, Allele.create(ref, true), Allele.create(alt));
    }

    private static List<List<Event>> eventSets(final Mq0RescueUtils.AlleleCounts counts, final List<Integer> positions, final int minAltReads,
                                               final int maxAlleles, final int maxEventSets) {
        return counts.eventSets(counts.selectAlleles(positions, minAltReads, maxAlleles), minAltReads, maxEventSets);
    }

    private static List<GATKRead> eventReads() {
        final List<GATKRead> reads = new ArrayList<>();
        // three reads carry G>A at 13 and T>C at 18 together
        reads.addAll(copies(3, "phased", 11, "GGAGGTTCTT", "10M"));
        // one read carries G>T at 13
        reads.add(read("other", 11, "GGTGGTTTTT", "10M"));
        // two reads carry a deletion of 16-17, anchored at 15
        reads.addAll(copies(2, "deletion", 11, "GGGGGTTT", "5M2D3M"));
        return reads;
    }

    @Test
    public void testEventSets() {
        final List<Integer> positions = List.of(13, 15, 18);
        final Event snp13 = event(13, "G", "A");
        final Event deletion = event(15, "GTT", "G");
        final Event snp18 = event(18, "T", "C");
        final Mq0RescueUtils.AlleleCounts counts = counts(eventReads());

        // one set per allele seen in at least 2 reads (not the G>T at 13), then the alleles carried together by 2 reads
        Assert.assertEquals(eventSets(counts, positions, 2, 100, 8), List.of(List.of(snp13), List.of(deletion), List.of(snp18), List.of(snp13, snp18)));
        // no more than the maximum number of multi-event sets
        Assert.assertEquals(eventSets(counts, positions, 2, 100, 0), List.of(List.of(snp13), List.of(deletion), List.of(snp18)));
        // no more than the maximum number of alleles, the most frequent kept, still in position order
        Assert.assertEquals(eventSets(counts, positions, 2, 2, 8), List.of(List.of(snp13), List.of(snp18), List.of(snp13, snp18)));
        // an allele that does not reach the minimum is not genotyped, even if it is the most frequent one at its position
        Assert.assertEquals(eventSets(counts, positions, 4, 100, 8), Collections.emptyList());
        // both alleles at a position that reach the minimum, in allele order
        Assert.assertEquals(eventSets(counts, List.of(13), 1, 100, 8), List.of(List.of(snp13), List.of(event(13, "G", "T"))));
        // only the given positions
        Assert.assertEquals(eventSets(counts, List.of(15), 2, 100, 8), List.of(List.of(deletion)));
        Assert.assertEquals(eventSets(counts, List.of(20), 2, 100, 8), Collections.emptyList());

        // an insertion is placed at its anchor base, left-aligned
        final Mq0RescueUtils.AlleleCounts insertions = counts(copies(2, "insertion", 11, "GGGGGTTTTTTAA", "8M1I4M"));
        Assert.assertEquals(eventSets(insertions, List.of(15), 2, 100, 8), List.of(List.of(event(15, "G", "GT"))));
        // an allele with a base other than A, C, G or T is not genotyped
        final Mq0RescueUtils.AlleleCounts nInsertions = counts(copies(2, "insertion", 11, "GGGGGNTTT", "5M1I3M"));
        Assert.assertEquals(nInsertions.altPositions(2, 0.1), List.of(15));
        Assert.assertEquals(eventSets(nInsertions, List.of(15), 2, 100, 8), Collections.emptyList());
    }

    @Test
    public void testEventSetsCountReadsCarryingAllAlleles() {
        // reads of different extents: 2 carry G>A at 13 and T>C at 18, 2 more carry those and A>G at 23 as well, and
        // 1 carries only T>C at 18 and A>G at 23
        final List<GATKRead> reads = new ArrayList<>();
        reads.addAll(copies(2, "short", 11, "GGAGGTTCTT", "10M"));
        reads.addAll(copies(2, "long", 11, "GGAGGTTCTTAAGAA", "15M"));
        reads.add(read("right", 16, "TTCTTAAGAA", "10M"));
        final Event snp13 = event(13, "G", "A");
        final Event snp18 = event(18, "T", "C");
        final Event snp23 = event(23, "A", "G");
        final Mq0RescueUtils.AlleleCounts counts = counts(reads);
        final List<Integer> positions = List.of(13, 18, 23);

        // a set counts every read carrying all of its alleles: 4 reads carry 13 and 18, 3 carry 18 and 23, 2 carry all 3
        Assert.assertEquals(eventSets(counts, positions, 3, 100, 8),
                List.of(List.of(snp13), List.of(snp18), List.of(snp23), List.of(snp13, snp18), List.of(snp18, snp23)));
        // larger sets first, so that alleles carried together are phased as fully as the reads allow
        Assert.assertEquals(eventSets(counts, positions, 2, 100, 8),
                List.of(List.of(snp13), List.of(snp18), List.of(snp23), List.of(snp13, snp18, snp23), List.of(snp13, snp18), List.of(snp18, snp23)));
        Assert.assertEquals(eventSets(counts, positions, 2, 100, 1), List.of(List.of(snp13), List.of(snp18), List.of(snp23), List.of(snp13, snp18, snp23)));

        // without the read carrying only 18 and 23, no read carries those two without 13
        final Mq0RescueUtils.AlleleCounts nested = counts(reads.subList(0, 4));
        Assert.assertEquals(eventSets(nested, positions, 2, 100, 8),
                List.of(List.of(snp13), List.of(snp18), List.of(snp23), List.of(snp13, snp18, snp23), List.of(snp13, snp18)));
    }

    private static final Allele REF_A = Allele.create("A", true);
    private static final Allele ALT_C = Allele.create("C");

    private static VariantContext call(final int pos, final double qual, final Allele... genotype) {
        return new VariantContextBuilder("test", CONTIG, pos, pos, List.of(REF_A, ALT_C, Allele.NON_REF_ALLELE))
                .log10PError(qual / -10.0)
                .genotypes(new GenotypeBuilder("sample", Arrays.asList(genotype)).make())
                .make();
    }

    private static VariantContext refBlock(final int pos) {
        return new VariantContextBuilder("test", CONTIG, pos, pos, List.of(REF_A, Allele.NON_REF_ALLELE))
                .genotypes(new GenotypeBuilder("sample", List.of(REF_A, REF_A)).make())
                .make();
    }

    private static VariantContext record(final int pos, final String ref, final List<String> alts, final int... genotype) {
        final List<Allele> alleles = new ArrayList<>();
        alleles.add(Allele.create(ref, true));
        alts.forEach(alt -> alleles.add(alt.equals("<NON_REF>") ? Allele.NON_REF_ALLELE : Allele.create(alt)));
        final List<Allele> gt = new ArrayList<>();
        for (final int i : genotype) {
            gt.add(i < 0 ? Allele.NO_CALL : alleles.get(i));
        }
        return new VariantContextBuilder("test", CONTIG, pos, pos + ref.length() - 1, alleles)
                .log10PError(-5).genotypes(new GenotypeBuilder("sample", gt).make()).make();
    }

    @Test
    public void testSupportsCalledAlleles() {
        final Mq0RescueUtils.AlleleCounts counts = counts(eventReads());
        Assert.assertTrue(counts.supportsCalledAlleles(record(13, "G", List.of("A", "<NON_REF>"), 0, 1)));
        Assert.assertTrue(counts.supportsCalledAlleles(record(13, "G", List.of("A"), 1, 1)));
        // an allele no MAPQ-0 read carries
        Assert.assertFalse(counts.supportsCalledAlleles(record(13, "G", List.of("C"), 0, 1)));
        // no called alternate allele
        Assert.assertFalse(counts.supportsCalledAlleles(record(13, "G", List.of("A"), 0, 0)));
        Assert.assertFalse(counts.supportsCalledAlleles(record(13, "G", List.of("A"), -1, -1)));
        // one called alternate allele that no read carries is enough to fail, an uncalled one is not
        Assert.assertFalse(counts.supportsCalledAlleles(record(13, "G", List.of("A", "C"), 1, 2)));
        Assert.assertTrue(counts.supportsCalledAlleles(record(13, "G", List.of("A", "C"), 0, 1)));
        // the deletion, with a padded reference allele, or placed further right in the T run
        Assert.assertTrue(counts.supportsCalledAlleles(record(15, "GTTT", List.of("GT", "<NON_REF>"), 0, 1)));
        Assert.assertTrue(counts.supportsCalledAlleles(record(16, "TTT", List.of("T"), 0, 1)));
        Assert.assertFalse(counts.supportsCalledAlleles(record(15, "GT", List.of("G"), 0, 1)));
        // padded multi-allelic records, whose called allele is a suffix of the reference allele
        Assert.assertTrue(counts.supportsCalledAlleles(record(15, "GTTTT", List.of("GTT", "G"), 0, 1)));
        Assert.assertFalse(counts.supportsCalledAlleles(record(13, "GGG", List.of("GG", "G"), 0, 1)));
        // an MNP, as --max-mnp-distance makes them, is carried if each of its substitutions is
        Assert.assertTrue(counts.supportsCalledAlleles(record(13, "GGGTTT", List.of("AGGTTC"), 0, 1)));
        Assert.assertFalse(counts.supportsCalledAlleles(record(13, "GGGTTT", List.of("AGGTTG"), 0, 1)));
        // the SNP at 18 in a record padded by a deletion at 18
        Assert.assertTrue(counts.supportsCalledAlleles(record(18, "TT", List.of("T", "CT"), 0, 2)));
    }

    @Test
    public void testIndelLeftAlignedPastOwnAllele() {
        // reads carrying T>C at 17 and a deletion aligned at 19, in the T run: left-aligned against the reference, the
        // deletion is anchored at 15, so the two alleles together would make a haplotype (GGGGGCTTTAA) the reads do not
        // carry (GGGGGTCTTAA). They count for each allele but not for multi-event sets.
        final List<GATKRead> shifted = copies(2, "shifted", 11, "GGGGGTCTTAA", "8M1D3M");
        final Mq0RescueUtils.AlleleCounts counts = counts(shifted);
        Assert.assertEquals(counts.count(allele(17, "T", "C")), 2);
        Assert.assertEquals(counts.count(allele(15, "GT", "G")), 2);
        Assert.assertEquals(eventSets(counts, List.of(15, 17), 2, 100, 8), List.of(List.of(event(15, "GT", "G")), List.of(event(17, "T", "C"))));
        // the same for an insertion aligned after the read's own SNV: an extra T after 18, left-aligned to after 15
        final Mq0RescueUtils.AlleleCounts insertion = counts(copies(2, "shiftedInsertion", 11, "GGGGGTCTTTTAA", "8M1I4M"));
        Assert.assertEquals(insertion.count(allele(15, "G", "GT")), 2);
        Assert.assertEquals(eventSets(insertion, List.of(15, 17), 2, 100, 8), List.of(List.of(event(15, "G", "GT")), List.of(event(17, "T", "C"))));
        // a SNV that the shift does not pass, at 13, is phased with the deletion
        final Mq0RescueUtils.AlleleCounts phased = counts(copies(2, "phased", 11, "GGAGGTTTTAA", "8M1D3M"));
        Assert.assertEquals(eventSets(phased, List.of(13, 15), 2, 100, 8),
                List.of(List.of(event(13, "G", "A")), List.of(event(15, "GT", "G")), List.of(event(13, "G", "A"), event(15, "GT", "G"))));
    }

    @Test
    public void testSupportOutsideSpan() {
        // an MNP starting at the last base of the span (30), whose second base (31) is outside it
        final Mq0RescueUtils.AlleleCounts mnp = counts(copies(2, "mnp", 21, "AAAAACCCCATGGGG", "15M"));
        Assert.assertTrue(mnp.supportsCalledAlleles(record(30, "CG", List.of("AT"), 0, 1)));
        Assert.assertFalse(mnp.supportsCalledAlleles(record(30, "CG", List.of("AC"), 0, 1)));
        // only the span counts for the trigger
        Assert.assertEquals(mnp.altPositions(2, 0.1), List.of(30));
        // a deletion called in the span that left-aligns to an anchor before it (10)
        final Mq0RescueUtils.AlleleCounts deletion = counts(copies(1, "deletion", 6, "CCCCCGGGGTTTTT", "7M1D7M"));
        Assert.assertEquals(deletion.count(allele(10, "CG", "C")), 1);
        Assert.assertTrue(deletion.supportsCalledAlleles(record(12, "GGG", List.of("GG"), 0, 1)));
        Assert.assertEquals(deletion.altPositions(1, 0.0), Collections.emptyList());
    }

    @Test
    public void testSupportsCalledInsertions() {
        // an extra T in the T run, aligned after 18, left-aligned to an insertion after 15
        final Mq0RescueUtils.AlleleCounts counts = counts(copies(2, "insertion", 11, "GGGGGTTTTTTAA", "8M1I4M"));
        Assert.assertTrue(counts.supportsCalledAlleles(record(15, "G", List.of("GT"), 0, 1)));
        // the same insertion placed further right in the T run, or with a padded reference allele
        Assert.assertTrue(counts.supportsCalledAlleles(record(19, "T", List.of("TT"), 0, 1)));
        Assert.assertTrue(counts.supportsCalledAlleles(record(15, "GT", List.of("GTT"), 0, 1)));
        // a different inserted base
        Assert.assertFalse(counts.supportsCalledAlleles(record(15, "G", List.of("GA"), 0, 1)));
    }

    @Test
    public void testSupportsCompoundAlleles() {
        // reads carrying a SNP (G>C at 15) and a deletion of 16 on the same anchor, as assembly writes them: GT to C
        final Mq0RescueUtils.AlleleCounts deletion = counts(copies(2, "snpDeletion", 11, "GGGGCTTTT", "5M1D4M"));
        Assert.assertTrue(deletion.supportsCalledAlleles(record(15, "GT", List.of("C"), 0, 1)));
        // with a padded reference allele, padded at the end or, in a multi-allelic record, at the start
        Assert.assertTrue(deletion.supportsCalledAlleles(record(15, "GTT", List.of("CT", "<NON_REF>"), 0, 1)));
        Assert.assertTrue(deletion.supportsCalledAlleles(record(14, "GGT", List.of("G", "GC"), 0, 2)));
        // the plain deletion is carried too, but a SNP the reads do not carry is not
        Assert.assertTrue(deletion.supportsCalledAlleles(record(15, "GT", List.of("G"), 0, 1)));
        Assert.assertFalse(deletion.supportsCalledAlleles(record(15, "GT", List.of("A"), 0, 1)));
        // a compound allele whose substitutions the reads do not all carry (T>A at 17 in each placement of the gap)
        Assert.assertFalse(deletion.supportsCalledAlleles(record(15, "GTT", List.of("CA"), 0, 1)));

        // the same haplotype aligned with the gap first: a deletion of 15 (left-aligned to after 10) and T>C at 16
        final Mq0RescueUtils.AlleleCounts gapFirst = counts(copies(2, "deletionSnp", 11, "GGGGCTTT", "4M1D3M"));
        Assert.assertTrue(gapFirst.supportsCalledAlleles(record(15, "GT", List.of("C"), 0, 1)));

        // a SNP and an insertion on the same anchor: G to CT, a G>C at 15 and a T inserted after it, or, aligned with
        // the gap first, a C inserted after 14 and G>T at 15
        final Mq0RescueUtils.AlleleCounts insertion = counts(copies(2, "snpInsertion", 11, "GGGGCTTTTTT", "5M1I5M"));
        Assert.assertTrue(insertion.supportsCalledAlleles(record(15, "G", List.of("CT"), 0, 1)));
        Assert.assertFalse(insertion.supportsCalledAlleles(record(15, "G", List.of("CA"), 0, 1)));
        final Mq0RescueUtils.AlleleCounts insertionFirst = counts(copies(2, "insertionSnp", 11, "GGGGCTTTTT", "4M1I5M"));
        Assert.assertTrue(insertionFirst.supportsCalledAlleles(record(15, "G", List.of("CT"), 0, 1)));

        // only called alleles count as compound
        Assert.assertTrue(deletion.hasCompoundAllele(record(15, "GT", List.of("C"), 0, 1)));
        Assert.assertTrue(deletion.hasCompoundAllele(record(15, "GTT", List.of("CA"), 0, 1)));
        Assert.assertTrue(deletion.hasCompoundAllele(record(14, "GGT", List.of("G", "GC"), 0, 2)));
        Assert.assertFalse(deletion.hasCompoundAllele(record(15, "GT", List.of("G"), 0, 1)));
        Assert.assertFalse(deletion.hasCompoundAllele(record(15, "GT", List.of("G", "C"), 0, 1)));
    }

    @Test
    public void testIsCall() {
        Assert.assertTrue(Mq0RescueUtils.isCall(call(100, 50, REF_A, ALT_C), 3.0));
        Assert.assertTrue(Mq0RescueUtils.isCall(call(100, 50, ALT_C, ALT_C), 3.0));
        Assert.assertTrue(Mq0RescueUtils.isCall(call(100, 3.0, REF_A, ALT_C), 3.0));
        Assert.assertFalse(Mq0RescueUtils.isCall(call(100, 2.9, REF_A, ALT_C), 3.0));
        Assert.assertFalse(Mq0RescueUtils.isCall(call(100, 50, REF_A, REF_A), 3.0));
        Assert.assertFalse(Mq0RescueUtils.isCall(call(100, 50, Allele.NO_CALL, Allele.NO_CALL), 3.0));
        Assert.assertFalse(Mq0RescueUtils.isCall(refBlock(100), 0.0));
        // a genotype whose only non-reference allele is <NON_REF>, or that is half no-call, calls no real allele
        Assert.assertFalse(Mq0RescueUtils.isCall(call(100, 50, REF_A, Allele.NON_REF_ALLELE), 3.0));
        Assert.assertFalse(Mq0RescueUtils.isCall(call(100, 50, REF_A, Allele.NO_CALL), 3.0));
        Assert.assertTrue(Mq0RescueUtils.isCall(call(100, 50, ALT_C, Allele.NO_CALL), 3.0));
        Assert.assertTrue(Mq0RescueUtils.hasRealAlt(call(100, 0, REF_A, REF_A)));
        Assert.assertFalse(Mq0RescueUtils.hasRealAlt(refBlock(100)));
    }

    @Test
    public void testAwayFrom() {
        final List<VariantContext> calls = List.of(call(205, 50, REF_A, ALT_C), record(300, "AAAAAAAAAAAAAAAAAAAAA", List.of("A"), 0, 1));
        // exactly 10 bp away is too close, 11 is not; a deletion's whole span counts
        Assert.assertEquals(Mq0RescueUtils.awayFrom(List.of(194, 195, 205, 215, 216, 289, 290, 330, 331), CONTIG, calls, 10),
                List.of(194, 216, 289, 331));
        Assert.assertEquals(Mq0RescueUtils.awayFrom(List.of(205), HEADER.getSequence(1).getSequenceName(), calls, 10), List.of(205));
    }

    @Test
    public void testCountOverlapping() {
        final GATKRead read = read("read", 100, "ACGTACGTAC", "10M");
        final GATKRead otherContig = ArtificialReadUtils.createArtificialRead(HEADER, "other", 1, 100, "ACGTACGTAC".getBytes(),
                new byte[10], "10M");
        final List<GATKRead> reads = List.of(read, otherContig);
        Assert.assertEquals(Mq0RescueUtils.countOverlapping(reads, CONTIG, 99), 0);
        Assert.assertEquals(Mq0RescueUtils.countOverlapping(reads, CONTIG, 100), 1);
        Assert.assertEquals(Mq0RescueUtils.countOverlapping(reads, CONTIG, 109), 1);
        Assert.assertEquals(Mq0RescueUtils.countOverlapping(reads, CONTIG, 110), 0);
        Assert.assertEquals(Mq0RescueUtils.countOverlapping(reads, HEADER.getSequence(1).getSequenceName(), 105), 1);
    }
}
