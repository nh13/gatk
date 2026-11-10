package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.Random;

/**
 * Checks {@link ReferenceConfidenceModel#sweepRead}, which walks a read's cigar an element at a time, against a
 * reference walk of the read a position at a time with the {@link AlignmentStateMachine} a pileup uses, on random
 * reads with clips, insertions, deletions, reference skips, low base qualities and short fragments. The two must add
 * the same evidence, in the same order, to every position, and reject the same malformed reads. The reference walk
 * uses the model's own per-element helpers, as the walk {@code sweepRead} replaced did, so it pins the walk itself.
 */
public final class ReferenceConfidenceModelSweepUnitTest extends GATKBaseTest {
    private static final String SAMPLE = "sample";
    private static final int CONTIG_LENGTH = 2000;
    private static final int PADDED_START = 400;
    private static final int PADDED_SIZE = 400;
    private static final byte DELETION_QUAL = 30;
    private static final int MAX_INDEL_SIZE = 10;
    private static final byte[] BASES = {'A', 'C', 'G', 'T'};
    private static final CigarOperator[] ZERO_LENGTH_OPERATORS = {CigarOperator.M, CigarOperator.I, CigarOperator.D, CigarOperator.N};

    private final SAMFileHeader header = makeHeader();

    private static SAMFileHeader makeHeader() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, CONTIG_LENGTH);
        final SAMReadGroupRecord rg = new SAMReadGroupRecord("rg");
        rg.setSample(SAMPLE);
        header.addReadGroup(rg);
        return header;
    }

    @DataProvider(name = "modes")
    public Object[][] modes() {
        return new Object[][]{{true, false}, {false, false}, {true, true}, {false, true}};
    }

    @Test(dataProvider = "modes")
    public void sweepingByCigarElementMatchesSteppingThroughThePileupStateMachine(final boolean useSoftClippedBases, final boolean flowBased) {
        final ReferenceConfidenceModel model = new ReferenceConfidenceModel(SampleList.singletonSampleList(SAMPLE), header,
                MAX_INDEL_SIZE, -1, DELETION_QUAL, useSoftClippedBases, flowBased);
        final Random random = new Random(useSoftClippedBases ? (flowBased ? 1 : 2) : (flowBased ? 3 : 4));
        final byte[] ref = randomBases(random, PADDED_SIZE);
        final int ploidy = 2;
        final double[][][] increments = ReferenceConfidenceModel.likelihoodIncrements(ploidy);

        int malformed = 0;
        int ranOffBases = 0;
        for (int region = 0; region < 300; region++) {
            final int globalRefOffset = random.nextInt(60);
            final int spanStart = PADDED_START + globalRefOffset;
            final int spanSize = 1 + random.nextInt(PADDED_SIZE - globalRefOffset - 60);
            final ReferenceConfidenceModel.SiteEvidence actual = new ReferenceConfidenceModel.SiteEvidence(spanSize, ploidy + 1);
            final Accumulators expected = new Accumulators(spanSize, ploidy);
            final int nReads = 1 + random.nextInt(60);
            for (int r = 0; r < nReads; r++) {
                final GATKRead read = randomRead(random, "read" + r, spanStart, spanSize, useSoftClippedBases);
                final RuntimeException failure = sweepBoth(model, read, spanStart, spanSize, ref, globalRefOffset, increments,
                        actual, expected, useSoftClippedBases, flowBased);
                if (failure instanceof UserException.MalformedRead) {
                    malformed++;
                } else if (failure != null) {
                    ranOffBases++;
                }
            }
            expected.assertSameAs(actual);
        }
        Assert.assertTrue(malformed > 0, "no malformed reads were generated");
        if (!flowBased) {
            Assert.assertEquals(ranOffBases, 0);
        }
    }

    /**
     * Sweeps a read with both walks, asserts that they fail alike, and if they fail starts both accumulators again
     * from the same state.
     *
     * @return the reference walk's failure, or null
     */
    private static RuntimeException sweepBoth(final ReferenceConfidenceModel model, final GATKRead read, final int spanStart, final int spanSize,
                                              final byte[] ref, final int globalRefOffset, final double[][][] increments,
                                              final ReferenceConfidenceModel.SiteEvidence actual, final Accumulators expected,
                                              final boolean useSoftClippedBases, final boolean flowBased) {
        final RuntimeException actualFailure = sweep(() -> model.sweepRead(read, spanStart, spanSize, ref, globalRefOffset,
                increments, actual));
        final RuntimeException expectedFailure = sweep(() -> referenceSweep(model, read, spanStart, spanSize, ref, globalRefOffset,
                increments, expected, useSoftClippedBases, flowBased));
        final String cigar = read.getCigar().toString();
        Assert.assertEquals(actualFailure == null ? null : actualFailure.getClass(),
                expectedFailure == null ? null : expectedFailure.getClass(), cigar);
        if (expectedFailure != null) {
            if (expectedFailure instanceof UserException.MalformedRead) {
                // the same rejection, for the same reason (an implicit IndexOutOfBoundsException's message
                // is not compared: the JIT drops it once the throwing code is hot)
                Assert.assertEquals(actualFailure.getMessage(), expectedFailure.getMessage(), cigar);
            } else {
                // only flow mode reads a base past a deletion, which can run off the end of the read
                Assert.assertTrue(flowBased && expectedFailure instanceof IndexOutOfBoundsException,
                        "unexpected failure for " + cigar + ": " + expectedFailure);
            }
            // A rejected read may have added part of its evidence; start both again from the same state.
            expected.copyTo(actual);
        }
        return expectedFailure;
    }

    @DataProvider(name = "zeroLengthAndPadCigars")
    public Object[][] zeroLengthAndPadCigars() {
        // the cigar, and whether a read with it that lies inside the span is rejected as malformed (outside flow mode)
        return new Object[][]{
                {"0D10M", false},
                {"10M0D", true},
                {"10M2D0M", false},
                {"3M0I2D5M", false},
                {"2P1D5M", true},
                {"0M2D10M", true},
                {"4S0D6M", false},
                {"5H2D5M", true},
                {"3M2P2D5M", false},
                {"5M0N5M", false},
        };
    }

    /**
     * Zero-length elements and pads add no positions, so a zero-length deletion neither starts nor ends a read with a
     * deletion, while a deletion after only pads or zero-length elements still starts one; both walks agree on each.
     */
    @Test(dataProvider = "zeroLengthAndPadCigars")
    public void zeroLengthAndPadElementsAreSkippedAsThePileupStateMachineSkipsThem(final String cigar, final boolean malformedInside) {
        final int ploidy = 2;
        final double[][][] increments = ReferenceConfidenceModel.likelihoodIncrements(ploidy);
        final Random random = new Random(cigar.hashCode());
        final byte[] ref = randomBases(random, PADDED_SIZE);
        final int globalRefOffset = 10;
        final int spanStart = PADDED_START + globalRefOffset;
        final int spanSize = 200;
        final int readLength = TextCigarCodec.decode(cigar).getReadLength();
        for (final boolean[] mode : new boolean[][]{{true, false}, {false, false}, {true, true}, {false, true}}) {
            final boolean useSoftClippedBases = mode[0];
            final boolean flowBased = mode[1];
            final ReferenceConfidenceModel model = new ReferenceConfidenceModel(SampleList.singletonSampleList(SAMPLE), header,
                    MAX_INDEL_SIZE, -1, DELETION_QUAL, useSoftClippedBases, flowBased);
            // inside the span, straddling its start and straddling its end
            for (final int start : new int[]{spanStart + 5, spanStart - 3, spanStart + spanSize - 4}) {
                final byte[] quals = new byte[readLength];
                Arrays.fill(quals, (byte) 30);
                final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "read", 0, start, randomBases(random, readLength), quals, cigar);
                if (!useSoftClippedBases) {
                    read.setAttribute(ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_START_TAG, read.getSoftStart());
                    read.setAttribute(ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_END_TAG, read.getSoftEnd());
                }
                final ReferenceConfidenceModel.SiteEvidence actual = new ReferenceConfidenceModel.SiteEvidence(spanSize, ploidy + 1);
                final Accumulators expected = new Accumulators(spanSize, ploidy);
                final RuntimeException failure = sweepBoth(model, read, spanStart, spanSize, ref, globalRefOffset, increments,
                        actual, expected, useSoftClippedBases, flowBased);
                expected.assertSameAs(actual);
                if (start == spanStart + 5 && !flowBased) {
                    Assert.assertEquals(failure instanceof UserException.MalformedRead, malformedInside, cigar);
                    if (!malformedInside) {
                        Assert.assertNull(failure, cigar);
                        Assert.assertTrue(Arrays.stream(actual.readCounts).sum() > 0, cigar);
                    }
                }
            }
        }
    }

    @DataProvider(name = "longCigars")
    public Object[][] longCigars() {
        return new Object[][]{{false}, {true}};
    }

    /**
     * A BAM record with more cigar elements than a BAM cigar holds keeps them in its CG tag, and reports the
     * placeholder cigar's length as its number of elements; both walks must still walk the whole cigar, and reject
     * it alike when it ends with a deletion.
     */
    @Test(dataProvider = "longCigars")
    public void aReadWhoseCigarIsInItsCGTagIsSweptWhole(final boolean endsWithDeletion) throws IOException {
        final SAMFileHeader longHeader = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1_000_000);
        longHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord rg = new SAMReadGroupRecord("rg");
        rg.setSample(SAMPLE);
        longHeader.addReadGroup(rg);

        // 17,000 x 2M1I1M1D, then 2M unless the read is to end with the deletion: 68,001 or 68,000 elements
        final StringBuilder cigarString = new StringBuilder();
        for (int k = 0; k < 17_000; k++) {
            cigarString.append("2M1I1M1D");
        }
        if (!endsWithDeletion) {
            cigarString.append("2M");
        }
        final SAMRecord record = new SAMRecord(longHeader);
        record.setReadName("long");
        record.setReferenceIndex(0);
        record.setAlignmentStart(1_000);
        record.setCigar(TextCigarCodec.decode(cigarString.toString()));
        final Random random = new Random(11);
        final int readLength = record.getCigar().getReadLength();
        record.setReadBases(randomBases(random, readLength));
        final byte[] quals = new byte[readLength];
        for (int i = 0; i < readLength; i++) {
            quals[i] = (byte) random.nextInt(45);
        }
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        record.setAttribute("RG", "rg");

        final File bam = createTempFile("longCigar", ".bam");
        try (SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(longHeader, true, bam)) {
            writer.addAlignment(record);
        }
        final GATKRead read;
        try (SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bam)) {
            read = new SAMRecordToGATKReadAdapter(reader.iterator().next());
        }
        final List<CigarElement> elements = read.getCigarElements();
        Assert.assertEquals(elements.size(), endsWithDeletion ? 68_000 : 68_001);

        final ReferenceConfidenceModel model = new ReferenceConfidenceModel(SampleList.singletonSampleList(SAMPLE), longHeader,
                MAX_INDEL_SIZE, -1, DELETION_QUAL, true, false);
        final int ploidy = 2;
        final double[][][] increments = ReferenceConfidenceModel.likelihoodIncrements(ploidy);
        // a span inside the read, and one holding all of it, which the walks reach the end of
        for (final int[] span : new int[][]{{1_100, 400}, {900, read.getEnd() - 900 + 50}}) {
            final int globalRefOffset = 20;
            final byte[] ref = randomBases(random, span[1] + 2 * globalRefOffset);
            final ReferenceConfidenceModel.SiteEvidence actual = new ReferenceConfidenceModel.SiteEvidence(span[1], ploidy + 1);
            final Accumulators expected = new Accumulators(span[1], ploidy);
            final RuntimeException failure = sweepBoth(model, read, span[0], span[1], ref, globalRefOffset, increments,
                    actual, expected, true, false);
            expected.assertSameAs(actual);
            final boolean whole = span[0] < read.getStart();
            Assert.assertEquals(failure instanceof UserException.MalformedRead, whole && endsWithDeletion);
            if (failure == null) {
                // every position of the read inside the span is either a counted base or one dropped for low quality
                Assert.assertTrue(Arrays.stream(actual.readCounts).sum() > span[1] / 2);
            }
        }
    }

    private static RuntimeException sweep(final Runnable sweep) {
        try {
            sweep.run();
            return null;
        } catch (final RuntimeException e) {
            // malformed reads, and in flow mode a read ending with a deletion inside the span also runs off its bases
            return e;
        }
    }

    private static byte[] randomBases(final Random random, final int n) {
        final byte[] bases = new byte[n];
        for (int i = 0; i < n; i++) {
            bases[i] = random.nextInt(50) == 0 ? (byte) 'N' : BASES[random.nextInt(4)];
        }
        return bases;
    }

    /** A read overlapping or near the span, with a random cigar that sometimes starts or ends with a deletion. */
    private GATKRead randomRead(final Random random, final String name, final int spanStart, final int spanSize, final boolean useSoftClippedBases) {
        final StringBuilder cigar = new StringBuilder();
        int readLength = 0;
        if (random.nextInt(8) == 0) {
            cigar.append(1 + random.nextInt(5)).append('H');
        }
        if (random.nextInt(3) == 0) {
            final int n = 1 + random.nextInt(15);
            cigar.append(n).append('S');
            readLength += n;
        }
        final int nElements = 1 + random.nextInt(7);
        CigarOperator previous = null;
        for (int e = 0; e < nElements; e++) {
            // Now and then a zero-length element, a deletion first or last included, or a pad, which both walks skip.
            if (random.nextInt(12) == 0) {
                cigar.append(0).append(ZERO_LENGTH_OPERATORS[random.nextInt(ZERO_LENGTH_OPERATORS.length)].toString());
            }
            if (random.nextInt(16) == 0) {
                cigar.append(1 + random.nextInt(3)).append('P');
            }
            CigarOperator op;
            do {
                final int pick = random.nextInt(20);
                op = pick < 12 ? CigarOperator.M : pick < 15 ? CigarOperator.I : pick < 18 ? CigarOperator.D
                        : pick < 19 ? CigarOperator.N : (random.nextBoolean() ? CigarOperator.EQ : CigarOperator.X);
            } while (op == previous);
            // Rarely let a read start or end with a deletion, which both walks must reject alike.
            if (op == CigarOperator.D && (e == 0 || e == nElements - 1) && random.nextInt(10) != 0) {
                op = CigarOperator.M;
                if (op == previous) {
                    continue;
                }
            }
            final int n = op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X ? 1 + random.nextInt(60) : 1 + random.nextInt(6);
            cigar.append(n).append(op.toString());
            if (op.consumesReadBases()) {
                readLength += n;
            }
            previous = op;
        }
        if (random.nextInt(12) == 0) {
            cigar.append("0D");
        }
        if (random.nextInt(3) == 0) {
            final int n = 1 + random.nextInt(15);
            cigar.append(n).append('S');
            readLength += n;
        }
        if (random.nextInt(8) == 0) {
            cigar.append(1 + random.nextInt(5)).append('H');
        }
        if (readLength == 0) {
            cigar.append("5M");
            readLength = 5;
        }
        final byte[] bases = randomBases(random, readLength);
        final byte[] quals = new byte[readLength];
        for (int i = 0; i < readLength; i++) {
            quals[i] = (byte) (random.nextInt(5) == 0 ? random.nextInt(8) : random.nextInt(45));
        }
        final int start = spanStart - 150 + random.nextInt(spanSize + 200);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, name, 0, Math.max(1, start), bases, quals, cigar.toString());
        read.setIsReverseStrand(random.nextBoolean());
        if (random.nextBoolean()) {
            // a pair whose fragment may be short enough for the read to run into the adaptor
            final int fragment = 20 + random.nextInt(200);
            read.setMatePosition(read.getContig(), read.isReverseStrand() ? Math.max(1, read.getEnd() - fragment + 1) : read.getStart());
            read.setMateIsReverseStrand(!read.isReverseStrand());
            read.setFragmentLength(read.isReverseStrand() ? -fragment : fragment);
        }
        if (!useSoftClippedBases) {
            final int softStart = read.getSoftStart();
            final int softEnd = read.getSoftEnd();
            read.setAttribute(ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_START_TAG, softStart + random.nextInt(10));
            read.setAttribute(ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_END_TAG, softEnd - random.nextInt(10));
        }
        return read;
    }

    /** Per-position results and counts of one span. */
    private static final class Accumulators {
        final RefVsAnyResult[] sites;
        final int[] readCounts;
        final int[] indelInformativeReads;

        Accumulators(final int spanSize, final int ploidy) {
            sites = new RefVsAnyResult[spanSize];
            for (int i = 0; i < spanSize; i++) {
                sites[i] = new RefVsAnyResult(ploidy + 1);
            }
            readCounts = new int[spanSize];
            indelInformativeReads = new int[spanSize];
        }

        void copyTo(final ReferenceConfidenceModel.SiteEvidence evidence) {
            for (int i = 0; i < sites.length; i++) {
                System.arraycopy(sites[i].genotypeLikelihoods, 0, evidence.likelihoods, i * evidence.likelihoodCount, evidence.likelihoodCount);
                evidence.refDepth[i] = sites[i].refDepth;
                evidence.nonRefDepth[i] = sites[i].nonRefDepth;
            }
            System.arraycopy(readCounts, 0, evidence.readCounts, 0, readCounts.length);
            System.arraycopy(indelInformativeReads, 0, evidence.indelInformativeReads, 0, indelInformativeReads.length);
        }

        void assertSameAs(final ReferenceConfidenceModel.SiteEvidence actual) {
            Assert.assertEquals(actual.readCounts, readCounts);
            Assert.assertEquals(actual.indelInformativeReads, indelInformativeReads);
            for (int i = 0; i < sites.length; i++) {
                final RefVsAnyResult actualSite = actual.toResult(i);
                Assert.assertEquals(actualSite.refDepth, sites[i].refDepth, "refDepth at " + i);
                Assert.assertEquals(actualSite.nonRefDepth, sites[i].nonRefDepth, "nonRefDepth at " + i);
                for (int k = 0; k < sites[i].genotypeLikelihoods.length; k++) {
                    // the same additions in the same order give the same bits
                    Assert.assertEquals(Double.doubleToRawLongBits(actualSite.genotypeLikelihoods[k]),
                            Double.doubleToRawLongBits(sites[i].genotypeLikelihoods[k]), "likelihood " + k + " at " + i);
                }
            }
        }
    }

    /** The position-at-a-time walk of a read with the pileup's state machine that {@code sweepRead} replaces. */
    private static void referenceSweep(final ReferenceConfidenceModel model, final GATKRead read, final int spanStart, final int spanSize,
                                       final byte[] ref, final int globalRefOffset, final double[][][] increments, final Accumulators acc,
                                       final boolean useSoftClippedBases, final boolean flowBasedModel) {
        final boolean skipOriginalSoftClips = !useSoftClippedBases;
        final int originalSoftStart = skipOriginalSoftClips ? read.getAttributeAsInteger(ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_START_TAG) : 0;
        final int originalSoftEnd = skipOriginalSoftClips ? read.getAttributeAsInteger(ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_END_TAG) : 0;

        final List<CigarElement> elements = read.getCigarElements();
        final AlignmentStateMachine state = new AlignmentStateMachine(read);
        int cigarElementIndex = -1;
        int alignedOffsetOfElement = 0;
        BitSet indelInformativeBases = null;

        while (state.stepForwardOnGenome() != null) {
            final int position = state.getGenomePosition();
            if (position < spanStart) {
                continue;
            }
            final int i = position - spanStart;
            if (i >= spanSize) {
                break;
            }
            final CigarElement element = state.getCurrentCigarElement();
            final CigarOperator op = element.getOperator();
            if (op == CigarOperator.N) {
                continue;
            }
            if (ReadUtils.isBaseInsideAdaptor(read, position)) {
                continue;
            }
            final int elementIndex = state.getCurrentCigarElementOffset();
            while (cigarElementIndex < elementIndex) {
                if (cigarElementIndex >= 0) {
                    alignedOffsetOfElement += ReferenceConfidenceModel.alignedLength(elements.get(cigarElementIndex));
                }
                cigarElementIndex++;
            }
            final int readOffset = state.getReadOffset();
            final boolean isDeletion = op == CigarOperator.D;
            final byte refBase = ref[i + globalRefOffset];

            final byte qual = isDeletion ? model.deletionQuality(read, readOffset, refBase, true) : read.getBaseQuality(readOffset);
            final boolean usable = !((qual <= ReferenceConfidenceModel.BASE_QUAL_THRESHOLD) && (flowBasedModel || !isDeletion))
                    && !(skipOriginalSoftClips && (originalSoftStart > position || originalSoftEnd < position));
            if (usable) {
                acc.readCounts[i]++;
                final boolean isAlt = isDeletion || read.getBase(readOffset) != refBase;
                // as the pileup path adds an element, at read weight 1
                ReferenceConfidenceModel.applyRefVsNonRefLikelihoodAndCount(increments, acc.sites[i], isAlt, qual, 1.0);
            }

            if (acc.indelInformativeReads[i] >= ReferenceConfidenceModel.MAX_N_INDEL_INFORMATIVE_READS) {
                continue;
            }
            final int offsetInElement = state.getOffsetIntoCurrentCigarElement();
            final boolean beforeIndel = offsetInElement == element.getLength() - 1
                    && (ReferenceConfidenceModel.nextOnGenomeOperatorIsDeletion(elements, elementIndex)
                        || ReferenceConfidenceModel.nextOperatorIs(elements, elementIndex, CigarOperator.I));
            if (!isDeletion && !beforeIndel) {
                final int alignedOffset = alignedOffsetOfElement + (ReferenceConfidenceModel.alignedLength(element) > 0 ? offsetInElement : 0);
                if (indelInformativeBases == null) {
                    indelInformativeBases = ReferenceConfidenceModel.indelInformativeBases(read, alignedOffset, ref, i + globalRefOffset, MAX_INDEL_SIZE);
                }
                if (indelInformativeBases.get(alignedOffset)) {
                    acc.indelInformativeReads[i]++;
                }
            }
        }
    }
}
