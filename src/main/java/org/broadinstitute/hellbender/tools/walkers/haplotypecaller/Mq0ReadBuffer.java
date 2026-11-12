package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.hash.Hashing;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;

/**
 * The MAPQ-0 reads HaplotypeCaller's rescue (see {@link Mq0RescueArgumentCollection}) sets aside as the read stream
 * is read, until the regions they overlap have been called. At most a given number of reads is kept per alignment
 * start: those whose names hash lowest, so that neither the choice nor the order in which they are returned depends on
 * the order in which reads with the same start arrive (a coordinate sort puts one strand first), and no random
 * generator is used.
 *
 * <p>Reads are forgotten from the front once they end before the region asked for. A long read holds back the
 * shorter ones after it until it ends, which only delays their release.</p>
 */
final class Mq0ReadBuffer {

    private final SAMSequenceDictionary dictionary;
    private final int maxReadsPerStart;
    // reads whose alignment start is complete, in coordinate order
    private final ArrayDeque<GATKRead> reads = new ArrayDeque<>();
    // the reads kept at the latest alignment start, the one to drop next first
    private final PriorityQueue<Keyed> latestStart = new PriorityQueue<>(Comparator.reverseOrder());
    private int peakSize = 0;
    // the contig index and start of the last read added, to check the coordinate order
    private int lastContig = -1;
    private int lastStart = 0;

    /** A read with the key that orders reads with the same start: a hash of the name, then the name. */
    private record Keyed(GATKRead read, int hash) implements Comparable<Keyed> {
        Keyed(final GATKRead read) {
            this(read, Hashing.murmur3_32_fixed().hashUnencodedChars(read.getName()).asInt());
        }

        @Override
        public int compareTo(final Keyed other) {
            final int c = Integer.compare(hash, other.hash);
            return c != 0 ? c : read.getName().compareTo(other.read.getName());
        }
    }

    /**
     * @param maxReadsPerStart the most reads to keep per alignment start, or 0 or less for no limit
     */
    Mq0ReadBuffer(final SAMSequenceDictionary dictionary, final int maxReadsPerStart) {
        this.dictionary = Utils.nonNull(dictionary);
        this.maxReadsPerStart = maxReadsPerStart;
    }

    /**
     * Adds a read. Reads must be added in coordinate order, as they come from the read stream.
     */
    void add(final GATKRead read) {
        final int contig = contigIndex(read.getContig());
        Utils.validate(contig > lastContig || (contig == lastContig && read.getStart() >= lastStart),
                () -> "MAPQ-0 read " + read + " is out of coordinate order");
        if (contig != lastContig || read.getStart() != lastStart) {
            // in key order, so that the order of the reads returned does not depend on their arrival order either
            latestStart.stream().sorted().forEach(keyed -> reads.add(keyed.read()));
            latestStart.clear();
        }
        lastContig = contig;
        lastStart = read.getStart();
        latestStart.add(new Keyed(read));
        if (maxReadsPerStart > 0 && latestStart.size() > maxReadsPerStart) {
            latestStart.poll();
        }
        peakSize = Math.max(peakSize, size());
    }

    /**
     * Returns the reads overlapping a region's padded span, first forgetting those that end before it. Regions must be
     * asked for in coordinate order, so that forgotten reads are never needed again.
     */
    List<GATKRead> overlapping(final SimpleInterval paddedSpan) {
        final int spanContig = contigIndex(paddedSpan.getContig());
        while (!reads.isEmpty() && endsBefore(reads.peekFirst(), spanContig, paddedSpan.getStart())) {
            reads.pollFirst();
        }
        final List<GATKRead> result = new ArrayList<>();
        for (final GATKRead r : reads) {
            if (!r.getContig().equals(paddedSpan.getContig()) || r.getStart() > paddedSpan.getEnd()) {
                return result;
            }
            if (r.overlaps(paddedSpan)) {
                result.add(r);
            }
        }
        latestStart.stream().sorted().map(Keyed::read).filter(r -> r.overlaps(paddedSpan)).forEach(result::add);
        return result;
    }

    private boolean endsBefore(final GATKRead read, final int spanContig, final int spanStart) {
        final int contig = contigIndex(read.getContig());
        return contig < spanContig || (contig == spanContig && read.getEnd() < spanStart);
    }

    private int contigIndex(final String contig) {
        final int index = dictionary.getSequenceIndex(contig);
        Utils.validate(index >= 0, () -> "contig " + contig + " is not in the sequence dictionary");
        return index;
    }

    /** @return the number of reads held */
    int size() {
        return reads.size() + latestStart.size();
    }

    /** @return the largest number of reads held at once */
    int peakSize() {
        return peakSize;
    }
}
