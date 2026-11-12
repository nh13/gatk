package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.hash.Hashing;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Unit tests for {@link Mq0ReadBuffer}, which holds HaplotypeCaller's MAPQ-0 rescue reads.
 */
public class Mq0ReadBufferUnitTest extends GATKBaseTest {

    private static final SAMFileHeader HEADER = ArtificialReadUtils.createArtificialSamHeader(2, 1, 100_000);
    private static final String CONTIG_1 = HEADER.getSequence(0).getSequenceName();
    private static final String CONTIG_2 = HEADER.getSequence(1).getSequenceName();

    private static GATKRead read(final String name, final int contig, final int start, final int length) {
        return ArtificialReadUtils.createArtificialRead(HEADER, name, contig, start, length);
    }

    private static Set<String> names(final List<GATKRead> reads) {
        return reads.stream().map(GATKRead::getName).collect(Collectors.toSet());
    }

    @Test
    public void testPerStartCap() {
        final List<GATKRead> sameStart = new ArrayList<>();
        for (int i = 0; i < 20; i++) {
            sameStart.add(read("read" + i, 0, 100, 50));
        }
        final Mq0ReadBuffer buffer = new Mq0ReadBuffer(HEADER.getSequenceDictionary(), 5);
        sameStart.forEach(buffer::add);
        // the same start on the other contig is a different start
        buffer.add(read("otherContig", 1, 100, 50));
        Assert.assertEquals(buffer.size(), 6);
        final List<String> kept = buffer.overlapping(new SimpleInterval(CONTIG_1, 1, 1000)).stream().map(GATKRead::getName).collect(Collectors.toList());
        // the reads whose names hash lowest (then by name), in that order
        final Comparator<String> byHash = Comparator.<String>comparingInt(name -> Hashing.murmur3_32_fixed().hashUnencodedChars(name).asInt())
                .thenComparing(Comparator.naturalOrder());
        Assert.assertEquals(kept, sameStart.stream().map(GATKRead::getName).sorted(byHash).limit(5).collect(Collectors.toList()));

        // neither which reads are kept nor their order depends on the order they arrive in (here, on the strand order
        // of a sorter), whether they are still at the latest start or not
        final List<GATKRead> reversed = new ArrayList<>(sameStart);
        Collections.reverse(reversed);
        final Mq0ReadBuffer reverseBuffer = new Mq0ReadBuffer(HEADER.getSequenceDictionary(), 5);
        reversed.forEach(reverseBuffer::add);
        Assert.assertEquals(reverseBuffer.overlapping(new SimpleInterval(CONTIG_1, 1, 1000)).stream().map(GATKRead::getName).collect(Collectors.toList()), kept);
        reverseBuffer.add(read("later", 0, 101, 50));
        Assert.assertEquals(reverseBuffer.overlapping(new SimpleInterval(CONTIG_1, 1, 100)).stream().map(GATKRead::getName).collect(Collectors.toList()), kept);

        // no cap
        final Mq0ReadBuffer uncapped = new Mq0ReadBuffer(HEADER.getSequenceDictionary(), 0);
        sameStart.forEach(uncapped::add);
        Assert.assertEquals(uncapped.overlapping(new SimpleInterval(CONTIG_1, 1, 1000)).size(), 20);
        Assert.assertEquals(uncapped.peakSize(), 20);
    }

    @Test
    public void testOverlapping() {
        final Mq0ReadBuffer buffer = new Mq0ReadBuffer(HEADER.getSequenceDictionary(), 0);
        buffer.add(read("a", 0, 100, 50));   // 100-149
        buffer.add(read("b", 0, 110, 200));  // 110-309
        buffer.add(read("c", 0, 150, 50));   // 150-199
        buffer.add(read("d", 0, 300, 50));   // 300-349
        buffer.add(read("e", 0, 301, 50));   // 301-350

        Assert.assertEquals(names(buffer.overlapping(new SimpleInterval(CONTIG_1, 149, 150))), Set.of("a", "b", "c"));
        // a read ending just before the span is forgotten, one starting just after it is not returned
        Assert.assertEquals(names(buffer.overlapping(new SimpleInterval(CONTIG_1, 200, 300))), Set.of("b", "d"));
        Assert.assertEquals(buffer.size(), 4);
        // the reads at the latest start are returned too, before any read with a later start arrives
        Assert.assertEquals(names(buffer.overlapping(new SimpleInterval(CONTIG_1, 340, 1000))), Set.of("d", "e"));
        Assert.assertEquals(buffer.size(), 2);
        // moving to the next contig forgets the reads of the first
        buffer.add(read("f", 1, 50, 50));
        Assert.assertEquals(names(buffer.overlapping(new SimpleInterval(CONTIG_2, 1, 60))), Set.of("f"));
        Assert.assertEquals(buffer.size(), 1);
        Assert.assertEquals(buffer.peakSize(), 5);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testOutOfOrder() {
        final Mq0ReadBuffer buffer = new Mq0ReadBuffer(HEADER.getSequenceDictionary(), 0);
        buffer.add(read("a", 0, 100, 50));
        buffer.add(read("b", 0, 99, 50));
    }
}
