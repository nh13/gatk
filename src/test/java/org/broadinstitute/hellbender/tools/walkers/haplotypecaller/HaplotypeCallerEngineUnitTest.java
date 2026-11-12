package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.iterators.ReadFilteringIterator;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class HaplotypeCallerEngineUnitTest extends GATKBaseTest {

    @Test
    public void testIsActive() throws IOException {
        final File testBam = new File(NA12878_20_21_WGS_bam);
        final Path reference = Paths.get(b37_reference_20_21);
        final SimpleInterval shardInterval = new SimpleInterval("20", 10000000, 10001000);
        final SimpleInterval paddedShardInterval = new SimpleInterval(shardInterval.getContig(), shardInterval.getStart() - 100, shardInterval.getEnd() + 100);
        final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();

        // We expect isActive() to return 1.0 for the sites below, and 0.0 for all other sites
        final List<SimpleInterval> expectedActiveSites = Arrays.asList(
                new SimpleInterval("20", 9999996, 9999996),
                new SimpleInterval("20", 9999997, 9999997),
                new SimpleInterval("20", 10000117, 10000117),
                new SimpleInterval("20", 10000211, 10000211),
                new SimpleInterval("20", 10000439, 10000439),
                new SimpleInterval("20", 10000598, 10000598),
                new SimpleInterval("20", 10000694, 10000694),
                new SimpleInterval("20", 10000758, 10000758),
                new SimpleInterval("20", 10001019, 10001019)
        );

        try (final ReadsDataSource reads = new ReadsPathDataSource(testBam.toPath());
             final ReferenceDataSource ref = new ReferenceFileSource(reference);
             final CachingIndexedFastaSequenceFile referenceReader = new CachingIndexedFastaSequenceFile(reference)) {

            final HaplotypeCallerEngine hcEngine = new HaplotypeCallerEngine(hcArgs, new AssemblyRegionArgumentCollection(), false, false, reads.getHeader(), referenceReader, new VariantAnnotatorEngine(new ArrayList<>(), hcArgs.dbsnp.dbsnp, hcArgs.comps, false, false));

            List<ReadFilter> hcFilters = HaplotypeCallerEngine.makeStandardHCReadFilters();
            hcFilters.forEach(filter -> filter.setHeader(reads.getHeader()));
            ReadFilter hcCombinedFilter = hcFilters.get(0);
            for ( int i = 1; i < hcFilters.size(); ++i ) {
                hcCombinedFilter = hcCombinedFilter.and(hcFilters.get(i));
            }
            final Iterator<GATKRead> readIter = new ReadFilteringIterator(reads.query(paddedShardInterval), hcCombinedFilter);

            final LocusIteratorByState libs = new LocusIteratorByState(readIter, DownsamplingMethod.NONE, ReadUtils.getSamplesFromHeader(reads.getHeader()), reads.getHeader(), true);

            libs.forEachRemaining(pileup -> {
                final SimpleInterval pileupInterval = new SimpleInterval(pileup.getLocation());
                final ReferenceContext pileupRefContext = new ReferenceContext(ref, pileupInterval);

                final ActivityProfileState isActiveResult = hcEngine.isActive(pileup, pileupRefContext, new FeatureContext((FeatureManager)null, pileupInterval));

                final double expectedIsActiveValue = expectedActiveSites.contains(pileupInterval) ? 1.0 : 0.0;
                Assert.assertEquals(isActiveResult.isActiveProb(), expectedIsActiveValue, "Wrong isActive probability for site " + pileupInterval);
            });
        }
    }

    private static HaplotypeCallerEngine engineFor(final SAMFileHeader header, final CachingIndexedFastaSequenceFile referenceReader) {
        final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        return new HaplotypeCallerEngine(hcArgs, new AssemblyRegionArgumentCollection(), false, false, header, referenceReader,
                new VariantAnnotatorEngine(new ArrayList<>(), hcArgs.dbsnp.dbsnp, hcArgs.comps, false, false));
    }

    /** A header over the given dictionary with one read group per sample, each read group named after its sample. */
    private static SAMFileHeader headerWithSamples(final SAMSequenceDictionary dictionary, final String... samples) {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(dictionary);
        for (final String sample : samples) {
            final SAMReadGroupRecord readGroup = new SAMReadGroupRecord(sample);
            readGroup.setSample(sample);
            header.addReadGroup(readGroup);
        }
        return header;
    }

    private static GATKRead readAt(final SAMFileHeader header, final String name, final int start, final int mappingQuality, final String sample) {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, name, 0, start, 50);
        read.setMappingQuality(mappingQuality);
        read.setReadGroup(sample);
        return read;
    }

    @Test
    public void filterNonPassingReadsRemovesEveryCopyOfADuplicatedRecord() throws IOException {
        try (final CachingIndexedFastaSequenceFile referenceReader = new CachingIndexedFastaSequenceFile(Paths.get(b37_reference_20_21))) {
            final SAMFileHeader header = headerWithSamples(referenceReader.getSequenceDictionary(), "sample1");
            final HaplotypeCallerEngine engine = engineFor(header, referenceReader);
            final AssemblyRegion region = new AssemblyRegion(new SimpleInterval("20", 10000000, 10001000), 0, header);
            final GATKRead lowMappingQuality = readAt(header, "dup", 10000100, 5, "sample1");
            final GATKRead equalCopy = lowMappingQuality.copy();
            final GATKRead passing = readAt(header, "ok", 10000100, 60, "sample1");
            region.addAll(Arrays.asList(lowMappingQuality, equalCopy, passing));
            Assert.assertEquals(equalCopy, lowMappingQuality);

            final List<GATKRead> filtered = engine.filterNonPassingReads(region);

            Assert.assertEquals(region.getReads(), Collections.singletonList(passing));
            Assert.assertEquals(filtered, Arrays.asList(lowMappingQuality, equalCopy));
        }
    }

    @Test
    public void removeReadsFromAllSamplesExceptRemovesEveryCopyOfADuplicatedRecord() throws IOException {
        try (final CachingIndexedFastaSequenceFile referenceReader = new CachingIndexedFastaSequenceFile(Paths.get(b37_reference_20_21))) {
            final SAMFileHeader header = headerWithSamples(referenceReader.getSequenceDictionary(), "sample1", "sample2");
            final HaplotypeCallerEngine engine = engineFor(header, referenceReader);
            final AssemblyRegion region = new AssemblyRegion(new SimpleInterval("20", 10000000, 10001000), 0, header);
            final GATKRead otherSample = readAt(header, "other", 10000100, 60, "sample2");
            final GATKRead equalCopy = otherSample.copy();
            final GATKRead kept = readAt(header, "keep", 10000100, 60, "sample1");
            region.addAll(Arrays.asList(otherSample, equalCopy, kept));

            engine.removeReadsFromAllSamplesExcept("sample1", region);

            Assert.assertEquals(region.getReads(), Collections.singletonList(kept));
        }
    }

    @Test
    public void testPassesGenotypingFilters() {
        // the filters genotyping applies, and that MAPQ-0 reads must pass to count for the rescue
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(2, 1, 1000);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "read", 0, 100, 50);
        read.setMappingQuality(0);
        read.setIsPaired(true);
        read.setMatePosition(header.getSequence(0).getSequenceName(), 300);
        Assert.assertTrue(HaplotypeCallerEngine.passesGenotypingFilters(read, 0, null));
        Assert.assertFalse(HaplotypeCallerEngine.passesGenotypingFilters(read, 20, null));
        read.setMatePosition(header.getSequence(1).getSequenceName(), 300);
        Assert.assertFalse(HaplotypeCallerEngine.passesGenotypingFilters(read, 0, null));
        final GATKRead shortRead = ArtificialReadUtils.createArtificialRead(header, "short", 0, 100, 5);
        Assert.assertFalse(HaplotypeCallerEngine.passesGenotypingFilters(shortRead, 0, null));
        // --keep-read-group: a read in another read group, or in none, is left out rather than failing
        final GATKRead noReadGroup = ArtificialReadUtils.createArtificialRead(header, "noReadGroup", 0, 100, 50);
        noReadGroup.setReadGroup(null);
        Assert.assertFalse(HaplotypeCallerEngine.passesGenotypingFilters(noReadGroup, 0, "rg"));
        noReadGroup.setReadGroup("other");
        Assert.assertFalse(HaplotypeCallerEngine.passesGenotypingFilters(noReadGroup, 0, "rg"));
        noReadGroup.setReadGroup("rg");
        Assert.assertTrue(HaplotypeCallerEngine.passesGenotypingFilters(noReadGroup, 0, "rg"));
    }
}
