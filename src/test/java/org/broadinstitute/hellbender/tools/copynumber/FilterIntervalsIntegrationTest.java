package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationMap;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.CopyNumberAnnotations;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Integration tests for {@link FilterIntervals}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class FilterIntervalsIntegrationTest extends CommandLineProgramTest {
    private static final File REFERENCE_FILE = new File(b37_reference_20_21);
    private static final SAMSequenceDictionary SEQUENCE_DICTIONARY = ReferenceDataSource.of(REFERENCE_FILE.toPath()).getSequenceDictionary();
    private static final LocatableMetadata LOCATABLE_METADATA = new SimpleLocatableMetadata(SEQUENCE_DICTIONARY);

    private static final AnnotatedIntervalCollection ANNOTATED_INTERVALS = new AnnotatedIntervalCollection(
            LOCATABLE_METADATA,
            Arrays.asList(
                    new AnnotatedInterval(new SimpleInterval("20", 1,	1),
                            new AnnotationMap(Arrays.asList(
                                    Pair.of(CopyNumberAnnotations.GC_CONTENT, 0.5),
                                    Pair.of(CopyNumberAnnotations.MAPPABILITY, 0.05),
                                    Pair.of(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT, 0.05)))),
                    new AnnotatedInterval(new SimpleInterval("20", 2,	2),
                            new AnnotationMap(Arrays.asList(
                                    Pair.of(CopyNumberAnnotations.GC_CONTENT, 0.5),
                                    Pair.of(CopyNumberAnnotations.MAPPABILITY, 0.5),
                                    Pair.of(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT, 0.95)))),
                    new AnnotatedInterval(new SimpleInterval("20", 3,	3),
                            new AnnotationMap(Arrays.asList(
                                    Pair.of(CopyNumberAnnotations.GC_CONTENT, 0.5),
                                    Pair.of(CopyNumberAnnotations.MAPPABILITY, 0.5),
                                    Pair.of(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT, 0.5)))),
                    new AnnotatedInterval(new SimpleInterval("20", 4,	4),
                            new AnnotationMap(Arrays.asList(
                                    Pair.of(CopyNumberAnnotations.GC_CONTENT, 0.05),
                                    Pair.of(CopyNumberAnnotations.MAPPABILITY, 0.5),
                                    Pair.of(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT, 0.5)))),
                    new AnnotatedInterval(new SimpleInterval("20", 5,	5),
                            new AnnotationMap(Arrays.asList(
                                    Pair.of(CopyNumberAnnotations.GC_CONTENT, 0.95),
                                    Pair.of(CopyNumberAnnotations.MAPPABILITY, 0.95),
                                    Pair.of(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT, 0.5))))));

    @DataProvider(name = "dataAnnotationBasedFilters")
    public Object[][] dataAnnotationBasedFilters() {
        final File annotatedIntervalsFile = createTempFile("annotated-intervals", ".tsv");
        ANNOTATED_INTERVALS.write(annotatedIntervalsFile);
        final File intervalsFile = createTempFile("intervals", ".interval_list");
        final IntervalList intervals = new IntervalList(ANNOTATED_INTERVALS.getMetadata().getSequenceDictionary());
        ANNOTATED_INTERVALS.getIntervals().forEach(i -> intervals.add(new Interval(i)));
        intervals.write(intervalsFile);

        return new Object[][]{
                //intervals file, annotated-intervals file,
                //min/max GC content, mix/max mappability, min/max seg-dupe content, expected array of indices of retained intervals
                {intervalsFile, annotatedIntervalsFile, 0., 1., 0., 1., 0., 1., Arrays.asList(0, 1, 2, 3, 4)},
                {intervalsFile, annotatedIntervalsFile, 0.1, 0.9, 0., 1., 0., 1., Arrays.asList(0, 1, 2)},
                {intervalsFile, annotatedIntervalsFile, 0., 1., 0.1, 0.9, 0., 1., Arrays.asList(1, 2, 3)},
                {intervalsFile, annotatedIntervalsFile, 0., 1., 0., 1., 0.1, 0.9, Arrays.asList(2, 3, 4)},
                {intervalsFile, annotatedIntervalsFile, 0.1, 0.9, 0.1, 0.9, 0., 1., Arrays.asList(1, 2)},
                {intervalsFile, annotatedIntervalsFile, 0.1, 0.9, 0., 1., 0.1, 0.9, Collections.singletonList(2)},
                {intervalsFile, annotatedIntervalsFile, 0., 1., 0.1, 0.9, 0.1, 0.9, Arrays.asList(2, 3)},
                {intervalsFile, annotatedIntervalsFile, 0.1, 0.9, 0.1, 0.9, 0.1, 0.9, Collections.singletonList(2)}};
    }

    @Test(dataProvider = "dataAnnotationBasedFilters")
    public void testAnnotationBasedFilters(final File intervalsFile,
                                           final File annotatedIntervalsFile,
                                           final double minimumGCContent,
                                           final double maximumGCContent,
                                           final double minimumMappability,
                                           final double maximumMappability,
                                           final double minimumSegmentalDuplicationContent,
                                           final double maximumSegmentalDuplicationContent,
                                           final List<Integer> expectedIndices) {
        final File outputFile = createTempFile("filter-intervals-test", ".interval_list");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, intervalsFile.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_LONG_NAME, annotatedIntervalsFile.getAbsolutePath())
                .addArgument(FilterIntervals.MINIMUM_GC_CONTENT_LONG_NAME, Double.toString(minimumGCContent))
                .addArgument(FilterIntervals.MAXIMUM_GC_CONTENT_LONG_NAME, Double.toString(maximumGCContent))
                .addArgument(FilterIntervals.MINIMUM_MAPPABILITY_LONG_NAME, Double.toString(minimumMappability))
                .addArgument(FilterIntervals.MAXIMUM_MAPPABILITY_LONG_NAME, Double.toString(maximumMappability))
                .addArgument(FilterIntervals.MINIMUM_SEGMENTAL_DUPLICATION_CONTENT_LONG_NAME, Double.toString(minimumSegmentalDuplicationContent))
                .addArgument(FilterIntervals.MAXIMUM_SEGMENTAL_DUPLICATION_CONTENT_LONG_NAME, Double.toString(maximumSegmentalDuplicationContent))
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        final IntervalList result = IntervalList.fromFile(outputFile);
        final IntervalList all = IntervalList.fromFile(intervalsFile);
        final List<Interval> allIntervals = all.getIntervals();
        final IntervalList expected = new IntervalList(all.getHeader().getSequenceDictionary());
        expectedIndices.stream().map(allIntervals::get).map(Interval::new).forEach(expected::add);
        Assert.assertEquals(result, expected);
        Assert.assertNotSame(result, expected);
    }

    @DataProvider(name = "dataCountBasedFilters")
    public Object[][] dataCountBasedFilters() {
        final int numSamples = 100;
        final int numIntervals = 110;
        final int numIntervalsBelowCountThreshold = 10;
        final int numUncorruptedSamples = 10;
        final int lowCountFilterCountThreshold = 5;
        final double epsilon = 1E-1;
        final double percentageOfSamples = 100. * (numSamples - numUncorruptedSamples) / numSamples - epsilon;

        final File intervalsFile = createTempFile("intervals", ".interval_list");
        final List<File> countFiles = new ArrayList<>(numSamples);
        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
            final String sampleName = String.format("sample-%d", sampleIndex);
            final SampleLocatableMetadata sampleMetadata = new SimpleSampleLocatableMetadata(sampleName, SEQUENCE_DICTIONARY);
            final List<SimpleCount> counts = new ArrayList<>(numIntervals);
            for (int intervalIndex = 0; intervalIndex < numIntervals; intervalIndex++) {
                final SimpleInterval interval = new SimpleInterval("20", intervalIndex + 1, intervalIndex + 1);
                if (intervalIndex < numIntervalsBelowCountThreshold && sampleIndex >= numUncorruptedSamples) {
                    //corrupt first numIntervalsBelowCountThreshold intervals in last (numSamples - numUncorruptedSamples) samples with low counts
                    counts.add(new SimpleCount(interval, lowCountFilterCountThreshold - 1));
                } else {
                    counts.add(new SimpleCount(interval, intervalIndex + lowCountFilterCountThreshold));
                }
            }
            final SimpleCountCollection sampleCounts = new SimpleCountCollection(sampleMetadata, counts);
            final File countFile = createTempFile(sampleName, ".tsv");
            sampleCounts.write(countFile);
            countFiles.add(countFile);

            if (sampleIndex == 0) {
                final IntervalList intervals = new IntervalList(sampleCounts.getMetadata().getSequenceDictionary());
                sampleCounts.getIntervals().forEach(i -> intervals.add(new Interval(i)));
                intervals.write(intervalsFile);
            }
        }

        return new Object[][]{
                //intervals file, count files, lowCountFilterCountThreshold, lowCountFilterPercentageOfSamples,
                //extremeCountFilterMinimumPercentile, extremeCountFilterMaximumPercentile, extremeCountFilterPercentageOfSamples,
                //expected array of indices of retained intervals
                {intervalsFile, countFiles, lowCountFilterCountThreshold, percentageOfSamples, 1., 99., percentageOfSamples,
                        IntStream.range(numIntervalsBelowCountThreshold + 1, numIntervalsBelowCountThreshold + 99).boxed().collect(Collectors.toList())},
                {intervalsFile, countFiles, lowCountFilterCountThreshold, percentageOfSamples, 1., 90., percentageOfSamples,
                        IntStream.range(numIntervalsBelowCountThreshold + 1, numIntervalsBelowCountThreshold + 90).boxed().collect(Collectors.toList())}};
    }

    @Test(dataProvider = "dataCountBasedFilters")
    public void testCountBasedFilters(final File intervalsFile,
                                      final List<File> countFiles,
                                      final int lowCountFilterCountThreshold,
                                      final double lowCountFilterPercentageOfSamples,
                                      final double extremeCountFilterMinimumPercentile,
                                      final double extremeCountFilterMaximumPercentile,
                                      final double extremeCountFilterPercentageOfSamples,
                                      final List<Integer> expectedIndices) {
        final File outputFile = createTempFile("filter-intervals-test", ".interval_list");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, intervalsFile.getAbsolutePath())
                .addArgument(FilterIntervals.LOW_COUNT_FILTER_COUNT_THRESHOLD_LONG_NAME, Integer.toString(lowCountFilterCountThreshold))
                .addArgument(FilterIntervals.LOW_COUNT_FILTER_PERCENTAGE_OF_SAMPLES_LONG_NAME, Double.toString(lowCountFilterPercentageOfSamples))
                .addArgument(FilterIntervals.EXTREME_COUNT_FILTER_MINIMUM_PERCENTILE_LONG_NAME, Double.toString(extremeCountFilterMinimumPercentile))
                .addArgument(FilterIntervals.EXTREME_COUNT_FILTER_MAXIMUM_PERCENTILE_LONG_NAME, Double.toString(extremeCountFilterMaximumPercentile))
                .addArgument(FilterIntervals.EXTREME_COUNT_FILTER_PERCENTAGE_OF_SAMPLES_LONG_NAME, Double.toString(extremeCountFilterPercentageOfSamples))
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addOutput(outputFile);
        countFiles.forEach(argsBuilder::addInput);
        runCommandLine(argsBuilder);
        final IntervalList result = IntervalList.fromFile(outputFile);
        final IntervalList all = IntervalList.fromFile(intervalsFile);
        final List<Interval> allIntervals = all.getIntervals();
        final IntervalList expected = new IntervalList(all.getHeader().getSequenceDictionary());
        expectedIndices.stream().map(allIntervals::get).map(Interval::new).forEach(expected::add);
        Assert.assertEquals(result, expected);
        Assert.assertNotSame(result, expected);
    }

    @DataProvider(name = "dataAllFilters")
    public Object[][] dataAllFilters() {
        final int numSamples = 100;
        final int numIntervals = 120;
        final int numIntervalsBelowCountThreshold = 10;
        final int numIntervalsFailingAnnotationFilters = 10;
        final int numUncorruptedSamples = 10;
        final int lowCountFilterCountThreshold = 5;
        final double epsilon = 1E-1;
        final double percentageOfSamples = 100. * (numSamples - numUncorruptedSamples) / numSamples - epsilon;

        final File intervalsFile = createTempFile("intervals", ".interval_list");
        final File annotatedIntervalsFile = createTempFile("annotated-intervals", ".tsv");
        final List<File> countFiles = new ArrayList<>(numSamples);
        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
            final String sampleName = String.format("sample-%d", sampleIndex);
            final SampleLocatableMetadata sampleMetadata = new SimpleSampleLocatableMetadata(sampleName, SEQUENCE_DICTIONARY);
            final List<SimpleCount> counts = new ArrayList<>(numIntervals);
            for (int intervalIndex = 0; intervalIndex < numIntervals; intervalIndex++) {
                final SimpleInterval interval = new SimpleInterval("20", intervalIndex + 1, intervalIndex + 1);
                if (intervalIndex < numIntervalsBelowCountThreshold && sampleIndex >= numUncorruptedSamples) {
                    //corrupt first numIntervalsBelowCountThreshold intervals in last (numSamples - numUncorruptedSamples) samples with low counts
                    counts.add(new SimpleCount(interval, lowCountFilterCountThreshold - 1));
                } else {
                    counts.add(new SimpleCount(interval, intervalIndex + lowCountFilterCountThreshold));
                }
            }
            final SimpleCountCollection sampleCounts = new SimpleCountCollection(sampleMetadata, counts);
            final File countFile = createTempFile(sampleName, ".tsv");
            sampleCounts.write(countFile);
            countFiles.add(countFile);

            if (sampleIndex == 0) {
                final List<SimpleInterval> sampleIntervals = sampleCounts.getIntervals();
                final IntervalList intervals = new IntervalList(sampleCounts.getMetadata().getSequenceDictionary());
                sampleIntervals.forEach(i -> intervals.add(new Interval(i)));
                intervals.write(intervalsFile);

                final AnnotationMap passingAnnotation = new AnnotationMap(Arrays.asList(
                        Pair.of(CopyNumberAnnotations.GC_CONTENT, 0.5),
                        Pair.of(CopyNumberAnnotations.MAPPABILITY, 0.5),
                        Pair.of(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT, 0.5)));
                final AnnotationMap failingAnnotation = new AnnotationMap(Arrays.asList(
                        Pair.of(CopyNumberAnnotations.GC_CONTENT, 0.05),
                        Pair.of(CopyNumberAnnotations.MAPPABILITY, 0.95),
                        Pair.of(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT, 0.05)));
                final List<AnnotatedInterval> annotatedIntervals = IntStream.range(0, numIntervals).boxed()
                        .map(i -> new AnnotatedInterval(new SimpleInterval(sampleIntervals.get(i)),
                            //last numIntervalsFailingAnnotationFilters intervals fail annotation-based filters
                            i >= numIntervals - numIntervalsFailingAnnotationFilters
                                    ? failingAnnotation
                                    : passingAnnotation))
                        .collect(Collectors.toList());
                new AnnotatedIntervalCollection(sampleMetadata, annotatedIntervals).write(annotatedIntervalsFile);
            }
        }

        return new Object[][]{
                //intervals file, annotated-intervals file, min/max GC content, mix/max mappability, min/max seg-dupe content,
                //count files, lowCountFilterCountThreshold, lowCountFilterPercentageOfSamples,
                //extremeCountFilterMinimumPercentile, extremeCountFilterMaximumPercentile, extremeCountFilterPercentageOfSamples,
                //expected array of indices of retained intervals
                {intervalsFile, annotatedIntervalsFile, 0.1, 0.9, 0.1, 0.9, 0.1, 0.9,
                        countFiles, lowCountFilterCountThreshold, percentageOfSamples, 1., 99., percentageOfSamples,
                        IntStream.range(numIntervalsBelowCountThreshold + 1, numIntervalsBelowCountThreshold + 99).boxed().collect(Collectors.toList())},
                {intervalsFile, annotatedIntervalsFile, 0.1, 0.9, 0.1, 0.9, 0.1, 0.9,
                        countFiles, lowCountFilterCountThreshold, percentageOfSamples, 1., 90., percentageOfSamples,
                        IntStream.range(numIntervalsBelowCountThreshold + 1, numIntervalsBelowCountThreshold + 90).boxed().collect(Collectors.toList())}};
    }

    @Test(dataProvider = "dataAllFilters")
    public void testAllFilters(final File intervalsFile,
                               final File annotatedIntervalsFile,
                               final double minimumGCContent,
                               final double maximumGCContent,
                               final double minimumMappability,
                               final double maximumMappability,
                               final double minimumSegmentalDuplicationContent,
                               final double maximumSegmentalDuplicationContent,
                               final List<File> countFiles,
                               final int lowCountFilterCountThreshold,
                               final double lowCountFilterPercentageOfSamples,
                               final double extremeCountFilterMinimumPercentile,
                               final double extremeCountFilterMaximumPercentile,
                               final double extremeCountFilterPercentageOfSamples,
                               final List<Integer> expectedIndices) {
        final File outputFile = createTempFile("filter-intervals-test", ".interval_list");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, intervalsFile.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_LONG_NAME, annotatedIntervalsFile.getAbsolutePath())
                .addArgument(FilterIntervals.MINIMUM_GC_CONTENT_LONG_NAME, Double.toString(minimumGCContent))
                .addArgument(FilterIntervals.MAXIMUM_GC_CONTENT_LONG_NAME, Double.toString(maximumGCContent))
                .addArgument(FilterIntervals.MINIMUM_MAPPABILITY_LONG_NAME, Double.toString(minimumMappability))
                .addArgument(FilterIntervals.MAXIMUM_MAPPABILITY_LONG_NAME, Double.toString(maximumMappability))
                .addArgument(FilterIntervals.MINIMUM_SEGMENTAL_DUPLICATION_CONTENT_LONG_NAME, Double.toString(minimumSegmentalDuplicationContent))
                .addArgument(FilterIntervals.MAXIMUM_SEGMENTAL_DUPLICATION_CONTENT_LONG_NAME, Double.toString(maximumSegmentalDuplicationContent))
                .addArgument(FilterIntervals.LOW_COUNT_FILTER_COUNT_THRESHOLD_LONG_NAME, Integer.toString(lowCountFilterCountThreshold))
                .addArgument(FilterIntervals.LOW_COUNT_FILTER_PERCENTAGE_OF_SAMPLES_LONG_NAME, Double.toString(lowCountFilterPercentageOfSamples))
                .addArgument(FilterIntervals.EXTREME_COUNT_FILTER_MINIMUM_PERCENTILE_LONG_NAME, Double.toString(extremeCountFilterMinimumPercentile))
                .addArgument(FilterIntervals.EXTREME_COUNT_FILTER_MAXIMUM_PERCENTILE_LONG_NAME, Double.toString(extremeCountFilterMaximumPercentile))
                .addArgument(FilterIntervals.EXTREME_COUNT_FILTER_PERCENTAGE_OF_SAMPLES_LONG_NAME, Double.toString(extremeCountFilterPercentageOfSamples))
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addOutput(outputFile);
        countFiles.forEach(argsBuilder::addInput);
        runCommandLine(argsBuilder);
        final IntervalList result = IntervalList.fromFile(outputFile);
        final IntervalList all = IntervalList.fromFile(intervalsFile);
        final List<Interval> allIntervals = all.getIntervals();
        final IntervalList expected = new IntervalList(all.getHeader().getSequenceDictionary());
        expectedIndices.stream().map(allIntervals::get).map(Interval::new).forEach(expected::add);
        Assert.assertEquals(result, expected);
        Assert.assertNotSame(result, expected);
    }
}