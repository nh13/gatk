package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;

public final class MultiSampleEdgeUnitTest extends GATKBaseTest {

    private class MultiplicityTestProvider {
        final List<Integer> countsPerSample;
        final int numSamplesPruning;
        public MultiplicityTestProvider(final List<Integer> countsPerSample, final int numSamplesPruning) {
            this.countsPerSample = countsPerSample;
            this.numSamplesPruning = numSamplesPruning;
        }
    }

    @DataProvider(name = "MultiplicityData")
    public Object[][] makeMultiplicityData() {
        List<Object[]> tests = new ArrayList<>();

        final List<Integer> countsPerSample = Arrays.asList(0, 1, 2, 3, 4, 5);
        for ( final int numSamplesPruning : Arrays.asList(1, 2, 3) ) {
            for ( final int nSamples : Arrays.asList(1, 2, 3, 4, 5)) {
                for ( final List<Integer> perm : Utils.makePermutations(countsPerSample, nSamples, false) ) {
                    tests.add(new Object[]{new MultiplicityTestProvider(perm, numSamplesPruning)});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    /**
     * Example testng test using MyDataProvider
     */
    @Test(dataProvider = "MultiplicityData")
    public void testMultiplicity(final MultiplicityTestProvider cfg) {
        final MultiSampleEdge edge = new MultiSampleEdge(false, 0, cfg.numSamplesPruning);
        Assert.assertEquals(edge.getMultiplicity(), 0);
        Assert.assertEquals(edge.getPruningMultiplicity(), 0);

        final MultiSampleEdge copy = edge.copy();
        Assert.assertEquals(edge.getCurrentSingleSampleMultiplicity(), copy.getCurrentSingleSampleMultiplicity());
        Assert.assertEquals(edge.getDotLabel(), copy.getDotLabel());
        Assert.assertEquals(edge.getPruningMultiplicity(), copy.getPruningMultiplicity());
        Assert.assertEquals(edge.getMultiplicity(), copy.getMultiplicity());
        Assert.assertEquals(edge.getClass(), copy.getClass());

        int total = 0;
        for ( int i = 0; i < cfg.countsPerSample.size(); i++ ) {
            int countForSample = 0;
            for ( int count = 0; count < cfg.countsPerSample.get(i); count++ ) {
                edge.incMultiplicity(1);
                total++;
                countForSample++;
                Assert.assertEquals(edge.getMultiplicity(), total);
                Assert.assertEquals(edge.getCurrentSingleSampleMultiplicity(), countForSample);
            }
            edge.flushSingleSampleMultiplicity();
        }

        ArrayList<Integer> counts = new ArrayList<>(cfg.countsPerSample);
        counts.add(0);
        Collections.sort(counts);
        final int prune = counts.get(Math.max(counts.size() - cfg.numSamplesPruning, 0));
        Assert.assertEquals(edge.getMultiplicity(), total);
        Assert.assertEquals(edge.getPruningMultiplicity(), prune);
    }

    /** An edge seen 3, 5 and 2 times in three samples, with a sample in progress and two reference path indexes. */
    private static MultiSampleEdge edgeWithSampleHistory() {
        final MultiSampleEdge edge = new MultiSampleEdge(true, 3, 2);
        edge.flushSingleSampleMultiplicity();
        edge.incMultiplicity(5);
        edge.flushSingleSampleMultiplicity();
        edge.incMultiplicity(2);
        edge.flushSingleSampleMultiplicity();
        edge.incMultiplicity(4);
        edge.addReferenceIndex(7);
        edge.addReferenceIndex(11);
        return edge;
    }

    @Test
    public void aDuplicateKeepsThePerSampleStateAndReferenceIndexes() {
        final MultiSampleEdge edge = edgeWithSampleHistory();
        final MultiSampleEdge duplicate = edge.duplicate();
        Assert.assertNotSame(duplicate, edge);
        Assert.assertEquals(duplicate.isRef(), edge.isRef());
        Assert.assertEquals(duplicate.getMultiplicity(), edge.getMultiplicity());
        Assert.assertEquals(duplicate.getPruningMultiplicity(), edge.getPruningMultiplicity());
        Assert.assertEquals(duplicate.getCurrentSingleSampleMultiplicity(), edge.getCurrentSingleSampleMultiplicity());
        Assert.assertEquals(duplicate.getReferencePathIndexes(), edge.getReferencePathIndexes());
        Assert.assertEquals(duplicate.getDotLabel(), edge.getDotLabel());
    }

    @Test
    public void aDuplicatesReferenceIndexesAreIndependentOfTheOriginals() {
        final MultiSampleEdge edge = edgeWithSampleHistory();
        final MultiSampleEdge duplicate = edge.duplicate();
        duplicate.addReferenceIndex(13);
        Assert.assertEquals(edge.getReferencePathIndexes(), List.of(7, 11));
        Assert.assertEquals(duplicate.getReferencePathIndexes(), List.of(7, 11, 13));
        edge.addReferenceIndex(17);
        Assert.assertEquals(edge.getReferencePathIndexes(), List.of(7, 11, 17));
        Assert.assertEquals(duplicate.getReferencePathIndexes(), List.of(7, 11, 13));

        // an edge with no indexes yet duplicates to one that can still take them, without giving them to the original
        final MultiSampleEdge bare = new MultiSampleEdge(false, 1, 1);
        final MultiSampleEdge bareDuplicate = bare.duplicate();
        bareDuplicate.addReferenceIndex(3);
        Assert.assertEquals(bareDuplicate.getReferencePathIndexes(), List.of(3));
        Assert.assertTrue(bare.getReferencePathIndexes().isEmpty());
    }

    @Test
    public void theReferenceIndexesAreReadOnlyWhetherOrNotAnyWereAdded() {
        final MultiSampleEdge edge = new MultiSampleEdge(false, 1, 1);
        Assert.assertThrows(UnsupportedOperationException.class, () -> edge.getReferencePathIndexes().add(1));
        edge.addReferenceIndex(2);
        Assert.assertThrows(UnsupportedOperationException.class, () -> edge.getReferencePathIndexes().add(1));
        Assert.assertEquals(edge.getReferencePathIndexes(), List.of(2));
    }

    @Test
    public void aDuplicateGoesOnToBehaveLikeTheOriginal() {
        final MultiSampleEdge edge = edgeWithSampleHistory();
        final MultiSampleEdge duplicate = edge.duplicate();
        for (final MultiSampleEdge e : List.of(edge, duplicate)) {
            e.flushSingleSampleMultiplicity();
            e.incMultiplicity(1);
            e.flushSingleSampleMultiplicity();
        }
        Assert.assertEquals(duplicate.getPruningMultiplicity(), edge.getPruningMultiplicity());
        Assert.assertEquals(duplicate.getMultiplicity(), edge.getMultiplicity());
    }

    @Test
    public void aCopyCountsItsWholeMultiplicityAsOneSample() {
        final MultiSampleEdge copy = edgeWithSampleHistory().copy();
        Assert.assertEquals(copy.getMultiplicity(), 14);
        Assert.assertEquals(copy.getPruningMultiplicity(), 14);
        Assert.assertEquals(copy.getCurrentSingleSampleMultiplicity(), 14);
        Assert.assertTrue(copy.getReferencePathIndexes().isEmpty());
    }

    /**
     * The pruning multiplicity follows a priority queue that holds the initial multiplicity and every flushed one,
     * dropping its lowest whenever it grows past the capacity, for any initial multiplicity and flush sequence.
     */
    @Test
    public void thePruningMultiplicityMatchesABoundedPriorityQueue() {
        final Random random = new Random(42);
        for (int trial = 0; trial < 2000; trial++) {
            final int capacity = 1 + random.nextInt(4);
            final int initial = random.nextInt(6);
            final MultiSampleEdge edge = new MultiSampleEdge(false, initial, capacity);
            final PriorityQueue<Integer> model = new PriorityQueue<>();
            model.add(initial);
            int current = initial;
            final int flushes = random.nextInt(8);
            for (int f = 0; f < flushes; f++) {
                final int increment = random.nextInt(7);
                edge.incMultiplicity(increment);
                current += increment;
                if (random.nextBoolean()) {
                    final MultiSampleEdge duplicate = edge.duplicate();
                    Assert.assertEquals(duplicate.getPruningMultiplicity(), (int) model.peek());
                }
                edge.flushSingleSampleMultiplicity();
                model.add(current);
                if (model.size() == capacity + 1) {
                    model.poll();
                }
                current = 0;
                Assert.assertEquals(edge.getPruningMultiplicity(), (int) model.peek());
            }
        }
    }
}
