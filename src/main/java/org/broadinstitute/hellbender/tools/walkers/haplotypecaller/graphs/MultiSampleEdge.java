package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Edge class for connecting nodes in the graph that tracks some per-sample information.
 * <p>
 * This class extends BaseEdge with the additional functionality of tracking the maximum
 * multiplicity seen within any single sample.  The workflow for using this class is:
 * </p>
 * <pre>
 * {@code
 *      MultiSampleEdge e = new MultiSampleEdge(ref, 1)
 *      e.incMultiplicity(1)              // total is 2, per sample is 2, max per sample is 1
 *      e.getPruningMultiplicity()        // = 1
 *      e.flushSingleSampleMultiplicity() // total is 2, per sample is 0, max per sample is 2
 *      e.getPruningMultiplicity()        // = 2
 *      e.incMultiplicity(3)              // total is 5, per sample is 3, max per sample is 2
 *      e.getPruningMultiplicity()        // = 2
 *      e.flushSingleSampleMultiplicity() // total is 5, per sample is 0, max per sample is 3
 *      e.getPruningMultiplicity()        // = 3
 * }
 * </pre>
 */
public final class MultiSampleEdge extends BaseEdge {
    // 2: single-sample multiplicities held in a sorted int array rather than a priority queue
    private static final long serialVersionUID = 2L;
    private int currentSingleSampleMultiplicity;
    private final int singleSampleCapacity;

    /**
     * The retained single-sample multiplicities in ascending order, in the first {@link #retainedCount} slots: the
     * {@code singleSampleCapacity} largest of the initial multiplicity and the flushed ones, held as primitives so an
     * edge costs one small array rather than a priority queue of boxed values. The extra slot lets a flushed value be
     * inserted in place, in order, before the lowest value is shifted out.
     */
    private final int[] retainedMultiplicities;
    private int retainedCount;

    /** Allocated on the first {@link #addReferenceIndex}, which only junction-tree assembly calls. */
    private List<Integer> referencePathIndexes = null;

    /**
     * Create a new MultiSampleEdge with weight multiplicity and, if isRef == true, indicates a path through the reference
     *
     * @param isRef indicates whether this edge is a path through the reference
     * @param multiplicity the number of observations of this edge in this sample
     * @param singleSampleCapacity the max number of samples to track edge multiplicities
     */
    public MultiSampleEdge(final boolean isRef, final int multiplicity, final int singleSampleCapacity) {
        super(isRef, multiplicity);

        Utils.validateArg( singleSampleCapacity > 0, () -> "singleSampleCapacity must be > 0 but found: " + singleSampleCapacity);
        retainedMultiplicities = new int[singleSampleCapacity + 1];
        retainedMultiplicities[0] = multiplicity;
        retainedCount = 1;
        currentSingleSampleMultiplicity = multiplicity;
        this.singleSampleCapacity = singleSampleCapacity;
    }

    /**
     * Create a new edge with this edge's reference flag, multiplicity and sample capacity, as if its whole
     * multiplicity had been seen in one sample: the per-sample multiplicities and reference path indexes are not
     * copied. {@link #duplicate()} copies them too.
     */
    @Override
    public MultiSampleEdge copy() {
        return new MultiSampleEdge(isRef(), getMultiplicity(), singleSampleCapacity);
    }

    /**
     * Create a new edge equal to this one in every respect: {@link #copy()}, plus the current and retained per-sample
     * multiplicities and a copy of the reference path indexes.
     */
    @Override
    public MultiSampleEdge duplicate() {
        final MultiSampleEdge duplicate = copy();
        duplicate.currentSingleSampleMultiplicity = currentSingleSampleMultiplicity;
        System.arraycopy(retainedMultiplicities, 0, duplicate.retainedMultiplicities, 0, retainedCount);
        duplicate.retainedCount = retainedCount;
        if (referencePathIndexes != null) {
            duplicate.referencePathIndexes = new ArrayList<>(referencePathIndexes);
        }
        return duplicate;
    }

    /**
     * update the single sample multiplicities by adding the current single sample multiplicity to the retained ones, and
     * reset the current single sample multiplicity to 0.
     */
    public void flushSingleSampleMultiplicity() {
        if (retainedCount == singleSampleCapacity && currentSingleSampleMultiplicity <= retainedMultiplicities[0]) {
            // the new value would be the lowest, or tie with it, and be dropped at once, leaving the same values
            currentSingleSampleMultiplicity = 0;
            return;
        }
        // insert in ascending order, after any equal values
        int i = retainedCount;
        while (i > 0 && retainedMultiplicities[i - 1] > currentSingleSampleMultiplicity) {
            retainedMultiplicities[i] = retainedMultiplicities[i - 1];
            i--;
        }
        retainedMultiplicities[i] = currentSingleSampleMultiplicity;
        retainedCount++;
        if (retainedCount == singleSampleCapacity + 1) {
            // remove the lowest multiplicity from the list
            System.arraycopy(retainedMultiplicities, 1, retainedMultiplicities, 0, singleSampleCapacity);
            retainedCount--;
        }
        currentSingleSampleMultiplicity = 0;
    }

    @Override
    public void incMultiplicity(final int incr) {
        super.incMultiplicity(incr);
        currentSingleSampleMultiplicity += incr;
    }

    @Override
    public int getPruningMultiplicity() {
        return retainedMultiplicities[0];
    }

    @Override
    public String getDotLabel() {
        return super.getDotLabel() + '/' + getPruningMultiplicity();
    }

    @VisibleForTesting
    int getCurrentSingleSampleMultiplicity() {
        return currentSingleSampleMultiplicity;
    }

    /**
     * Records that this edge is on the reference path at index {@code i}. The index list is allocated on the first call.
     */
    public void addReferenceIndex(int i) {
        if (referencePathIndexes == null) {
            referencePathIndexes = new ArrayList<>(2);
        }
        referencePathIndexes.add(i);
    }

    /**
     * @return the indexes of the reference path this edge is on, in the order they were added, as a read-only view
     *         (empty if none were added)
     */
    public List<Integer> getReferencePathIndexes() {
        return referencePathIndexes == null ? Collections.emptyList() : Collections.unmodifiableList(referencePathIndexes);
    }
}
