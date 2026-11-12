package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Merges the records of HaplotypeCaller's MAPQ-0 rescue (see {@link Mq0RescueArgumentCollection}) into its regular
 * output, so that the result is one coordinate-sorted VCF or GVCF in which no rescue record overlaps another record.
 *
 * <p>Regions are added in coordinate order, each with its regular records (as the regular calling path returns them:
 * calls and, in reference-confidence mode, one reference-confidence record per base for the rest, and possibly
 * non-call records such as GVCF hom-ref records or the monomorphic records of --output-mode EMIT_ALL_ACTIVE_SITES)
 * and its rescue records. A rescue record is decided once every regular record within the rescue distance of it is
 * known. It is dropped</p>
 * <ul>
 *     <li>if it overlaps a regular call (see {@link Mq0RescueUtils#isCall}), which wins;</li>
 *     <li>if it overlaps a regular non-call record that starts before it, which may already have been written;</li>
 *     <li>if a regular call lies within the rescue distance of it;</li>
 *     <li>or if it overlaps an earlier rescue record.</li>
 * </ul>
 * <p>Otherwise it takes the place of the reference-confidence and non-call records over its span. A non-call record
 * that starts before it is not replaced, since that would leave the bases before the rescue record without a record.
 * Regular records are released once no rescue record can still be placed before them, so they come out in the same
 * order as without rescue, and a GVCF writer downstream splits its reference blocks around the rescue records.</p>
 *
 * <p>So the rescue records written depend on the regular records the output holds. In GVCF and BP_RESOLUTION modes it
 * holds regular records that VCF mode does not write: non-call records (hom-ref or no-call records at the sites the
 * genotyper considered, such as a 0/0 deletion) and calls with QUAL below --standard-min-confidence-threshold-for-calling.
 * A rescue record that such a call overlaps or lies near, or that such a non-call record starting before it overlaps,
 * is written in VCF mode only.</p>
 */
final class Mq0RescueMerger {

    private final SAMSequenceDictionary dictionary;
    private final int distance;
    private final double minCallQual;

    // regular records not yet released, with the rescue records merged among them, in coordinate order
    private final List<VariantContext> pending = new ArrayList<>();
    // the last region and regular record added, to check that they come in coordinate order
    private Locatable lastRegion = null;
    private VariantContext lastRegular = null;
    // rescue records not yet decided, in coordinate order
    private final ArrayDeque<VariantContext> held = new ArrayDeque<>();
    // regular records other than reference-confidence ones that may overlap or lie near a held or future rescue record
    private final ArrayDeque<VariantContext> regularRecords = new ArrayDeque<>();
    private VariantContext lastMerged = null;

    private long merged = 0;
    private long replaced = 0;
    private long droppedOverlap = 0;
    private long droppedDistance = 0;

    /**
     * @param distance a rescue record is dropped if a regular call lies within this many bp
     * @param minCallQual the minimum QUAL for a regular record to count as a call
     */
    Mq0RescueMerger(final SAMSequenceDictionary dictionary, final int distance, final double minCallQual) {
        Utils.validateArg(distance >= 0, "distance must not be negative");
        this.dictionary = Utils.nonNull(dictionary);
        this.distance = distance;
        this.minCallQual = minCallQual;
    }

    /**
     * @return whether a record is a reference-confidence record: its only alternate allele is {@code <NON_REF>}
     */
    static boolean isReferenceConfidence(final VariantContext vc) {
        return vc.getNAlleles() == 2 && vc.getAlternateAllele(0).equals(Allele.NON_REF_ALLELE);
    }

    /**
     * Adds a region's records.
     *
     * @param region the region, which must not start before any region added before it
     * @param regular the region's regular records, in coordinate order
     * @param rescue the region's rescue records, in coordinate order, each starting within the region
     * @return the records that can now be written, in coordinate order
     */
    List<VariantContext> add(final Locatable region, final List<VariantContext> regular, final List<VariantContext> rescue) {
        Utils.validate(lastRegion == null || IntervalUtils.compareLocatables(lastRegion, region, dictionary) <= 0,
                () -> "region " + region + " added after region " + lastRegion);
        lastRegion = region;
        for (final VariantContext vc : regular) {
            final VariantContext previous = lastRegular;
            Utils.validate(previous == null || compare(previous, vc) <= 0, () -> "regular record " + vc + " added after " + previous);
            lastRegular = vc;
            pending.add(vc);
            if (!isReferenceConfidence(vc)) {
                regularRecords.add(vc);
            }
        }
        for (final VariantContext vc : rescue) {
            Utils.validate(vc.getContig().equals(region.getContig()) && vc.getStart() >= region.getStart() && vc.getStart() <= region.getEnd(),
                    () -> "rescue record " + vc + " does not start in region " + region);
            Utils.validate(held.isEmpty() || compare(held.peekLast(), vc) <= 0, () -> "rescue record out of order: " + vc);
            held.add(vc);
        }
        // decide the rescue records that no regular record still to come can lie near
        while (!held.isEmpty() && (!held.peek().getContig().equals(region.getContig())
                || (long) held.peek().getEnd() + distance <= region.getEnd())) {
            decide(held.poll());
        }
        // forget regular records too far upstream of every held and future rescue record to matter
        final int keepFrom = held.isEmpty() ? region.getEnd() + 1 : held.peek().getStart();
        while (!regularRecords.isEmpty() && (!regularRecords.peek().getContig().equals(region.getContig())
                || (long) regularRecords.peek().getEnd() + distance < keepFrom)) {
            regularRecords.poll();
        }
        return release(held.peek());
    }

    /**
     * @return the records not yet written, after deciding every rescue record still held
     */
    List<VariantContext> finish() {
        while (!held.isEmpty()) {
            decide(held.poll());
        }
        return release(null);
    }

    private void decide(final VariantContext rescue) {
        boolean nearCall = false;
        for (final VariantContext vc : regularRecords) {
            if (!vc.getContig().equals(rescue.getContig())) {
                continue;
            }
            final boolean call = Mq0RescueUtils.isCall(vc, minCallQual);
            if (vc.overlaps(rescue) && (call || vc.getStart() < rescue.getStart())) {
                droppedOverlap++;
                return;
            }
            nearCall |= call && vc.withinDistanceOf(rescue, distance);
        }
        if (nearCall) {
            droppedDistance++;
            return;
        }
        if (lastMerged != null && lastMerged.overlaps(rescue)) {
            droppedOverlap++;
            return;
        }
        merge(rescue);
        lastMerged = rescue;
        merged++;
    }

    /**
     * Puts a rescue record among the pending records, in place of the records it covers: reference-confidence records
     * and regular records that are not calls, all starting within it.
     */
    private void merge(final VariantContext rescue) {
        for (final Iterator<VariantContext> it = pending.iterator(); it.hasNext(); ) {
            final VariantContext vc = it.next();
            if (vc.overlaps(rescue)) {
                Utils.validate(vc.getStart() >= rescue.getStart() && !Mq0RescueUtils.isCall(vc, minCallQual),
                        () -> "unexpected record " + vc + " under rescue record " + rescue);
                if (!isReferenceConfidence(vc)) {
                    regularRecords.remove(vc);
                    replaced++;
                }
                it.remove();
            }
        }
        // after the removals, so that the index is one of the remaining records or the end
        int index = pending.size();
        for (int i = 0; i < pending.size(); i++) {
            if (compare(pending.get(i), rescue) > 0) {
                index = i;
                break;
            }
        }
        pending.add(index, rescue);
    }

    /**
     * Removes and returns the pending records that start before {@code firstHeld}, or all of them if it is null.
     */
    private List<VariantContext> release(final VariantContext firstHeld) {
        int n = 0;
        while (n < pending.size() && (firstHeld == null || compare(pending.get(n), firstHeld) < 0)) {
            n++;
        }
        final List<VariantContext> released = new ArrayList<>(pending.subList(0, n));
        pending.subList(0, n).clear();
        return released;
    }

    /** Compares the starts of two records by contig (in dictionary order) and position. */
    private int compare(final VariantContext a, final VariantContext b) {
        return IntervalUtils.compareLocatables(start(a), start(b), dictionary);
    }

    private static SimpleInterval start(final VariantContext vc) {
        return new SimpleInterval(vc.getContig(), vc.getStart(), vc.getStart());
    }

    /** @return the number of rescue records merged into the output */
    long merged() {
        return merged;
    }

    /** @return the number of regular non-call records that merged rescue records took the place of */
    long replaced() {
        return replaced;
    }

    /** @return the number of rescue records dropped because they overlap a regular record or an earlier rescue record */
    long droppedOverlap() {
        return droppedOverlap;
    }

    /** @return the number of rescue records dropped because a regular call lies within the rescue distance */
    long droppedDistance() {
        return droppedDistance;
    }
}
