package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.collect.Comparators;
import htsjdk.samtools.CigarElement;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

/**
 * Helpers for HaplotypeCaller's MAPQ-0 rescue (see {@link Mq0RescueArgumentCollection}): counting the alleles that
 * MAPQ-0 reads carry, finding where they carry one often enough to genotype the region again and which alleles to
 * genotype, and checking which rescue calls those reads support.
 */
final class Mq0RescueUtils {

    private Mq0RescueUtils() { }

    // the most distinct sets of alleles carried by single reads that multi-event sets are made from, which bounds the
    // work where many MAPQ-0 reads carry many alleles
    private static final int MAX_READ_SETS = 64;

    /**
     * An allele a read carries, as a reference and an alternate allele at a position. A substitution covers one base;
     * an indel is placed at the base before it (its anchor) and left-aligned against the reference, so the same indel
     * has one representation however the aligner placed it. The left-alignment ignores the other alleles the read
     * carries: see {@link AlleleCounts} for what that means for phasing.
     */
    record ReadAllele(int position, String ref, String alt) implements Comparable<ReadAllele> {
        private static final Comparator<ReadAllele> ORDER = Comparator.comparingInt(ReadAllele::position)
                .thenComparing(ReadAllele::ref).thenComparing(ReadAllele::alt);

        @Override
        public int compareTo(final ReadAllele other) {
            return ORDER.compare(this, other);
        }

        /** @return the last reference position the allele covers */
        int end() {
            return position + ref.length() - 1;
        }

        boolean overlaps(final ReadAllele other) {
            return position <= other.end() && other.position <= end();
        }

        /** @return whether both alleles are made of A, C, G and T only */
        boolean isRegular() {
            return BaseUtils.isAllRegularBases(ref.getBytes()) && BaseUtils.isAllRegularBases(alt.getBytes());
        }

        Event toEvent(final String contig) {
            return new Event(contig, position, Allele.create(ref, true), Allele.create(alt, false));
        }
    }

    /**
     * Reference bases of a window, with the left-alignment of indels.
     */
    private record Reference(byte[] bases, int start, int end) {
        static Reference of(final ReferenceContext referenceContext) {
            return new Reference(referenceContext.getBases(), referenceContext.getWindow().getStart(), referenceContext.getWindow().getEnd());
        }

        byte base(final int position) {
            return (byte) Character.toUpperCase(bases[position - start]);
        }

        String bases(final int from, final int length) {
            return new String(bases, from - start, length).toUpperCase();
        }

        /**
         * @return a deletion of the {@code length} bases after {@code anchor}, moved as far left as the reference allows
         * (within the window), or null if it does not fit in the window
         */
        ReadAllele deletion(int anchor, final int length) {
            if (anchor < start || anchor + length > end) {
                return null;
            }
            while (anchor > start && base(anchor) == base(anchor + length)) {
                anchor--;
            }
            return new ReadAllele(anchor, bases(anchor, length + 1), bases(anchor, 1));
        }

        /**
         * @return an insertion of {@code inserted} after {@code anchor}, moved as far left as the reference allows
         * (within the window), or null if the anchor is outside the window
         */
        ReadAllele insertion(int anchor, String inserted) {
            if (anchor < start || anchor > end) {
                return null;
            }
            inserted = inserted.toUpperCase();
            while (anchor > start && base(anchor) == inserted.charAt(inserted.length() - 1)) {
                inserted = (char) base(anchor) + inserted.substring(0, inserted.length() - 1);
                anchor--;
            }
            final String anchorBase = bases(anchor, 1);
            return new ReadAllele(anchor, anchorBase, anchorBase + inserted);
        }

        /**
         * @return the allele {@code alt} of a record at {@code position} with reference allele {@code ref}, trimmed to
         * its minimal representation and, if it is an indel, left-aligned; or null if it is not a plain allele
         */
        ReadAllele normalize(final int position, final Allele ref, final Allele alt) {
            if (!isRealAllele(alt)) {
                return null;
            }
            // trim the common suffix, keeping at least one base of each allele (Event does not, for padded
            // multi-allelic records such as REF TTT ALT TT)
            String r = ref.getBaseString().toUpperCase();
            String a = alt.getBaseString().toUpperCase();
            while (r.length() > 1 && a.length() > 1 && r.charAt(r.length() - 1) == a.charAt(a.length() - 1)) {
                r = r.substring(0, r.length() - 1);
                a = a.substring(0, a.length() - 1);
            }
            if (r.length() == a.length() || r.charAt(0) != a.charAt(0) || Math.min(r.length(), a.length()) != 1) {
                return new ReadAllele(position, r, a);
            }
            return r.length() > a.length() ? deletion(position, r.length() - 1) : insertion(position, a.substring(1));
        }
    }

    /**
     * The alleles a set of reads carries over a span, as found by walking each read's CIGAR once. A base is counted
     * only if its quality is at least the minimum; indels are counted whatever the qualities of the bases around them.
     * Reads are counted, not fragments: the two mates of a short fragment that both cover a position count twice, as
     * they do in HaplotypeCaller's DP. Depths, trigger positions and the alleles each read carries (for phasing) are
     * kept for the span only; alleles the reads carry elsewhere in the reference window are counted too, but only for
     * {@link #supportsCalledAlleles}, so that a call near the end of the span (an MNP, or an indel that left-aligns to
     * before the span) is not rejected for a part of it outside the span.
     *
     * <p>Each indel is left-aligned against the reference alone, as the genotyper's calls are, so a read and a call
     * carrying the same indel give the same allele. If the shift moves an indel past another allele the same read
     * carries (say, an aligner put a deletion to the right of the read's SNV in a homopolymer), the shifted indel and
     * that allele together make a haplotype the read does not carry. Such a read still counts for each of its alleles,
     * but not for the multi-event sets of {@link #eventSets}.</p>
     *
     * <p>The counting is written out rather than built on {@link org.broadinstitute.hellbender.utils.pileup.PileupBasedAlleles}
     * or read pileups, which give only alleles that pass their own activity and quality gates, without raw counts,
     * the alleles each read carries (needed for phasing) or left-aligned indels. The indels are left-aligned here
     * against the window's bases, which is cheaper per read than building a VariantContext for
     * {@link org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils#leftAlignAndTrim} or the allele
     * arrays {@link org.broadinstitute.hellbender.utils.read.AlignmentUtils#normalizeAlleles} takes, and never reaches
     * outside the window.</p>
     */
    static final class AlleleCounts {
        private final SimpleInterval span;
        private final Reference reference;
        // the reads with a base or a deletion counted at each position of the span
        private final int[] depth;
        // the number of reads carrying each allele anchored in the span
        private final Map<ReadAllele, Integer> counts = new HashMap<>();
        // the number of reads carrying each allele anchored outside the span, for the support check only
        private final Map<ReadAllele, Integer> outsideCounts = new HashMap<>();
        // the alleles anchored in the span that each read carries, in order, for the reads carrying any whose indels
        // were not left-aligned past another of their alleles
        private final List<List<ReadAllele>> perRead = new ArrayList<>();

        AlleleCounts(final Collection<GATKRead> reads, final SimpleInterval span, final ReferenceContext referenceContext, final int minBaseQuality) {
            this.span = Utils.nonNull(span);
            this.reference = Reference.of(referenceContext);
            this.depth = new int[span.size()];
            for (final GATKRead read : reads) {
                if (read.overlaps(span)) {
                    addRead(read, minBaseQuality);
                }
            }
        }

        private boolean inSpan(final int position) {
            return position >= span.getStart() && position <= span.getEnd();
        }

        private void addRead(final GATKRead read, final int minBaseQuality) {
            final byte[] bases = read.getBasesNoCopy();
            final byte[] quals = read.getBaseQualitiesNoCopy();
            final List<ReadAllele> carried = new ArrayList<>(2);
            final List<ReadAllele> outside = new ArrayList<>(0);
            boolean phaseable = true;
            int refPos = read.getStart();
            int readPos = 0;
            for (final CigarElement e : read.getCigarElements()) {
                final int len = e.getLength();
                switch (e.getOperator()) {
                    case M: case EQ: case X:
                        for (int k = 0; k < len; k++) {
                            final int pos = refPos + k;
                            if (quals[readPos + k] < minBaseQuality) {
                                continue;
                            }
                            final boolean in = inSpan(pos);
                            if (in) {
                                depth[pos - span.getStart()]++;
                            }
                            final byte b = (byte) Character.toUpperCase(bases[readPos + k]);
                            if (pos >= reference.start() && pos <= reference.end()) {
                                final byte r = reference.base(pos);
                                if (BaseUtils.isRegularBase(b) && BaseUtils.isRegularBase(r) && b != r) {
                                    (in ? carried : outside).add(new ReadAllele(pos, String.valueOf((char) r), String.valueOf((char) b)));
                                }
                            }
                        }
                        refPos += len;
                        readPos += len;
                        break;
                    case D:
                        for (int k = 0; k < len; k++) {
                            if (inSpan(refPos + k)) {
                                depth[refPos + k - span.getStart()]++;
                            }
                        }
                        phaseable &= addIndel(carried, outside, reference.deletion(refPos - 1, len), refPos - 1);
                        refPos += len;
                        break;
                    case I:
                        if (readPos + len <= bases.length) {
                            phaseable &= addIndel(carried, outside, reference.insertion(refPos - 1, new String(bases, readPos, len)), refPos - 1);
                        }
                        readPos += len;
                        break;
                    case N:
                        refPos += len;
                        break;
                    case S:
                        readPos += len;
                        break;
                    default:
                        break;
                }
            }
            if (!carried.isEmpty()) {
                carried.sort(Comparator.naturalOrder());
                if (phaseable) {
                    perRead.add(carried);
                }
                for (final ReadAllele allele : new LinkedHashSet<>(carried)) {
                    counts.merge(allele, 1, Integer::sum);
                }
            }
            for (final ReadAllele allele : new LinkedHashSet<>(outside)) {
                outsideCounts.merge(allele, 1, Integer::sum);
            }
        }

        /**
         * Adds a read's left-aligned indel, as aligned at {@code alignedAnchor}, to the alleles it carries in or
         * outside the span.
         *
         * @return false if the left-alignment moved it past another allele the read carries in the span
         */
        private boolean addIndel(final List<ReadAllele> carried, final List<ReadAllele> outside, final ReadAllele indel, final int alignedAnchor) {
            if (indel == null) {
                return true;
            }
            (inSpan(indel.position()) ? carried : outside).add(indel);
            return carried.stream().noneMatch(other -> other != indel && other.position() > indel.position() && other.position() <= alignedAnchor);
        }

        /** @return the number of reads with a base of sufficient quality, or a deletion, at a position of the span */
        int depth(final int position) {
            return depth[position - span.getStart()];
        }

        /** @return the number of reads carrying {@code allele}, in the span or elsewhere in the reference window */
        int count(final ReadAllele allele) {
            return counts.getOrDefault(allele, 0) + outsideCounts.getOrDefault(allele, 0);
        }

        /**
         * @return the positions of the span, in order, at which at least {@code minAltReads} reads carry the same
         * allele, and those reads make up at least {@code minAltFraction} of the reads counted there
         */
        List<Integer> altPositions(final int minAltReads, final double minAltFraction) {
            final Map<Integer, Integer> best = new TreeMap<>();
            counts.forEach((allele, n) -> best.merge(allele.position(), n, Math::max));
            return best.entrySet().stream()
                    .filter(en -> en.getValue() >= minAltReads && en.getValue() >= minAltFraction * Math.max(depth[en.getKey() - span.getStart()], 1))
                    .map(Map.Entry::getKey)
                    .collect(Collectors.toList());
        }

        /**
         * @return the alleles at {@code positions} with only A, C, G and T that at least {@code minAltReads} reads carry,
         * at most {@code maxAlleles} of them (the most frequent): the alleles {@link #eventSets} is given
         */
        Set<ReadAllele> selectAlleles(final Collection<Integer> positions, final int minAltReads, final int maxAlleles) {
            final Set<Integer> wanted = new HashSet<>(positions);
            return counts.entrySet().stream()
                    .filter(en -> wanted.contains(en.getKey().position()) && en.getValue() >= minAltReads && en.getKey().isRegular())
                    .sorted(Map.Entry.<ReadAllele, Integer>comparingByValue().reversed().thenComparing(Map.Entry.comparingByKey()))
                    .limit(maxAlleles)
                    .map(Map.Entry::getKey)
                    .collect(Collectors.toSet());
        }

        /**
         * The event sets to insert into the reference haplotype for the alleles {@link #selectAlleles} selected at the
         * trigger positions. First, one single-event set per selected allele, in position order. Then up to
         * {@code maxEventSets} sets of two or more of those alleles, not overlapping, that at least {@code minAltReads}
         * reads each carry all of (a read may carry more), largest first and then most frequent, so that where reads
         * carry many alleles the haplotypes phase them as fully as the reads allow. The sets tried are those the reads
         * carry and the alleles any two of those have in common, so that reads of different extents add up.
         *
         * @return the event sets, or an empty list if no allele is selected
         */
        List<List<Event>> eventSets(final Set<ReadAllele> selected, final int minAltReads, final int maxEventSets) {
            final List<List<Event>> eventSets = selected.stream().sorted()
                    .map(allele -> List.of(allele.toEvent(span.getContig())))
                    .collect(Collectors.toList());

            // the selected alleles each read carries, if two or more and not overlapping; the most common of these sets
            final Map<List<ReadAllele>, Integer> readSets = new HashMap<>();
            for (final List<ReadAllele> carried : perRead) {
                final List<ReadAllele> set = carried.stream().filter(selected::contains).distinct().collect(Collectors.toList());
                if (set.size() >= 2 && !hasOverlap(set)) {
                    readSets.merge(set, 1, Integer::sum);
                }
            }
            final List<List<ReadAllele>> distinct = readSets.entrySet().stream()
                    .sorted(Map.Entry.<List<ReadAllele>, Integer>comparingByValue().reversed().thenComparing(Map.Entry::getKey, LEXICOGRAPHIC))
                    .limit(MAX_READ_SETS)
                    .map(Map.Entry::getKey)
                    .collect(Collectors.toList());
            if (distinct.isEmpty() || maxEventSets == 0) {
                return eventSets;
            }
            // candidates: those sets and the pairwise intersections of them, counted by the reads carrying all of them,
            // overlapping or not; a read with fewer than two selected alleles carries no candidate
            final List<Set<ReadAllele>> perReadSets = perRead.stream()
                    .map(carried -> carried.stream().filter(selected::contains).collect(Collectors.toSet()))
                    .filter(set -> set.size() >= 2)
                    .collect(Collectors.toList());
            final Set<List<ReadAllele>> candidates = new HashSet<>(distinct);
            for (int i = 0; i < distinct.size(); i++) {
                for (int j = i + 1; j < distinct.size(); j++) {
                    final List<ReadAllele> common = new ArrayList<>(distinct.get(i));
                    common.retainAll(distinct.get(j));
                    if (common.size() >= 2) {
                        candidates.add(common);
                    }
                }
            }
            final Map<List<ReadAllele>, Integer> carriedBy = new HashMap<>();
            for (final List<ReadAllele> candidate : candidates) {
                final int n = (int) perReadSets.stream().filter(carried -> carried.containsAll(candidate)).count();
                if (n >= minAltReads) {
                    carriedBy.put(candidate, n);
                }
            }
            carriedBy.entrySet().stream()
                    .sorted(Comparator.<Map.Entry<List<ReadAllele>, Integer>>comparingInt(en -> en.getKey().size()).reversed()
                            .thenComparing(Map.Entry.<List<ReadAllele>, Integer>comparingByValue().reversed())
                            .thenComparing(Map.Entry::getKey, LEXICOGRAPHIC))
                    .limit(maxEventSets)
                    .forEach(en -> eventSets.add(en.getKey().stream().map(allele -> allele.toEvent(span.getContig())).collect(Collectors.toList())));
            return eventSets;
        }

        /**
         * @return whether every alternate allele that a called genotype of {@code vc} carries is carried by at least one
         * of the reads, and there is at least one such allele. Each allele is checked on its own (an MNP base by base),
         * so a pair of alleles is accepted even if no single read carries both.
         */
        boolean supportsCalledAlleles(final VariantContext vc) {
            final Set<Allele> called = new LinkedHashSet<>();
            for (final Genotype g : vc.getGenotypes()) {
                for (final Allele a : g.getAlleles()) {
                    if (a.isNonReference() && isRealAllele(a)) {
                        called.add(a);
                    }
                }
            }
            if (called.isEmpty()) {
                return false;
            }
            for (final Allele alt : called) {
                final ReadAllele allele = reference.normalize(vc.getStart(), vc.getReference(), alt);
                if (allele == null || substitutions(allele).stream().anyMatch(a -> count(a) == 0)) {
                    return false;
                }
            }
            return true;
        }

        /**
         * @return the single-base substitutions an allele of equal-length reference and alternate alleles (an MNP, as
         * --max-mnp-distance makes them) is made of, which is how reads carry it; or the allele itself otherwise
         */
        private static List<ReadAllele> substitutions(final ReadAllele allele) {
            if (allele.ref().length() != allele.alt().length() || allele.ref().length() == 1) {
                return List.of(allele);
            }
            final List<ReadAllele> result = new ArrayList<>();
            for (int i = 0; i < allele.ref().length(); i++) {
                if (allele.ref().charAt(i) != allele.alt().charAt(i)) {
                    result.add(new ReadAllele(allele.position() + i, allele.ref().substring(i, i + 1), allele.alt().substring(i, i + 1)));
                }
            }
            return result;
        }
    }

    private static boolean hasOverlap(final List<ReadAllele> sorted) {
        for (int i = 1; i < sorted.size(); i++) {
            if (sorted.get(i).overlaps(sorted.get(i - 1))) {
                return true;
            }
        }
        return false;
    }

    private static final Comparator<Iterable<ReadAllele>> LEXICOGRAPHIC = Comparators.lexicographical(Comparator.<ReadAllele>naturalOrder());

    /**
     * @return the positions of {@code positions} that are more than {@code distance} bp from every call of {@code calls}
     */
    static List<Integer> awayFrom(final List<Integer> positions, final String contig, final Collection<VariantContext> calls, final int distance) {
        return positions.stream()
                .filter(pos -> calls.stream().noneMatch(vc -> vc.withinDistanceOf(new SimpleInterval(contig, pos, pos), distance)))
                .collect(Collectors.toList());
    }

    /** @return the number of reads (not fragments) overlapping {@code contig:pos} */
    static int countOverlapping(final Collection<GATKRead> reads, final String contig, final int pos) {
        final SimpleInterval site = new SimpleInterval(contig, pos, pos);
        return (int) reads.stream().filter(read -> read.overlaps(site)).count();
    }

    /**
     * @return whether a record is a call: it has a non-symbolic alternate allele, a called genotype carrying a
     * non-reference allele, and QUAL at least {@code minQual}. GVCF reference blocks, hom-ref and no-call records are not.
     */
    static boolean isCall(final VariantContext vc, final double minQual) {
        if (!hasRealAlt(vc) || vc.getPhredScaledQual() < minQual) {
            return false;
        }
        for (final Genotype g : vc.getGenotypes()) {
            for (final Allele a : g.getAlleles()) {
                if (a.isNonReference() && isRealAllele(a)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * @return whether a record has an alternate allele that is not symbolic or a spanning deletion, i.e. it is not a
     * reference-confidence record
     */
    static boolean hasRealAlt(final VariantContext vc) {
        return vc.getAlternateAlleles().stream().anyMatch(Mq0RescueUtils::isRealAllele);
    }

    /** @return whether an allele is called, not symbolic and not a spanning deletion */
    private static boolean isRealAllele(final Allele a) {
        return !a.isNoCall() && !a.isSymbolic() && !Allele.SPAN_DEL.equals(a);
    }
}
