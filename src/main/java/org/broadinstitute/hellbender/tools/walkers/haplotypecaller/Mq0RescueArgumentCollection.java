package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.Serializable;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Arguments for HaplotypeCaller's MAPQ-0 rescue.
 *
 * <p>The rescue has been validated in DRAGEN mode only, so {@code --mq0-rescue} requires --dragen-mode (or
 * --dragen-378-concordance-mode), and HaplotypeCaller fails if it is given without one of them.</p>
 *
 * <p>With {@code --mq0-rescue}, reads that fail the mapping-quality read filter only because their MAPQ is 0 are kept
 * aside instead of being discarded (at most --max-reads-per-alignment-start per alignment start, held until the regions
 * they overlap are called). Each region is first called exactly as without the flag. Where the MAPQ-0 reads carry the
 * same allele often enough away from the region's regular calls, the region is genotyped a second time with all reads,
 * against haplotypes made as {@code --mq0-rescue-mode} says: by default the reference haplotype with those alleles
 * inserted into it, or else the region assembled again. That genotyping uses likelihood engines that do not cap base
 * qualities at the read's MAPQ, for all of its reads, as DRAGEN mode does for every call, even if the cap is turned
 * back on. Records from it that call a non-reference genotype, with QUAL at least {@code --mq0-rescue-min-call-qual},
 * whose called alleles MAPQ-0 reads carry, are merged into the output with the FILTER
 * {@value GATKVCFConstants#MQ0_RESCUE_FILTER_NAME} and the INFO fields {@value GATKVCFConstants#MQ0_RESCUE_KEY} (a
 * flag), {@value GATKVCFConstants#MQ0_RESCUE_DEPTH_KEY} and {@value GATKVCFConstants#MQ0_RESCUE_MQ0_DEPTH_KEY}. A
 * rescue record is dropped if it overlaps a regular call or an earlier rescue record, if a regular call lies within
 * {@code --mq0-rescue-distance} bp, if it overlaps a regular non-call record that starts before it, or if it overlaps a
 * record given with --alleles that is force-called (the regular record force-called there wins, whatever its genotype);
 * a rescue call with a '*' allele is not made into a rescue record at all. Otherwise a rescue record takes the place of
 * the reference-confidence and non-call records over its span. The other records are the same as without the flag. If
 * the read filters let MAPQ-0 reads through, none is set aside for the rescue and the flag has no effect (with a
 * warning).</p>
 *
 * <p>The rescue records of GVCF (or BP_RESOLUTION) and VCF output differ in three ways. Rescue calls with QUAL below
 * --standard-min-confidence-threshold-for-calling are made in GVCF mode only. GVCF output holds regular records VCF
 * output leaves out (calls with QUAL below that threshold, and non-call records such as a 0/0 deletion), which drop
 * the rescue records they overlap or, for calls, lie near. And a rescue call carrying a '*' allele, for an upstream
 * deletion the rescue genotyper saw, is dropped; the two modes keep different deletions, so this can happen in one
 * and not the other. So, but for those drops, a rescue record of the GVCF with QUAL at least that threshold is in
 * the VCF too, while the VCF can hold rescue records the GVCF does not.</p>
 *
 * <p>Read counts here count reads, not fragments: overlapping mates count twice. GenotypeGVCFs and CombineGVCFs keep
 * {@value GATKVCFConstants#MQ0_RESCUE_DEPTH_KEY} and {@value GATKVCFConstants#MQ0_RESCUE_MQ0_DEPTH_KEY} but drop the
 * FILTER and the {@value GATKVCFConstants#MQ0_RESCUE_KEY} flag, and GenotypeGVCFs on a GenomicsDB workspace keeps none
 * of them, so after joint calling a rescue call is told apart only by those two depths, which only rescue records carry,
 * if at all.</p>
 */
public final class Mq0RescueArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String MQ0_RESCUE_LONG_NAME = "mq0-rescue";
    public static final String MQ0_RESCUE_DISTANCE_LONG_NAME = "mq0-rescue-distance";
    public static final String MQ0_RESCUE_MIN_ALT_READS_LONG_NAME = "mq0-rescue-min-alt-reads";
    public static final String MQ0_RESCUE_MIN_ALT_FRACTION_LONG_NAME = "mq0-rescue-min-alt-fraction";
    public static final String MQ0_RESCUE_MIN_CALL_QUAL_LONG_NAME = "mq0-rescue-min-call-qual";
    public static final String MQ0_RESCUE_MAX_EVENT_SETS_LONG_NAME = "mq0-rescue-max-event-sets";
    public static final String MQ0_RESCUE_MODE_LONG_NAME = "mq0-rescue-mode";

    /** The arguments other than --mq0-rescue that only take effect with it. */
    static final List<String> RESCUE_ONLY_ARGUMENTS = List.of(MQ0_RESCUE_DISTANCE_LONG_NAME, MQ0_RESCUE_MIN_ALT_READS_LONG_NAME,
            MQ0_RESCUE_MIN_ALT_FRACTION_LONG_NAME, MQ0_RESCUE_MIN_CALL_QUAL_LONG_NAME, MQ0_RESCUE_MAX_EVENT_SETS_LONG_NAME,
            MQ0_RESCUE_MODE_LONG_NAME);

    /** The arguments that only take effect in INJECT mode. */
    static final List<String> INJECT_ONLY_ARGUMENTS = List.of(MQ0_RESCUE_MAX_EVENT_SETS_LONG_NAME);

    private static final int DEFAULT_DISTANCE = 10;
    private static final double DEFAULT_MIN_CALL_QUAL = 3.0;

    @Argument(fullName = MQ0_RESCUE_LONG_NAME, optional = true,
            doc = "Keep MAPQ-0 reads aside and, where they carry alt evidence away from the regular calls, genotype the region " +
                  "again with all reads; add the non-reference calls with QUAL at least --" + MQ0_RESCUE_MIN_CALL_QUAL_LONG_NAME + " (" +
                  DEFAULT_MIN_CALL_QUAL + ") whose alleles MAPQ-0 reads carry to the output, with FILTER " + GATKVCFConstants.MQ0_RESCUE_FILTER_NAME +
                  " and INFO " + GATKVCFConstants.MQ0_RESCUE_KEY + ". A rescue call is left out if it has a '*' allele, if it overlaps a " +
                  "regular call, an earlier rescue call, a regular record starting before it or a force-called --alleles record, or if a " +
                  "regular call lies within --" + MQ0_RESCUE_DISTANCE_LONG_NAME + " (" + DEFAULT_DISTANCE + ") bp of it. No effect if the " +
                  "read filters let MAPQ-0 reads through. Requires --dragen-mode (or --dragen-378-concordance-mode): the rescue is " +
                  "validated only in DRAGEN mode. Not supported by RampedHaplotypeCaller")
    public boolean enabled = false;

    @Hidden
    @Argument(fullName = MQ0_RESCUE_DISTANCE_LONG_NAME, optional = true, minValue = 0, maxValue = 1000,
            doc = "MAPQ-0 rescue: a rescue record is not written if a regular call lies within this many bp of it, and positions " +
                  "this close to a regular call in the same region do not make the region be genotyped again")
    public int distance = DEFAULT_DISTANCE;

    @Hidden
    @Argument(fullName = MQ0_RESCUE_MIN_ALT_READS_LONG_NAME, optional = true, minValue = 1,
            doc = "MAPQ-0 rescue: genotype a region again only if at some position at least this many MAPQ-0 reads (reads, not " +
                  "fragments) carry the same alt base or indel of A, C, G and T. In INJECT mode those alleles are the ones inserted; in " +
                  "ASSEMBLE mode this only decides whether the region is assembled again. In either mode a rescue record needs only one " +
                  "MAPQ-0 read carrying each allele it calls")
    public int minAltReads = 2;

    @Hidden
    @Argument(fullName = MQ0_RESCUE_MIN_ALT_FRACTION_LONG_NAME, optional = true, minValue = 0, maxValue = 1,
            doc = "MAPQ-0 rescue: genotype a region again only if at some position the MAPQ-0 reads carrying the same allele are at " +
                  "least this fraction of the MAPQ-0 reads with a base of sufficient quality, or a deletion, there. In INJECT mode those " +
                  "alleles are the ones inserted; in ASSEMBLE mode this only decides whether the region is assembled again")
    public double minAltFraction = 0.1;

    @Hidden
    @Argument(fullName = MQ0_RESCUE_MIN_CALL_QUAL_LONG_NAME, optional = true, minValue = 0,
            doc = "MAPQ-0 rescue: a record, regular or rescue, counts as a call only if a genotype is called non-reference and " +
                  "its QUAL is at least this; only regular calls keep rescue records away, and only rescue calls are written. " +
                  "Without reference confidence output, rescue calls below --standard-min-confidence-threshold-for-calling are not made at all")
    public double minCallQual = DEFAULT_MIN_CALL_QUAL;

    @Hidden
    @Argument(fullName = MQ0_RESCUE_MAX_EVENT_SETS_LONG_NAME, optional = true, minValue = 0,
            doc = "MAPQ-0 rescue, INJECT mode: also make one haplotype for each of up to this many sets of two or more of the alleles that at least --" +
                  MQ0_RESCUE_MIN_ALT_READS_LONG_NAME + " MAPQ-0 reads each carry all of, so that phased alleles are not split across haplotypes")
    public int maxEventSets = 8;

    /** How the haplotypes of a rescue call are made. */
    public enum Mode {
        /**
         * No assembly: the haplotypes are the reference haplotype, the reference haplotype with each allele the MAPQ-0
         * reads carry at the trigger positions inserted into it (alleles of A, C, G and T that at least
         * --mq0-rescue-min-alt-reads of them carry, at most --max-num-haplotypes-in-population of them), and one with
         * each set of those alleles the same reads carry together. Pileup detection (--pileup-detection), if on, adds
         * its haplotypes as in regular calling. No assembly graph is built.
         */
        INJECT,
        /**
         * A region is genotyped again under the same condition as in INJECT mode (some trigger position has an allele
         * INJECT would insert), so both modes genotype the same regions again; --mq0-rescue-min-alt-reads and
         * --mq0-rescue-min-alt-fraction do no more than decide that. The region is then assembled again with all
         * reads, as the regular calling path assembles it, and the assembled haplotypes are genotyped. Every call
         * starting in the region span, at a trigger position or not, is kept, dropped and merged by the same rules as
         * INJECT's (see the class documentation), so it is kept if at least one MAPQ-0 read carries each allele it
         * calls. So ASSEMBLE can call alleles at other positions, and combinations, that INJECT does not make,
         * including where the regular calls found too little evidence: since rescue genotyping scores reads with MAPQ
         * above 0 at full base quality, such a call can rest mostly on those reads. Pileup detection
         * (--pileup-detection) works from the regular reads' pileups in both modes. ASSEMBLE is not a superset of
         * INJECT: alleles at a low fraction, which INJECT forces into haplotypes, may not be assembled. A compound
         * allele it calls (an indel and substitutions together, say CT to G) counts as carried by MAPQ-0 reads if they
         * carry its indel and substitutions, however the aligner placed the gap. Assembly makes it costlier than
         * INJECT per region genotyped again.
         */
        ASSEMBLE
    }

    @Advanced
    @Argument(fullName = MQ0_RESCUE_MODE_LONG_NAME, optional = true,
            doc = "MAPQ-0 rescue: how the haplotypes of a region genotyped again are made: by inserting the MAPQ-0 reads' alleles " +
                  "into the reference haplotype (INJECT, the default, cheaper) or by assembling the region again (ASSEMBLE, more " +
                  "costly, finds alleles INJECT does not but may miss some it finds). Has no effect without --" + MQ0_RESCUE_LONG_NAME)
    public Mode mode = Mode.INJECT;

    /**
     * @param parser the parser that parsed these arguments
     * @return the long names of the arguments of {@link #RESCUE_ONLY_ARGUMENTS} given on the command line, in that order
     */
    static List<String> givenRescueOnlyArguments(final CommandLineArgumentParser parser) {
        return givenArguments(parser, RESCUE_ONLY_ARGUMENTS);
    }

    /**
     * @param parser the parser that parsed these arguments
     * @return the long names of the INJECT-only arguments given on the command line
     */
    static List<String> givenInjectOnlyArguments(final CommandLineArgumentParser parser) {
        return givenArguments(parser, INJECT_ONLY_ARGUMENTS);
    }

    private static List<String> givenArguments(final CommandLineArgumentParser parser, final List<String> names) {
        return names.stream().filter(name -> {
            final NamedArgumentDefinition definition = parser.getNamedArgumentDefinitionByAlias(name);
            return definition != null && definition.getHasBeenSet();
        }).collect(Collectors.toList());
    }
}
