package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

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
 * against the reference haplotype with those alleles inserted into it. That genotyping uses likelihood engines that do
 * not cap base qualities at the read's MAPQ, for all of its reads, as DRAGEN mode does for every call, even if the cap
 * is turned back on. Records from it that call a non-reference genotype, with QUAL at least {@code
 * --mq0-rescue-min-call-qual}, whose called alleles MAPQ-0 reads carry, are merged into the output with the FILTER
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

    /** The arguments other than --mq0-rescue that only take effect with it. */
    static final List<String> RESCUE_ONLY_ARGUMENTS = List.of(MQ0_RESCUE_DISTANCE_LONG_NAME, MQ0_RESCUE_MIN_ALT_READS_LONG_NAME,
            MQ0_RESCUE_MIN_ALT_FRACTION_LONG_NAME, MQ0_RESCUE_MIN_CALL_QUAL_LONG_NAME, MQ0_RESCUE_MAX_EVENT_SETS_LONG_NAME);

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
            doc = "MAPQ-0 rescue: genotype a region again only at positions where at least this many MAPQ-0 reads (reads, not " +
                  "fragments) carry the same alt base or indel, and genotype only those alleles")
    public int minAltReads = 2;

    @Hidden
    @Argument(fullName = MQ0_RESCUE_MIN_ALT_FRACTION_LONG_NAME, optional = true, minValue = 0, maxValue = 1,
            doc = "MAPQ-0 rescue: genotype a region again only at positions where the MAPQ-0 reads carrying the same allele are at " +
                  "least this fraction of the MAPQ-0 reads with a base of sufficient quality, or a deletion, there")
    public double minAltFraction = 0.1;

    @Hidden
    @Argument(fullName = MQ0_RESCUE_MIN_CALL_QUAL_LONG_NAME, optional = true, minValue = 0,
            doc = "MAPQ-0 rescue: a record, regular or rescue, counts as a call only if a genotype is called non-reference and " +
                  "its QUAL is at least this; only regular calls keep rescue records away, and only rescue calls are written. " +
                  "Without reference confidence output, rescue calls below --standard-min-confidence-threshold-for-calling are not made at all")
    public double minCallQual = DEFAULT_MIN_CALL_QUAL;

    @Hidden
    @Argument(fullName = MQ0_RESCUE_MAX_EVENT_SETS_LONG_NAME, optional = true, minValue = 0,
            doc = "MAPQ-0 rescue: also make one haplotype for each of up to this many sets of two or more of the alleles that at least --" +
                  MQ0_RESCUE_MIN_ALT_READS_LONG_NAME + " MAPQ-0 reads each carry all of, so that phased alleles are not split across haplotypes")
    public int maxEventSets = 8;

    /**
     * @param parser the parser that parsed these arguments
     * @return the long names of the arguments of {@link #RESCUE_ONLY_ARGUMENTS} given on the command line, in that order
     */
    static List<String> givenRescueOnlyArguments(final CommandLineArgumentParser parser) {
        return RESCUE_ONLY_ARGUMENTS.stream().filter(name -> {
            final NamedArgumentDefinition definition = parser.getNamedArgumentDefinitionByAlias(name);
            return definition != null && definition.getHasBeenSet();
        }).collect(Collectors.toList());
    }
}
