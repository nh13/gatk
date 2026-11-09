package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.genotyper.PloidyModel;
import org.broadinstitute.hellbender.tools.walkers.variantutils.PosteriorProbabilitiesUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Code for estimating the reference confidence
 *
 * This code can estimate the probability that the data for a single sample is consistent with a
 * well-determined REF/REF diploid genotype.
 *
 */
public class ReferenceConfidenceModel {

    // Read attributes holding the original soft-clip start and end, set only when soft clipping is reverted
    public static final String ORIGINAL_SOFTCLIP_START_TAG = "os";
    public static final String ORIGINAL_SOFTCLIP_END_TAG = "oe";

    private final int indelInformativeDepthIndelSize;
    private final int numRefSamplesForPrior;
    private final byte refModelDeletionQuality;
    private final boolean useSoftClippedBases;
    private final boolean flowBasedModel;

    @VisibleForTesting
    protected static final String NON_REF_ALLELE_DESCRIPTION = "Represents any possible alternative allele not already represented at this location by REF and ALT";

    private final PosteriorProbabilitiesUtils.PosteriorProbabilitiesOptions options;

    /**
     * Surrogate quality score for no base calls.
     * <p>
     * This is the quality assigned to deletion (so without its own base-call quality) pile-up elements,
     * when assessing the confidence on the hom-ref call at that site.
     * </p>
     */
    public static final byte REF_MODEL_DELETION_QUAL = 30;

    /**
     * Base calls with quality threshold lower than this number won't be considered when assessing the
     * confidence on the hom-ref call.
     */
    protected static final byte BASE_QUAL_THRESHOLD = 6;

    /**
     * Only base calls with quality strictly greater than this constant,
     * will be considered high quality if they are part of a soft-clip.
     */
    private static final byte HQ_BASE_QUALITY_SOFTCLIP_THRESHOLD = 28;

    //TODO change this: https://github.com/broadinstitute/gsa-unstable/issues/1108
    protected static final int MAX_N_INDEL_INFORMATIVE_READS = 40; // more than this is overkill because GQs are capped at 99 anyway

    private static final int INITIAL_INDEL_LK_CACHE_PLOIDY_CAPACITY = 20;
    private static GenotypeLikelihoods[][] indelPLCache = new GenotypeLikelihoods[INITIAL_INDEL_LK_CACHE_PLOIDY_CAPACITY + 1][];

    /**
     * Indel error rate for the indel model used to assess the confidence on the hom-ref call.
     */
    private static final double INDEL_ERROR_RATE = -4.5; // 10^-4.5 indel errors per bp

    /**
     * Phred scaled qual value that corresponds to the {@link #INDEL_ERROR_RATE indel error rate}.
     */
    private static final byte INDEL_QUAL = (byte) Math.round(INDEL_ERROR_RATE * -10.0);

    /**
     * No indel likelihood (ref allele) used in the indel model to assess the confidence on the hom-ref call.
     */
    private static final double NO_INDEL_LIKELIHOOD = QualityUtils.qualToProbLog10(INDEL_QUAL);

    /**
     * Indel likelihood (alt. allele) used in the indel model to assess the confidence on the hom-ref call.
     */
    private static final double INDEL_LIKELIHOOD = QualityUtils.qualToErrorProbLog10(INDEL_QUAL);
    private static final int IDX_HOM_REF = 0;

    /**
     * Per-element genotype likelihood increments by ploidy: entry [ploidy][alt ? 1 : 0][qual & 0xff] holds, in
     * genotype order, the amounts a pileup element of that base quality adds to the likelihoods at read weight 1,
     * for every quality the quality tables define (up to {@link QualityUtils#MAX_QUAL}).
     * Each element's increments depend only on whether it is alt and on its quality, so they are computed once per
     * ploidy with the same expressions an element would evaluate, and adding them in element order gives the same
     * doubles.
     */
    private static final ConcurrentHashMap<Integer, double[][][]> LIKELIHOOD_INCREMENTS_BY_PLOIDY = new ConcurrentHashMap<>();

    @VisibleForTesting
    static double[][][] likelihoodIncrements(final int ploidy) {
        return LIKELIHOOD_INCREMENTS_BY_PLOIDY.computeIfAbsent(ploidy, ReferenceConfidenceModel::computeLikelihoodIncrements);
    }

    private static double[][][] computeLikelihoodIncrements(final int ploidy) {
        final int likelihoodCount = ploidy + 1;
        final double log10Ploidy = Math.log10(ploidy);
        final double[][][] increments = new double[2][QualityUtils.MAX_QUAL + 1][likelihoodCount];
        for (int alt = 0; alt < 2; alt++) {
            for (int qualIndex = 0; qualIndex <= QualityUtils.MAX_QUAL; qualIndex++) {
                final byte qual = (byte) qualIndex;
                final double referenceLikelihood;
                final double nonRefLikelihood;
                if (alt == 1) {
                    nonRefLikelihood = QualityUtils.qualToProbLog10(qual);
                    referenceLikelihood = QualityUtils.qualToErrorProbLog10(qual) + MathUtils.LOG10_ONE_THIRD;
                } else {
                    referenceLikelihood = QualityUtils.qualToProbLog10(qual);
                    nonRefLikelihood = QualityUtils.qualToErrorProbLog10(qual) + MathUtils.LOG10_ONE_THIRD;
                }
                final double[] entry = increments[alt][qualIndex];
                // Homozygous likelihoods don't need the logSum trick.
                entry[0] = referenceLikelihood + log10Ploidy;
                entry[likelihoodCount - 1] = nonRefLikelihood + log10Ploidy;
                // Heterozygous likelihoods need the logSum trick:
                for (int i = 1, j = likelihoodCount - 2; i < likelihoodCount - 1; i++, j--) {
                    entry[i] = MathUtils.approximateLog10SumLog10(referenceLikelihood + Math.log10(j), nonRefLikelihood + Math.log10(i));
                }
            }
        }
        return increments;
    }

    /**
     * Options related to posterior probability calcs
     */
    private static final boolean useInputSamplesAlleleCounts = false;  //by definition ref-conf will be single-sample; inputs should get ignored but let's be explicit
    private static final boolean useMLEAC = true;
    private static final boolean ignoreInputSamplesForMissingVariants = true;
    private static final boolean useFlatPriorsForIndels = false;


    /**
     * Create a new ReferenceConfidenceModel
     *
     * @param samples the list of all samples we'll be considering with this model
     * @param header the SAMFileHeader describing the read information (used for debugging)
     * @param indelInformativeDepthIndelSize the max size of indels to consider when calculating indel informative depths
     * @param refDelQual reference model deletion quality (30 - Ilmn, 20 - JB)
     * @param useSoftClippedBases should the soft clipped bases be used as evidence against reference
     */
    public ReferenceConfidenceModel(final SampleList samples,
                                    final SAMFileHeader header,
                                    final int indelInformativeDepthIndelSize,
                                    final int numRefForPrior,
                                    final byte refDelQual,
                                    final boolean useSoftClippedBases,
                                    final boolean flowBasedModel) {
        Utils.nonNull(samples, "samples cannot be null");
        Utils.validateArg( samples.numberOfSamples() > 0, "samples cannot be empty");
        Utils.nonNull(header, "header cannot be empty");
        //TODO: code and comment disagree -- which is right?
        Utils.validateArg( indelInformativeDepthIndelSize >= 0, () -> "indelInformativeDepthIndelSize must be >= 1 but got " + indelInformativeDepthIndelSize);


        this.indelInformativeDepthIndelSize = indelInformativeDepthIndelSize;
        this.numRefSamplesForPrior = numRefForPrior;
        this.options = new PosteriorProbabilitiesUtils.PosteriorProbabilitiesOptions(HomoSapiensConstants.SNP_HETEROZYGOSITY,
                HomoSapiensConstants.INDEL_HETEROZYGOSITY, useInputSamplesAlleleCounts, useMLEAC, ignoreInputSamplesForMissingVariants,
                useFlatPriorsForIndels);
        this.refModelDeletionQuality = refDelQual;
        this.useSoftClippedBases = useSoftClippedBases;
        this.flowBasedModel = flowBasedModel;
    }

    /**
     * Get the VCF header lines to include when emitting reference confidence values via {@link #calculateRefConfidence}.
     * @return a non-null set of VCFHeaderLines
     */
    public Set<VCFHeaderLine> getVCFHeaderLines() {
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>();
        headerLines.add(new VCFSimpleHeaderLine(GATKVCFConstants.SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE_NAME, NON_REF_ALLELE_DESCRIPTION));
        return headerLines;
    }

    public List<VariantContext> calculateRefConfidence(final Haplotype refHaplotype,
                                                       final Collection<Haplotype> calledHaplotypes,
                                                       final SimpleInterval paddedReferenceLoc,
                                                       final AssemblyRegion activeRegion,
                                                       final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                       final PloidyModel ploidyModel,
                                                       final List<VariantContext> variantCalls) {
        return calculateRefConfidence(refHaplotype, calledHaplotypes, paddedReferenceLoc, activeRegion, readLikelihoods,
                ploidyModel, variantCalls, false, Collections.emptyList());
    }

    /**
     * Calculate the reference confidence for a single sample given the its read data
     *
     * Returns a list of variant contexts, one for each position in the {@code activeRegion.getLoc()}, each containing
     * detailed information about the certainty that the sample is hom-ref for each base in the region.
     *
     *
     *
     * @param refHaplotype the reference haplotype, used to get the reference bases across activeRegion.getLoc()
     * @param calledHaplotypes a list of haplotypes that segregate in this region, for realignment of the reads in the
     *                         readLikelihoods, corresponding to each reads best haplotype.  Must contain the refHaplotype.
     * @param paddedReferenceLoc the location of refHaplotype (which might be larger than activeRegion.getLoc())
     * @param activeRegion the active region we want to get the reference confidence over
     * @param readLikelihoods a map from a single sample to its PerReadAlleleLikelihoodMap for each haplotype in calledHaplotypes
     * @param ploidyModel indicate the ploidy of each sample in {@code stratifiedReadMap}.
     * @param variantCalls calls made in this region.  The return result will contain any variant call in this list in the
     *                     correct order by genomic position, and any variant in this list will stop us emitting a ref confidence
     *                     under any position it covers (for snps and insertions that is 1 bp, but for deletions its the entire ref span)
     * @return an ordered list of variant contexts that spans activeRegion.getLoc() and includes both reference confidence
     *         contexts as well as calls from variantCalls if any were provided
     */
    public List<VariantContext> calculateRefConfidence(final Haplotype refHaplotype,
                                                       final Collection<Haplotype> calledHaplotypes,
                                                       final SimpleInterval paddedReferenceLoc,
                                                       final AssemblyRegion activeRegion,
                                                       final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                       final PloidyModel ploidyModel,
                                                       final List<VariantContext> variantCalls,
                                                       final boolean applyPriors,
                                                       final List<VariantContext> VCpriors) {
        Utils.nonNull(refHaplotype, "refHaplotype cannot be null");
        Utils.nonNull(calledHaplotypes, "calledHaplotypes cannot be null");
        Utils.validateArg(calledHaplotypes.contains(refHaplotype), "calledHaplotypes must contain the refHaplotype");
        Utils.nonNull(paddedReferenceLoc, "paddedReferenceLoc cannot be null");
        Utils.nonNull(activeRegion, "activeRegion cannot be null");
        Utils.nonNull(readLikelihoods, "readLikelihoods cannot be null");
        Utils.validateArg(readLikelihoods.numberOfSamples() == 1, () -> "readLikelihoods must contain exactly one sample but it contained " + readLikelihoods.numberOfSamples());
        Utils.validateArg( refHaplotype.length() == activeRegion.getPaddedSpan().size(), () -> "refHaplotype " + refHaplotype.length() + " and activeRegion location size " + activeRegion.getSpan().size() + " are different");
        Utils.nonNull(ploidyModel, "the ploidy model cannot be null");
        final int ploidy = ploidyModel.samplePloidy(0); // the first sample = the only sample in reference-confidence mode.

        final SimpleInterval refSpan = activeRegion.getSpan();
        final byte[] ref = refHaplotype.getBases();
        final List<ReferenceConfidenceResult> siteResults = calculateSiteResults(readLikelihoods, activeRegion, ref, ploidy);
        final List<VariantContext> results = new ArrayList<>(refSpan.size());
        final String sampleName = readLikelihoods.getSample(0);

        final int globalRefOffset = refSpan.getStart() - activeRegion.getPaddedSpan().getStart();
        final int refSpanSize = refSpan.size();
        for (int i = 0; i < refSpanSize; i++) {
            final int position = refSpan.getStart() + i;
            final Locatable curPos = new SimpleInterval(refSpan.getContig(), position, position);

            final VariantContext overlappingSite = GATKVariantContextUtils.getOverlappingVariantContext(curPos, variantCalls);
            final List<VariantContext> currentPriors = VCpriors.isEmpty() ? Collections.emptyList() : getMatchingPriors(curPos, overlappingSite, VCpriors);
            if (overlappingSite != null && overlappingSite.getStart() == curPos.getStart()) {
                if (applyPriors) {
                    results.add(PosteriorProbabilitiesUtils.calculatePosteriorProbs(overlappingSite, currentPriors,
                            numRefSamplesForPrior, options));
                } else {
                    results.add(overlappingSite);
                }
            } else {
                // otherwise emit a reference confidence variant context
                results.add(makeReferenceConfidenceVariantContext(ploidy, ref[i + globalRefOffset], sampleName, siteResults.get(i), curPos, applyPriors, currentPriors));
            }
        }

        return results;
    }

    /**
     * Computes the reference confidence result of every position in the region's active span, in order.
     *
     * Reads are visited one at a time in coordinate order, and each read adds its evidence to every span position it
     * aligns to, so a position accumulates its reads in the order a pileup at that position would list them and the
     * likelihood sums are the same as a per-position pileup would produce.
     *
     * @param readLikelihoods the single-sample read likelihoods holding the reads to use
     * @param activeRegion the region whose active span is evaluated
     * @param ref the reference bases of the region's padded span
     * @param ploidy the sample ploidy
     * @return one result per position of {@code activeRegion.getSpan()}
     */
    protected List<ReferenceConfidenceResult> calculateSiteResults(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                                   final AssemblyRegion activeRegion,
                                                                   final byte[] ref,
                                                                   final int ploidy) {
        final SimpleInterval span = activeRegion.getSpan();
        final int spanStart = span.getStart();
        final int spanSize = span.size();
        final int globalRefOffset = spanStart - activeRegion.getPaddedSpan().getStart();
        final int likelihoodCount = ploidy + 1;
        final double log10Ploidy = Math.log10(ploidy);
        final double[][][] increments = likelihoodIncrements(ploidy);

        final RefVsAnyResult[] sites = new RefVsAnyResult[spanSize];
        for (int i = 0; i < spanSize; i++) {
            sites[i] = new RefVsAnyResult(likelihoodCount);
        }
        final int[] readCounts = new int[spanSize];
        final int[] indelInformativeReads = new int[spanSize];

        final List<GATKRead> reads = new ArrayList<>(readLikelihoods.sampleEvidence(0));
        reads.sort(new ReadCoordinateComparator(activeRegion.getHeader()));
        for (final GATKRead read : reads) {
            if (read.getEnd() >= spanStart && read.getStart() <= span.getEnd()) {
                sweepRead(read, spanStart, spanSize, ref, globalRefOffset, increments, sites, readCounts, indelInformativeReads);
            }
        }

        final List<ReferenceConfidenceResult> results = new ArrayList<>(spanSize);
        for (int i = 0; i < spanSize; i++) {
            final RefVsAnyResult site = sites[i];
            final double denominator = readCounts[i] * log10Ploidy;
            for (int k = 0; k < likelihoodCount; k++) {
                site.genotypeLikelihoods[k] -= denominator;
            }
            applyIndelRefConfidence(ploidy, indelInformativeReads[i], site);
            results.add(site);
        }
        return results;
    }

    /**
     * Adds one read's evidence to every span position it aligns to.
     *
     * The alignment is walked as the state machine a pileup uses walks it, so the bases, deletions and cigar context
     * seen here are those a pileup element would carry, but a cigar element at a time and over the read's base and
     * quality arrays rather than a step per position. Reference skips are not evidence, nor are bases inside the
     * adaptor of a short fragment. A base contributes to the SNP likelihoods only above the base quality threshold,
     * while every aligned base that is not part of, or immediately before, an indel contributes to the
     * indel-informative count. A read that starts with a deletion, or that ends with one before passing the end of
     * the span, is rejected as the state machine rejects it.
     */
    @VisibleForTesting
    void sweepRead(final GATKRead read, final int spanStart, final int spanSize, final byte[] ref, final int globalRefOffset,
                           final double[][][] increments,
                           final RefVsAnyResult[] sites, final int[] readCounts, final int[] indelInformativeReads) {
        // When soft-clipped bases are not evidence, region finalization reverts every read's soft clips and records
        // where they were in these tags, so every read reaching here carries them.
        final boolean skipOriginalSoftClips = !useSoftClippedBases;
        final int originalSoftStart = skipOriginalSoftClips ? getOriginalSoftStart(read) : 0;
        final int originalSoftEnd = skipOriginalSoftClips ? getOriginalSoftEnd(read) : 0;
        final int adaptorBoundary = ReadUtils.adaptorBoundaryForBaseChecks(read);
        final boolean hasAdaptorBoundary = adaptorBoundary != ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY;
        final boolean reverseStrand = read.isReverseStrand();
        final byte[] bases = read.getBasesNoCopy();
        final byte[] quals = read.getBaseQualitiesNoCopy();
        final int spanEnd = spanStart + spanSize; // exclusive

        // The offset of the last read base consumed and the last reference position stepped onto, as the state
        // machine tracks them: a deletion's positions keep the offset of the base before it.
        int readOffset = -1;
        int position = read.getStart() - 1;
        // Reference-aligned read offset (insertions collapsed, deletions padded) at which the current cigar element starts.
        int alignedOffsetOfElement = 0;
        BitSet indelInformativeBases = null;

        // The cigar's own list, as the state machine walks: for a BAM record whose cigar is stored in its CG tag,
        // numCigarElements() reports the placeholder cigar's length.
        final List<CigarElement> elements = read.getCigarElements();
        final int nElements = elements.size();
        for (int elementIndex = 0; elementIndex < nElements; elementIndex++) {
            final CigarElement element = elements.get(elementIndex);
            final CigarOperator op = element.getOperator();
            final int length = element.getLength();
            if (length == 0 || op == CigarOperator.H || op == CigarOperator.P) {
                continue;
            }
            if (op == CigarOperator.I || op == CigarOperator.S) {
                readOffset += length;
                alignedOffsetOfElement += alignedLength(element);
                continue;
            }
            final boolean isDeletion = op == CigarOperator.D;
            if (isDeletion && readOffset < 0) {
                throw ReadUtils.readStartsWithDeletionError(read);
            }
            final int firstPosition = position + 1;
            if (firstPosition >= spanEnd) {
                return;
            }
            final int lastPosition = position + length;
            if (op != CigarOperator.N) {
                // a deletion's positions all carry the offset of the base before it; a match's step through the read
                final int firstReadOffset = isDeletion ? readOffset : readOffset + 1;
                final int readOffsetStep = isDeletion ? 0 : 1;
                final int stop = Math.min(lastPosition, spanEnd - 1);
                for (int p = Math.max(firstPosition, spanStart); p <= stop; p++) {
                    if (hasAdaptorBoundary && ReadUtils.isPastAdaptorBoundary(reverseStrand, adaptorBoundary, p)) {
                        continue;
                    }
                    final int i = p - spanStart;
                    final int offsetInElement = p - firstPosition;
                    final int baseOffset = firstReadOffset + readOffsetStep * offsetInElement;
                    final byte refBase = ref[i + globalRefOffset];

                    final byte qual = isDeletion ? deletionQuality(read, baseOffset, refBase, true) : quals[baseOffset];
                    final boolean usable = !((qual <= BASE_QUAL_THRESHOLD) && (flowBasedModel || !isDeletion))
                            && !(skipOriginalSoftClips && (originalSoftStart > p || originalSoftEnd < p));
                    if (usable) {
                        readCounts[i]++;
                        final boolean isAlt = isDeletion || bases[baseOffset] != refBase;
                        applyRefVsNonRefLikelihoodAndCount(increments, sites[i], isAlt, qual, 1.0);
                    }

                    // A position's indel-informative count is only used up to MAX_N_INDEL_INFORMATIVE_READS, so a position that
                    // has reached it is skipped. In deep pileups most positions are full before most reads arrive, so most
                    // reads never compute their indel-informative bases; a read that does computes them from the first
                    // position that still counts, and bits from there on do not depend on where the computation starts.
                    if (isDeletion || indelInformativeReads[i] >= MAX_N_INDEL_INFORMATIVE_READS) {
                        continue;
                    }
                    final boolean beforeIndel = offsetInElement == length - 1
                            && (nextOnGenomeOperatorIsDeletion(elements, elementIndex) || nextOperatorIs(elements, elementIndex, CigarOperator.I));
                    if (!beforeIndel) {
                        final int alignedOffset = alignedOffsetOfElement + offsetInElement;
                        if (indelInformativeBases == null) {
                            indelInformativeBases = indelInformativeBases(read, alignedOffset, ref, i + globalRefOffset, indelInformativeDepthIndelSize);
                        }
                        if (indelInformativeBases.get(alignedOffset)) {
                            indelInformativeReads[i]++;
                        }
                    }
                }
                if (!isDeletion) {
                    readOffset += length;
                }
            }
            if (lastPosition >= spanEnd) {
                return;
            }
            position = lastPosition;
            alignedOffsetOfElement += alignedLength(element);
        }
        if (nElements > 0 && elements.get(nElements - 1).getOperator() == CigarOperator.D) {
            throw ReadUtils.readEndsWithDeletionError(read);
        }
    }

    /**
     * @return the number of reference-aligned read offsets a cigar element spans: soft clips and reference-consuming
     *         operators count, insertions and hard clips do not
     */
    @VisibleForTesting
    static int alignedLength(final CigarElement element) {
        final CigarOperator op = element.getOperator();
        return op.consumesReferenceBases() || op == CigarOperator.S ? element.getLength() : 0;
    }

    /**
     * @return whether the cigar element right after the one at {@code elementIndex} exists and has operator {@code op}
     */
    @VisibleForTesting
    static boolean nextOperatorIs(final List<CigarElement> elements, final int elementIndex, final CigarOperator op) {
        return elementIndex + 1 < elements.size() && elements.get(elementIndex + 1).getOperator() == op;
    }

    /**
     * @return whether the next cigar element after the one at {@code elementIndex} that consumes reference bases is a
     *         deletion, looking past clips, insertions and pads
     */
    @VisibleForTesting
    static boolean nextOnGenomeOperatorIsDeletion(final List<CigarElement> elements, final int elementIndex) {
        final int n = elements.size();
        for (int k = elementIndex + 1; k < n; k++) {
            final CigarOperator op = elements.get(k).getOperator();
            if (op == CigarOperator.D) {
                return true;
            } else if (op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X) {
                return false;
            }
        }
        return false;
    }

    /**
     * Builds the reference-confidence variant context for one position: a hom-ref genotype of the given ploidy over the
     * reference base and the non-ref symbolic allele, carrying the site's AD, DP, PL and GQ.
     *
     * @param ploidy the sample ploidy
     * @param refBase the reference base at the position
     * @param sampleName the sample the genotype belongs to
     * @param homRefCalc the finished reference confidence result for the position
     * @param curPos the position
     * @param applyPriors whether to fold the given priors into the genotype's posteriors
     * @param VCpriors the priors at the position, used only when applyPriors is set
     * @return the variant context for the position
     */
    public VariantContext makeReferenceConfidenceVariantContext(final int ploidy,
                                                                final byte refBase,
                                                                final String sampleName,
                                                                final ReferenceConfidenceResult homRefCalc,
                                                                final Locatable curPos,
                                                                final boolean applyPriors,
                                                                final List<VariantContext> VCpriors) {
        // Assume infinite population on a single sample.
        final Allele refAllele = Allele.create(refBase, true);
        final List<Allele> refSiteAlleles = Arrays.asList(refAllele, Allele.NON_REF_ALLELE);
        final VariantContextBuilder vcb = new VariantContextBuilder("HC", curPos.getContig(), curPos.getStart(), curPos.getStart(), refSiteAlleles);
        final GenotypeBuilder gb = new GenotypeBuilder(sampleName, GATKVariantContextUtils.homozygousAlleleList(refAllele, ploidy));
        gb.AD(homRefCalc.getAD());
        gb.DP(homRefCalc.getDP());
        addGenotypeData(homRefCalc, gb);
        if(!applyPriors) {
            return vcb.genotypes(gb.make()).make();
        }
        else {
            return PosteriorProbabilitiesUtils.calculatePosteriorProbs(vcb.genotypes(gb.make()).make(), VCpriors, numRefSamplesForPrior, options);
            //TODO FIXME: after new-qual refactoring, these should be static calls to AF calculator
        }
    }

    /**
     * Combines a site's SNP genotype likelihoods with the indel likelihoods implied by its indel-informative read count
     * into the site's final PLs.
     */
    @VisibleForTesting
    void applyIndelRefConfidence(final int ploidy, final int nIndelInformativeReads, final RefVsAnyResult homRefCalc) {
        final GenotypeLikelihoods snpGLs = GenotypeLikelihoods.fromLog10Likelihoods(homRefCalc.getGenotypeLikelihoodsCappedByHomRefLikelihood());
        final GenotypeLikelihoods indelGLs = getIndelPLs(ploidy,nIndelInformativeReads);

        // now that we have the SNP and indel GLs, we take the one with the least confidence,
        // as this is the most conservative estimate of our certainty that we are hom-ref.
        // For example, if the SNP PLs are 0,10,100 and the indel PLs are 0,100,1000
        // we are very certain that there's no indel here, but the SNP confidence imply that we are
        // far less confident that the ref base is actually the only thing here.  So we take 0,10,100
        // as our GLs for the site.
        final GenotypeLikelihoods leastConfidenceGLs = getGLwithWorstGQ(indelGLs, snpGLs);

        homRefCalc.finalPhredScaledGenotypeLikelihoods = leastConfidenceGLs.getAsPLs();
    }

    public void addGenotypeData(final ReferenceConfidenceResult result, final GenotypeBuilder gb) {
        final int[] pls = ((RefVsAnyResult)result).finalPhredScaledGenotypeLikelihoods;
        gb.PL(pls);
        gb.GQ(GATKVariantContextUtils.calculateGQFromPLs(pls));
    }

    /**
     * Get the GenotypeLikelihoods with the least strong corresponding GQ value
     * @param gl1 first to consider (cannot be null)
     * @param gl2 second to consider (cannot be null)
     * @return gl1 or gl2, whichever has the worst GQ
     */
    @VisibleForTesting
    GenotypeLikelihoods getGLwithWorstGQ(final GenotypeLikelihoods gl1, final GenotypeLikelihoods gl2) {
        if (getGQForHomRef(gl1) > getGQForHomRef(gl2)) {
            return gl1;
        } else {
            return gl2;
        }
    }

    private double getGQForHomRef(final GenotypeLikelihoods gls){
        return GenotypeLikelihoods.getGQLog10FromLikelihoods(IDX_HOM_REF, gls.getAsVector());
    }

    /**
     * Get indel PLs corresponding to seeing N nIndelInformativeReads at this site
     *
     * @param nInformativeReads the number of reads that inform us about being ref without an indel at this site
     * @param ploidy the requested ploidy.
     * @return non-null GenotypeLikelihoods given N
     */
    @VisibleForTesting
    GenotypeLikelihoods getIndelPLs(final int ploidy, final int nInformativeReads) {
        return indelPLCache(ploidy, nInformativeReads > MAX_N_INDEL_INFORMATIVE_READS ? MAX_N_INDEL_INFORMATIVE_READS : nInformativeReads);
    }

    private GenotypeLikelihoods indelPLCache(final int ploidy, final int nInformativeReads) {
        return initializeIndelPLCache(ploidy)[nInformativeReads];
    }

    private GenotypeLikelihoods[] initializeIndelPLCache(final int ploidy) {

        if (indelPLCache.length <= ploidy) {
            indelPLCache = Arrays.copyOf(indelPLCache, ploidy << 1);
        }

        if (indelPLCache[ploidy] != null) {
            return indelPLCache[ploidy];
        }

        final double denominator =  - Math.log10(ploidy);
        final GenotypeLikelihoods[] result = new GenotypeLikelihoods[MAX_N_INDEL_INFORMATIVE_READS + 1];

        //Note: an array of zeros is the right answer for result[0].
        result[0] = GenotypeLikelihoods.fromLog10Likelihoods(new double[ploidy + 1]);
        for( int nInformativeReads = 1; nInformativeReads <= MAX_N_INDEL_INFORMATIVE_READS; nInformativeReads++ ) {
            final double[] PLs = new double[ploidy + 1];
            PLs[0] = nInformativeReads * NO_INDEL_LIKELIHOOD;
            for (int altCount = 1; altCount <= ploidy; altCount++) {
                final double refLikelihoodAccum = NO_INDEL_LIKELIHOOD + Math.log10(ploidy - altCount);
                final double altLikelihoodAccum = INDEL_LIKELIHOOD + Math.log10(altCount);
                PLs[altCount] = nInformativeReads * (MathUtils.approximateLog10SumLog10(refLikelihoodAccum ,altLikelihoodAccum) + denominator);
            }
            result[nInformativeReads] = GenotypeLikelihoods.fromLog10Likelihoods(PLs);
        }
        indelPLCache[ploidy] = result;
        return result;
    }

    /**
     * Calculate the genotype likelihoods for the sample in pileup for being hom-ref contrasted with being ref vs. alt
     *
     * @param ploidy target sample ploidy.
     * @param pileup the read backed pileup containing the data we want to evaluate
     * @param refBase the reference base at this pileup position
     * @param minBaseQual the min base quality for a read in the pileup at the pileup position to be included in the calculation
     * @param hqSoftClips running average data structure (can be null) to collect information about the number of high quality soft clips
     * @return a RefVsAnyResult genotype call.
     */
    public ReferenceConfidenceResult calcGenotypeLikelihoodsOfRefVsAny(final int ploidy,
                                                                       final ReadPileup pileup,
                                                                       final byte refBase,
                                                                       final byte minBaseQual,
                                                                       final MathUtils.RunningAverage hqSoftClips,
                                                                       final boolean readsWereRealigned,
                                                                       final double altReadWeight) {

        final int likelihoodCount = ploidy + 1;
        final double log10Ploidy = Math.log10(ploidy);
        final double[][][] increments = likelihoodIncrements(ploidy);

        final RefVsAnyResult result = new RefVsAnyResult(likelihoodCount);
        int readCount = 0;
        for (final PileupElement p : pileup) {
            //note that reference confidence model is used both in active region detection (readsWereRealigned=false) and
            //after the assembly. Only in the latter case flow based model is used in getDeletionQuality, in
            //active region detection we for now use more sensitive reference confidence model where any deletion is
            //a strong evidence against reference.
            final byte qual = countedQuality(p, refBase, readsWereRealigned);
            if (isDroppedByBaseQuality(p, qual, minBaseQual)) {
                continue;
            }
            if (!useSoftClippedBases && readsWereRealigned){
                int loc = pileup.getLocation().getStart();
                //skip bases that were originally softclipped
                if ((getOriginalSoftStart(p.getRead()) > loc) || (getOriginalSoftEnd(p.getRead()) < loc)){
                    continue;
                }
            }

            readCount++;
            applyPileupElementRefVsNonRefLikelihoodAndCount(refBase, increments, result, p, qual, hqSoftClips, readsWereRealigned, altReadWeight);
        }
        final double denominator = readCount * log10Ploidy;
        for (int i = 0; i < likelihoodCount; i++) {
            result.genotypeLikelihoods[i] -= denominator;
        }
        return result;
    }

    public ReferenceConfidenceResult calcGenotypeLikelihoodsOfRefVsAny(final int ploidy,
                                                                       final ReadPileup pileup,
                                                                       final byte refBase,
                                                                       final byte minBaseQual,
                                                                       final MathUtils.RunningAverage hqSoftClips,
                                                                       final boolean readsWereRealigned) {
        return calcGenotypeLikelihoodsOfRefVsAny(ploidy, pileup, refBase, minBaseQual, hqSoftClips, readsWereRealigned, 1.0);
    }

    /**
     * Whether any element that this class's {@link #calcGenotypeLikelihoodsOfRefVsAny} counts before assembly (with
     * {@code readsWereRealigned} false) is evidence for a non-reference allele, i.e. whether that method would report
     * a non-zero non-reference depth. The two share the element filter and the alt test, so they agree exactly.
     *
     * @param pileup the pileup at the site
     * @param refBase the reference base at the site
     * @param minBaseQual bases with this quality or lower are not counted; outside flow mode deletions always are
     * @return true if a counted element is alt by {@link #isAltBeforeAssembly}
     */
    final boolean hasAltEvidenceBeforeAssembly(final ReadPileup pileup, final byte refBase, final byte minBaseQual) {
        for (final PileupElement p : pileup) {
            if (!isDroppedByBaseQuality(p, countedQuality(p, refBase, false), minBaseQual) && isAltBeforeAssembly(p, refBase)) {
                return true;
            }
        }
        return false;
    }

    /** The quality an element is counted with: its base quality, or for a deletion {@link #getDeletionQuality}. */
    private byte countedQuality(final PileupElement p, final byte refBase, final boolean readsWereRealigned) {
        return p.isDeletion() ? getDeletionQuality(p, refBase, readsWereRealigned) : p.getQual();
    }

    /**
     * Whether the base-quality filter drops an element counted with quality {@code qual}: any element at or below
     * {@code minBaseQual}, except that outside flow mode every deletion is kept.
     */
    private boolean isDroppedByBaseQuality(final PileupElement p, final byte qual, final byte minBaseQual) {
        return qual <= minBaseQual && (flowBasedModel || !p.isDeletion());
    }

    private int getOriginalSoftStart(GATKRead read) {
        if (!read.hasAttribute(ORIGINAL_SOFTCLIP_START_TAG)){
            throw new GATKException("Attempt to read soft clip start that was not saved");
        } else {
            return read.getAttributeAsInteger(ORIGINAL_SOFTCLIP_START_TAG);
        }
    }

    private int getOriginalSoftEnd(GATKRead read) {
        if (!read.hasAttribute(ORIGINAL_SOFTCLIP_END_TAG)){
            throw new GATKException("Attempt to read soft clip end that was not saved");
        } else {
            return read.getAttributeAsInteger(ORIGINAL_SOFTCLIP_END_TAG);
        }
    }


    /**
     *  in flow based read the deletion quality can be found as the quality of the nearby base
     *  assuming that the deletion is hmer deletion
     *  we use this only in the case of using reference confidence model after the assembly/realignment.
     *  when reference confidence model is used for active region detection we use standard reference
     *  confidence where any deletion has a constant quality
     *
     */
    private byte getDeletionQuality(PileupElement p, byte refBase, final boolean notInIsActive) {
        return deletionQuality(p.getRead(), p.getOffset(), refBase, notInIsActive);
    }

    /**
     * The quality of a deletion at a position whose read offset is that of the base before the deletion. Under the
     * flow-based model, outside active region detection, a deletion followed by a base equal to {@code refBase} is taken
     * as a left-aligned homopolymer deletion and gets that base's quality; every other deletion gets the constant
     * reference model deletion quality.
     */
    @VisibleForTesting
    byte deletionQuality(final GATKRead read, final int offset, final byte refBase, final boolean notInIsActive) {
        if (flowBasedModel && notInIsActive){
            if (read.getBase(offset + 1) == refBase){ // if hmer indel - assume that deletion is left aligned
                return read.getBaseQuality(offset + 1);
            }
        }
        return refModelDeletionQuality;
    }

    private void applyPileupElementRefVsNonRefLikelihoodAndCount(final byte refBase, final double[][][] increments, final RefVsAnyResult result, final PileupElement element, final byte qual, final MathUtils.RunningAverage hqSoftClips, final boolean readsWereRealigned, final double altReadWeight) {
        final boolean isAlt = readsWereRealigned ? isAltAfterAssembly(element, refBase) : isAltBeforeAssembly(element, refBase);
        applyRefVsNonRefLikelihoodAndCount(increments, result, isAlt, qual, altReadWeight);
        if (isAlt && hqSoftClips != null && element.isNextToSoftClip()) {
            hqSoftClips.add(AlignmentUtils.countHighQualitySoftClips(element.getRead(), HQ_BASE_QUALITY_SOFTCLIP_THRESHOLD));
        }
    }

    /**
     * Adds one element's evidence to a site: counts it in the ref or non-ref depth and adds its genotype likelihood
     * increments, weighted by {@code altReadWeight} if it is an alt element and by 1 otherwise.
     */
    @VisibleForTesting
    static void applyRefVsNonRefLikelihoodAndCount(final double[][][] increments, final RefVsAnyResult result, final boolean isAlt, final byte qual, final double altReadWeight) {
        if (isAlt) {
            result.nonRefDepth++;
        } else {
            result.refDepth++;
        }
        final double readWeight = isAlt ? altReadWeight : 1.0;
        final double[] entry = increments[isAlt ? 1 : 0][qual & 0xff];
        final double[] genotypeLikelihoods = result.genotypeLikelihoods;
        for (int k = 0; k < entry.length; k++) {
            genotypeLikelihoods[k] += readWeight * entry[k];
        }
    }

    public static boolean isAltBeforeAssembly(final PileupElement element, final byte refBase){
        return element.getBase() != refBase || element.isDeletion() || element.isBeforeDeletionStart()
                || element.isAfterDeletionEnd() || element.isBeforeInsertion() || element.isAfterInsertion() || element.isNextToSoftClip();
    }

    protected static boolean isAltAfterAssembly(final PileupElement element, final byte refBase){
        return element.getBase() != refBase || element.isDeletion(); //we shouldn't have soft clips after assembly
    }


    /**
     * Note that we don't have to match alleles because the PosteriorProbabilitesUtils will take care of that
     * @param curPos position of interest for genotyping
     * @param call (may be null)
     * @param priorList priors within the current ActiveRegion
     * @return prior VCs representing the same variant position as call
     */
    private List<VariantContext> getMatchingPriors(final Locatable curPos, final VariantContext call, final List<VariantContext> priorList) {
        final int position = call != null ? call.getStart() : curPos.getStart();
        final List<VariantContext> matchedPriors = new ArrayList<>(priorList.size());
        // NOTE: a for loop is used here because this method ends up being called per-pileup, per-read and using a loop instead of streaming saves runtime
        final int priorsListSize = priorList.size();
        for (int i = 0; i < priorsListSize; i++) {
            if (position == priorList.get(i).getStart()) {
                matchedPriors.add(priorList.get(i));
            }
        }
        return matchedPriors;
    }

    /**
     * Compute the sum of mismatching base qualities for readBases aligned to refBases at readStart / refStart
     * assuming no insertions or deletions in the read w.r.t. the reference
     *
     * @param readBases non-null bases of the read
     * @param readQuals non-null quals of the read
     * @param readStart the starting position of the read (i.e., that aligns it to a position in the reference)
     * @param refBases the reference bases
     * @param refStart the offset into refBases that aligns to the readStart position in readBases
     * @return an array containing the sum of quality scores for readBases that mismatch following this base and their corresponding ref base for each read base in readBases
     */
    private static int[] calculateBaselineMMQualities(final byte[] readBases,
                                final byte[] readQuals,
                                final int readStart,
                                final byte[] refBases,
                                final int refStart) {
        final int n = Math.min(readBases.length - readStart, refBases.length - refStart);
        int[] results = new int[n];
        int sum = 0;

        // Note that we start this loop at the end based on the principle that in order to calculate the number of mismatches remaining
        // between the read and the reference after the nth base, one can simply first calculate the remaining mismatches for the n + 1th
        // base first and so on.
        for ( int i = n - 1; i >= 0; i-- ) {
            final byte readBase = readBases[readStart + i];
            final byte refBase  = refBases[refStart + i];
            if (isMismatchAndNotAnAlignmentGap(readBase, refBase)) {
                sum += readQuals[readStart + i];
            }
            results[i] = sum;
        }

        return results;
    }

    /**
     * Compute whether a read is informative to eliminate an indel of size <= maxIndelSize segregating at readStart/refStart
     *
     * For each base this method determines if there are any plausible indels of size <= maxIndelSize that start at that
     * base. The method returns true if no indels were found that align as well or better than the rest of this read
     * compared to the reference.
     *
     * The result covers every read offset from readStart to the end of the read: bit i of the returned set is true when
     * the read has no plausible indel at reference-aligned offset i. Bits at offsets before readStart are not meaningful.
     * The bits at and after readStart do not depend on which offset the computation was anchored at, so one call per
     * read serves every later position of that read. This holds for base qualities below 128; the mismatch sums that
     * decide it are over signed bytes.
     *
     * Positions <= maxIndelSize from the end of the provided read/ref are always false.
     *
     * @param read the read
     * @param readStart the 0-based index with respect to @{param}refBases where the read starts (this is the "IGV View" offset for the read)
     * @param refBases the reference bases
     * @param refStart the 0-based offset into refBases that aligns to the readStart position in readBases
     * @param maxIndelSize the max indel size to consider for the read to be informative
     * @return the set of reference-aligned read offsets at which the read rules out an indel of size <= maxIndelSize
     */
    @VisibleForTesting
    static BitSet indelInformativeBases(final GATKRead read,
                                        final int readStart,
                                        final byte[] refBases,
                                        final int refStart,
                                        final int maxIndelSize) {
        Utils.validate(readStart >= 0, "readStart must >= 0");
        Utils.validate(refStart >= 0, "refStart must >= 0");
        BitSet informativeBases = new BitSet(read.getLength());

        // Check that we aren't so close to the end of the end of the read that we don't have to compute anything more
        if ( !(read.getLength() - readStart < maxIndelSize) && !(refBases.length - refStart < maxIndelSize) ) {
            //TODO this should be removed, see https://github.com/broadinstitute/gatk/issues/5646 to track its progress
            final int secondaryReadBreakPosition = read.getLength() - maxIndelSize;

            // We are safe to use the faster no-copy versions of getBases and getBaseQualities here,
            // since we're not modifying the returned arrays in any way. This makes a small difference
            // in the HaplotypeCaller profile, since this method is a major hotspot.
            final Pair<byte[], byte[]> readBasesAndBaseQualities = AlignmentUtils.getBasesAndBaseQualitiesAlignedOneToOne(read);  //calls getBasesNoCopy if CIGAR is all match
            final byte[] readBases = readBasesAndBaseQualities.getLeft();
            final byte[] readQualities = readBasesAndBaseQualities.getRight();

            // Need to check for closeness to the end of the read again as the array size may be different than read.Len() due to deletions in the cigar
            if (readBases.length - readStart > maxIndelSize) {

                // Compute where the end of marking would have been given the above two break conditions so we can stop marking there
                final int lastReadBaseToMarkAsIndelRelevant;
                final boolean referenceWasShorter;
                if (readBases.length < refBases.length - refStart + readStart + 1) {
                    // If the read ends first, then we don't mark the last maxIndelSize bases from it as relevant
                    lastReadBaseToMarkAsIndelRelevant = readBases.length - maxIndelSize;
                    referenceWasShorter = false;
                } else {
                    // If the reference ends first, then we don't mark the last maxIndelSize bases from it as relevant
                    lastReadBaseToMarkAsIndelRelevant = refBases.length - refStart + readStart - maxIndelSize + 1;
                    referenceWasShorter = true;
                }


                // Compute the absolute baseline sum against which to test
                final int[] baselineMisMatchSums = calculateBaselineMMQualities(readBases, readQualities, readStart, refBases, refStart);

                // consider each indel size up to max in term, checking if an indel that deletes either the ref bases (deletion)
                // or read bases (insertion) would fit as well as the origin baseline sum of mismatching quality scores. These scores
                // are computed starting from the last base in the read/reference that would be offset by the indel and compared against
                // the mismatch cost for the same base of the reference. Once the sum of mismatch qualities counting from the back for
                // one indel size exceeds the global indel mismatch cost, the code stops as it will never find a better mismatch value.
                for (int indelSize = 1; indelSize <= maxIndelSize; indelSize++) {
                    // Computing mismatches corresponding to a deletion
                    traverseEndOfReadForIndelMismatches(informativeBases,
                            readStart,
                            readBases,
                            readQualities,
                            lastReadBaseToMarkAsIndelRelevant,
                            secondaryReadBreakPosition,
                            refStart,
                            refBases,
                            baselineMisMatchSums,
                            indelSize,
                            false);

                    // Computing mismatches corresponding to an insertion
                    traverseEndOfReadForIndelMismatches(informativeBases,
                            readStart,
                            readBases,
                            readQualities,
                            lastReadBaseToMarkAsIndelRelevant,
                            secondaryReadBreakPosition,
                            refStart,
                            refBases,
                            baselineMisMatchSums,
                            indelSize,
                            true);
                }


                // Flip the bases at the front of the read (the ones not within maxIndelSize of the end as those are never informative)
                // These must be flipped because thus far we have marked reads for which there were plausible indels with a true value in
                // the bitset. This method returns false for cases where we have discovered plausible indels so we must flip them. This
                // is done in part to preserve a sensible default behavior for bases not considered by this approach.
                if ( lastReadBaseToMarkAsIndelRelevant <= secondaryReadBreakPosition) {
                    informativeBases.flip(0, lastReadBaseToMarkAsIndelRelevant);
                    // Resolve the fact that the old approach would always mark the last base examined as being indel uninformative when the reference
                    // ends first despite it corresponding to a comparison of zero bases against the read
                    if (referenceWasShorter) {
                        informativeBases.set(lastReadBaseToMarkAsIndelRelevant - 1, false);
                    }
                } else {
                    informativeBases.flip(0, secondaryReadBreakPosition + 1);
                }

            }
        }
        return informativeBases;
    }

    /**
     * Helper method responsible for read-end traversal. This method will handle both insertions and deletions,
     * indicated by setting the insertion parameter.
     *
     * Given the array of sums baselineMMSums, this method will start from the back of the read and reference and sum the
     * quality score of all mismatching bases to the reference. If the score is equal to or lower than the baseline sum and
     * the base being examined is before lastReadBaseToMarkAsIndelRelevant, then this method will store a true into the informativeBases
     * Bitset for that particular read. If at any point the sum for a given size insertion/deletion exceeds the global cost
     * of all aligned mismatches to the reference with no indels added (the first position in baselineMMSums) the process will
     * end prematurely so as to avoid comparing any additional bases beyond what is necessary.
     *
     * It is expected that only the bases between readStart and lastReadBaseToMarkAsIndelRelevant in the bitset will be set to true
     * by this method if they are ambiguous about an indel of the given size. We then flip these values later in the process because
     * an ambiguous indel positions in the read actually return false in indelInformativeBases.
     *
     * NOTE: This method examines overhanging bases to the reference/read if they do not end at the same position.
     *       (eg. if the reference ends 20 bases after the read does and you are looking at a deletion of size 5, the first
     *       base compared will be the last base of the read and the 15th from last base on the reference)
     *
     * @param informativeBases ReadBases indexed bitset into which to store the results
     * @param readStart Offset of first comparison base into the read
     * @param readBases Read bases aligned to be indexed by reference base
     * @param readQuals Read qualities aligned to be indexed by reference base
     * @param lastReadBaseToMarkAsIndelRelevant Final base in the read that is valid to store results for based on closeness to the edge of the indel.
     *                                          This method may compare bases beyond this point but it will not mark them as being relevant in the output
     *                                          unless the read base being compared lies before this index.
     * @param secondaryReadBreakPosition Break position to compare based on the read.Length() ending position compared to the readBases.length
     * @param refStart Starting base in the reference array to consider
     * @param refBases Array of reference bases to compare
     * @param baselineMMSums Array of mismatch scores for each position on the read. (NOTE this array should not be mutated
     *                       by this method as it is shared between calls to this method)
     * @param indelSize size of offset between reference and read bases to consider
     * @param insertion whether to compute offsets for an insertion (otherwise treats the offset as a deletion)
     */
    private static void traverseEndOfReadForIndelMismatches(final BitSet informativeBases, final int readStart, final byte[] readBases, final byte[] readQuals, final int lastReadBaseToMarkAsIndelRelevant,  final int secondaryReadBreakPosition, final int refStart, final byte[] refBases,  final int[] baselineMMSums, final int indelSize, final boolean insertion) {
        final int globalMismatchCostForReadAlignedToReference = baselineMMSums[0];
        int baseQualitySum = 0;

        // Compute how many bases forward we should compare taking into account reference/read overhang
        final int insertionLength = !insertion ? 0 : indelSize;
        final int deletionLength = insertion ? 0 : indelSize;

        // Based on the offsets and the indelSize we are considering, how many bases until we fall off the end of the read/reference arrays?
        final int numberOfBasesToDirectlyCompare = Math.min(readBases.length - readStart - insertionLength,
                refBases.length - refStart - deletionLength);

        for (int readOffset = numberOfBasesToDirectlyCompare + insertionLength - 1,
             refOffset = numberOfBasesToDirectlyCompare + deletionLength - 1;
             readOffset >= 0 && refOffset >= 0;
             readOffset--, refOffset--) {

            // Calculate the real base offset for the read:
            final byte readBase = readBases[readStart + readOffset];
            final byte refBase = refBases[refStart + refOffset];
            if (isMismatchAndNotAnAlignmentGap(readBase, refBase)) {
                baseQualitySum += readQuals[readStart + readOffset];
                if (baseQualitySum > globalMismatchCostForReadAlignedToReference) { // abort early if we are over our global mismatch cost
                    break;
                }
            }
            // The hypothetical "readOffset" that corresponds to the comparison we are currently making
            int siteOfRealComparisonPoint = Math.min(readOffset, refOffset);

            // If it's a real character and the cost isn't greater than the non-indel cost, label it as uninformative
            if (readBases[readStart + siteOfRealComparisonPoint] != AlignmentUtils.GAP_CHARACTER &&
                    // Use less than here because lastReadBaseToMarkAsIndelRelevant is the exclusive site where we flip bases later on.
                    readStart + siteOfRealComparisonPoint < lastReadBaseToMarkAsIndelRelevant &&
                    // Resolving the edge case involving read.getLength() disagreeing with the realigned indel length
                    readStart + siteOfRealComparisonPoint <= secondaryReadBreakPosition &&
                    baselineMMSums[siteOfRealComparisonPoint] >= baseQualitySum) {
                informativeBases.set(readStart + siteOfRealComparisonPoint, true); // Label with true here because we flip these results later
            }
        }
    }

    // Are these two bases different (including IUPAC bases) and does the read not correspond to a deletion on the reference
    private static boolean isMismatchAndNotAnAlignmentGap(byte readBase, byte refBase) {
        return !Nucleotide.intersect(readBase, refBase) && (readBase != AlignmentUtils.GAP_CHARACTER);
    }

    /**
     * Create a reference haplotype for an active region
     *
     * @param activeRegion the active region
     * @param refBases the ref bases
     * @param paddedReferenceLoc the location spanning of the refBases -- can be longer than activeRegion.getLocation()
     * @return a reference haplotype
     */
    public static Haplotype createReferenceHaplotype(final AssemblyRegion activeRegion, final byte[] refBases, final SimpleInterval paddedReferenceLoc) {
        Utils.nonNull(activeRegion, "null region");
        Utils.nonNull(refBases, "null refBases");
        Utils.nonNull(paddedReferenceLoc, "null paddedReferenceLoc");

        final int alignmentStart = activeRegion.getPaddedSpan().getStart() - paddedReferenceLoc.getStart();
        if ( alignmentStart < 0 ) {
            throw new IllegalStateException("Bad alignment start in createReferenceHaplotype " + alignmentStart);
        }
        final Haplotype refHaplotype = new Haplotype(refBases, true);
        refHaplotype.setGenomeLocation(activeRegion.getPaddedSpan());
        refHaplotype.setAlignmentStartHapwrtRef(alignmentStart);
        final Cigar c = new Cigar();
        c.add(new CigarElement(refHaplotype.getBases().length, CigarOperator.M));
        refHaplotype.setCigar(c);
        return refHaplotype;
    }
}
