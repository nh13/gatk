package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.Arrays;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Assigns genotypes from posteriors computed with each sample's own genotype priors, as recorded in
 * {@value GATKVCFConstants#GENOTYPE_PRIOR_PROVENANCE_KEY} (see {@link GenotypePriorProvenance}), instead of from the
 * likelihoods alone.
 *
 * <p>A sample whose priors do not apply to the site (no priors, or an allele its record did not have) keeps its
 * likelihood-based genotype. The site QUAL should be computed from posteriors only if every sample with informative
 * likelihoods got posteriors ({@link Result#allInformativeSamplesHavePosteriors()}); otherwise the posteriors of a
 * subset of samples would stand in for the whole cohort.</p>
 *
 * <p>Keeps counts of how each genotype was assigned, for {@link #logSummary}. Not thread-safe.</p>
 */
public final class InputGenotypePriors {

    /**
     * The FORMAT fields that carry a sample's input priors or posteriors. {@link #assign} removes them from every
     * genotype, and recomputes GP and PG where it assigns a genotype from posteriors.
     */
    public static final Set<String> PRIOR_FORMAT_KEYS = Set.of(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY,
            GATKVCFConstants.DRAGEN_GENOTYPE_PRIOR_KEY, GATKVCFConstants.GENOTYPE_PRIOR_KEY, VCFConstants.GENOTYPE_POSTERIORS_KEY);

    /** How a genotype with informative likelihoods was assigned. */
    enum Outcome {
        /** From posteriors with the sample's priors. */
        POSTERIOR("posteriors with the sample's priors", null),
        /** From likelihoods: the genotype carries no prior provenance. */
        NO_PRIORS("likelihoods: no usable priors for the sample (e.g. a reference block, priors dropped upstream, " +
                "or priors without provenance)", null),
        /** From likelihoods: the priors were written for a record starting elsewhere, e.g. a spanning deletion. */
        DIFFERENT_START("likelihoods: priors written for a record starting elsewhere (e.g. a spanning deletion)",
                GenotypePriorProvenance.Inapplicability.DIFFERENT_START),
        /** From likelihoods: the priors were written for another ploidy. */
        DIFFERENT_PLOIDY("likelihoods: priors written for another ploidy", GenotypePriorProvenance.Inapplicability.DIFFERENT_PLOIDY),
        /** From likelihoods: the site's reference allele does not extend the record's. */
        DIFFERENT_REFERENCE("likelihoods: site reference allele does not extend the record's",
                GenotypePriorProvenance.Inapplicability.DIFFERENT_REFERENCE),
        /** From likelihoods: the site has an allele the sample's record did not have. */
        NEW_ALLELE("likelihoods: site has an allele the sample's record did not have", GenotypePriorProvenance.Inapplicability.NEW_ALLELE),
        /** From likelihoods: the site keeps the symbolic non-reference allele. */
        NON_REF_AT_SITE("likelihoods: site keeps the <NON_REF> allele", null);

        private final String description;
        private final GenotypePriorProvenance.Inapplicability inapplicability;

        Outcome(final String description, final GenotypePriorProvenance.Inapplicability inapplicability) {
            this.description = description;
            this.inapplicability = inapplicability;
        }

        static Outcome of(final GenotypePriorProvenance.Inapplicability inapplicability) {
            return Arrays.stream(values()).filter(o -> o.inapplicability == inapplicability).findFirst()
                    .orElseThrow(() -> new IllegalArgumentException("unexpected inapplicability " + inapplicability));
        }
    }

    private final Map<Outcome, Long> genotypeCounts = new EnumMap<>(Outcome.class);
    private long sitesWithPosteriorQual = 0;
    private long sitesWithLikelihoodQual = 0;
    private long sitesWherePriorsRemovedAlternateCalls = 0;

    /**
     * Re-assigns genotypes from posteriors where the sample's priors apply.
     *
     * <p>Every genotype loses the prior fields copied from its input ({@value GATKVCFConstants#GENOTYPE_PRIOR_PROVENANCE_KEY},
     * PRI, PG and GP). A genotype with informative likelihoods and applicable priors gets its genotype, GQ, GP and PG
     * from {@link GATKVariantContextUtils#makeGenotypeCallFromPriors}; the others keep their genotype and GQ.</p>
     *
     * <p>The priors can leave no sample with an alternate allele, in which case the caller drops the site unless it
     * emits non-variant sites; such sites are counted.</p>
     *
     * @param contig the site's contig, for error messages
     * @param siteStart the start of the site
     * @param siteAlleles the site's alleles, reference first, to which the genotypes' likelihoods refer
     * @param genotypes the site's genotypes, already assigned from likelihoods
     * @return the new genotypes and whether every genotype with informative likelihoods was assigned from posteriors
     */
    Result assign(final String contig, final int siteStart, final List<Allele> siteAlleles, final GenotypesContext genotypes) {
        Utils.nonEmpty(siteAlleles, "siteAlleles");
        Utils.nonNull(genotypes, "genotypes");
        final boolean siteHasNonRef = siteAlleles.contains(Allele.NON_REF_ALLELE);
        final GenotypesContext result = GenotypesContext.create(genotypes.size());
        boolean allInformativeHavePosteriors = true;
        boolean anyInformative = false;
        for (final Genotype genotype : genotypes) {
            if (!hasInformativeLikelihoods(genotype)) {
                result.add(withoutPriorAttributes(genotype));
                continue;
            }
            final GenotypeBuilder builder = new GenotypeBuilder(withoutPriorAttributes(genotype));
            anyInformative = true;
            final Outcome outcome = siteHasNonRef ? Outcome.NON_REF_AT_SITE : assignFromPriors(builder, genotype, contig, siteStart, siteAlleles);
            genotypeCounts.merge(outcome, 1L, Long::sum);
            allInformativeHavePosteriors &= outcome == Outcome.POSTERIOR;
            result.add(builder.make());
        }
        final boolean posteriorQual = anyInformative && allInformativeHavePosteriors;
        if (posteriorQual) {
            sitesWithPosteriorQual++;
        } else {
            sitesWithLikelihoodQual++;
        }
        if (hasAlternateCall(genotypes) && !hasAlternateCall(result)) {
            sitesWherePriorsRemovedAlternateCalls++;
        }
        return new Result(result, posteriorQual);
    }

    private static Outcome assignFromPriors(final GenotypeBuilder builder, final Genotype genotype, final String contig,
                                            final int siteStart, final List<Allele> siteAlleles) {
        final GenotypePriorProvenance provenance = GenotypePriorProvenance.fromGenotype(genotype, contig, siteStart);
        if (provenance == null) {
            return Outcome.NO_PRIORS;
        }
        final GenotypePriorProvenance.Mapping mapping = provenance.mapTo(siteStart, siteAlleles, genotype.getPloidy());
        if (!mapping.isApplicable()) {
            return Outcome.of(mapping.getInapplicability());
        }
        GATKVariantContextUtils.makeGenotypeCallFromPriors(genotype.getPloidy(), builder, genotype.getLikelihoods().getAsVector(),
                siteAlleles, mapping.getLog10Priors());
        return Outcome.POSTERIOR;
    }

    private static boolean hasInformativeLikelihoods(final Genotype genotype) {
        return genotype.hasLikelihoods() && genotype.getPloidy() > 0 && GATKVariantContextUtils.isInformative(genotype.getLikelihoods().getAsVector());
    }

    private static boolean hasAlternateCall(final GenotypesContext genotypes) {
        return genotypes.stream().anyMatch(g -> g.getAlleles().stream().anyMatch(a -> a.isCalled() && a.isNonReference()));
    }

    /** @return the genotype without {@link #PRIOR_FORMAT_KEYS}, or the genotype itself if it has none of them */
    private static Genotype withoutPriorAttributes(final Genotype genotype) {
        if (PRIOR_FORMAT_KEYS.stream().noneMatch(genotype::hasExtendedAttribute)) {
            return genotype;
        }
        final Map<String, Object> attributes = new HashMap<>(genotype.getExtendedAttributes());
        attributes.keySet().removeAll(PRIOR_FORMAT_KEYS);
        return new GenotypeBuilder(genotype).noAttributes().attributes(attributes).make();
    }

    /** @return the number of genotypes with informative likelihoods assigned each way so far */
    Map<Outcome, Long> getGenotypeCounts() {
        return new EnumMap<>(genotypeCounts);
    }

    /** @return the number of sites genotyped so far whose QUAL may be computed from posteriors */
    long getSitesWithPosteriorQual() { return sitesWithPosteriorQual; }

    /** @return the number of sites genotyped so far whose QUAL must stay likelihood-based */
    long getSitesWithLikelihoodQual() { return sitesWithLikelihoodQual; }

    /** @return the number of sites genotyped so far where the priors left no sample with an alternate allele */
    long getSitesWherePriorsRemovedAlternateCalls() { return sitesWherePriorsRemovedAlternateCalls; }

    /** Logs how genotypes and site QUALs were computed. */
    public void logSummary(final Logger logger) {
        Utils.nonNull(logger, "logger");
        logger.info(String.format("Input genotype priors: of the sites genotyped with alternate alleles, QUAL came from posteriors "
                        + "at %d and from likelihoods at %d (where some sample with informative likelihoods had no applicable priors)",
                sitesWithPosteriorQual, sitesWithLikelihoodQual));
        if (sitesWherePriorsRemovedAlternateCalls > 0) {
            logger.info(String.format("Input genotype priors: at %d of those sites the priors left no sample with an alternate allele, "
                    + "so they are written only with non-variant sites", sitesWherePriorsRemovedAlternateCalls));
        }
        final String counts = Arrays.stream(Outcome.values())
                .filter(genotypeCounts::containsKey)
                .map(o -> String.format("%d from %s", genotypeCounts.get(o), o.description))
                .collect(Collectors.joining("; "));
        logger.info("Input genotype priors: genotypes with informative likelihoods were assigned " + (counts.isEmpty() ? "nowhere" : counts));
    }

    /** Genotypes re-assigned by {@link #assign}, and whether the site QUAL may be computed from their posteriors. */
    public static final class Result {
        private final GenotypesContext genotypes;
        private final boolean allInformativeSamplesHavePosteriors;

        Result(final GenotypesContext genotypes, final boolean allInformativeSamplesHavePosteriors) {
            this.genotypes = genotypes;
            this.allInformativeSamplesHavePosteriors = allInformativeSamplesHavePosteriors;
        }

        public GenotypesContext getGenotypes() { return genotypes; }

        /** @return true if at least one genotype had informative likelihoods and all such genotypes got posteriors */
        public boolean allInformativeSamplesHavePosteriors() { return allInformativeSamplesHavePosteriors; }
    }
}
